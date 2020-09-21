from __future__ import division
from builtins import range
import re
import logging
import numpy as np
from astropy.io import ascii
from uncertainties import ufloat, unumpy
import scipy.odr as ODR
from scipy.stats import linregress
from astropy.stats import sigma_clip
from scipy import stats

def compute_average_abundance(starname, w=False, alias='test'):
    filemoog = open('./output/%s_out.test' % alias, 'r')
    flag, nfailed = 0, 0
    ep = dif = rw = final_Fe = -99

    abund = {'ab':{'FeI':[], 'FeII':[]},\
             'lines':{'FeI':[], 'FeII':[]},\
             'wave':{'FeI':[], 'FeII': []},\
             'EW':{'FeI':[], 'FeII':[]},\
             'err_EW':{'FeI':[], 'FeII':[]},\
             'err_ab': {'FeI':[], 'FeII':[]}}

    names = ['FeI', 'FeII']
    names_file = ['Fe I', 'Fe II']

    for line in filemoog:
        line = line.strip()
        m = re.search(r'OH NO! ANOTHER FAILED ITERATION!', line)
        if m:
            nfailed += 1
        m = re.search(r'CANNOT DECIDE ON A LINE WAVELENGTH STEP FOR', line)
        if m:
            nfailed += 1

        for p in range(2):
            m = re.search(r'Abundance Results for Species (' + names_file[p] + '\s)\.*', line)
            if m:
                flag = p+1

        m = re.search(r'[a-z]', line)
        if m == None:
            m = re.search(r'[\d]', line)
            if m:
                abund['lines'][names[flag - 1]].append(line)

        m = re.search(r'E.P. correlation', line)
        if m and flag == 1:
            ep = float(line.split()[4])

        m = re.search(r'R.W. correlation', line)
        if m and flag == 1:
            rw = float(line.split()[4])

        m = re.search(r'average abundance', line)
        if m:
            abund['ab'][names[flag-1]] = float(line.split()[3])

    vals_moog = [abund['ab']['FeI'], abund['ab']['FeII'], ep, rw]
    
    filename = './EW/' + starname + '.txt'
    filelines = ascii.read(filename, include_names = ('col1', 'col3', 'col4', 'col5'))
    file_wave = filelines['col1']
    file_ew = filelines['col3']
    file_e_ew = np.maximum(filelines['col4'], filelines['col5'])
    
    for p in names:
        
        ab = np.array([float(fe.split()[6]) for fe in abund['lines'][p]])
        wave = np.array([float(fe.split()[0]) for fe in abund['lines'][p]])
        e_ab = np.array([float(fe.split()[7]) for fe in abund['lines'][p]])
        ew = np.array([file_ew[int(np.where(file_wave == wave[i])[0])] for i in range(len(wave))])
        e_ew = np.array([file_e_ew[int(np.where(file_wave == wave[i])[0])] for i in range(len(wave))])
        
        abund['ab'][p] = ab
        abund['wave'][p] = wave
        abund['EW'][p] = ew
        abund['err_ab'][p] = e_ab
        abund['err_EW'][p] = e_ew
        
        if p == 'FeI':
            ep_list = np.array([float(fe.split()[2]) for fe in abund['lines'][p]])
            rw_list = np.array([float(fe.split()[5]) for fe in abund['lines'][p]])
            

    filemoog.close()

    del names, names_file, filelines, file_wave, file_ew, file_e_ew, ab, wave, e_ab, ew, e_ew

    return abund, ep_list, rw_list, vals_moog
    

def obtain_errors(starname, t_moog, xmet_moog, logg_moog, vt_moog, \
                  err_init_vals=None, \
                  hold=[], use_Tc='no', use_vt='no', alias='test',
                  read_mode='linearregression'):
    
    logging.info('Computation of errors for giants.')
    abund, ep_list, rw_list, vals_moog = compute_average_abundance(starname, alias=alias)
    if ep_list.size == 0 or rw_list.size == 0:
        logging.warning('No information found in output file. '\
                        'Please check for next iteration.')
        del abund, ep_list, rw_list
        return 0.1, 0.1, 150., 150., 0.5, 0.5
        
    if read_mode == 'linearregression':
        ab = abund['ab']['FeI']
        iclip = sigma_clip(ab, maxiters=1)
        isort_ep = np.argsort(ep_list[~iclip.mask])
        ep,_,_,_,err_ep = stats.linregress(ep_list[~iclip.mask][isort_ep], ab[~iclip.mask][isort_ep])
        isort_rw = np.argsort(rw_list[~iclip.mask])
        rw,_,_,_,err_rw = stats.linregress(rw_list[~iclip.mask][isort_rw], ab[~iclip.mask][isort_rw])
        
        m = unumpy.uarray(ab[~iclip.mask], 1.0*np.ones(ab[~iclip.mask].size))
        dif = np.mean(m).n-np.median(abund['ab']['FeII'])
        err_dif = np.mean(m).s
        err_abI = np.mean(m).s# np.std(ab[~iclip.mask])
        
    elif read_mode == 'odr':
        wave = abund['wave']['FeI']
        ew = abund['EW']['FeI']
        err_ew = abund['err_EW']['FeI']
        ab = abund['ab']['FeI']
        
        isort_ep = np.argsort(ep_list)
        isort_rw = np.argsort(rw_list)
        weights = 1./(ew_err/ew)
        weights = weights/np.sum(weights)
        
        func = ODR.polynomial(1)
        mydata = ODR.Data(x=ep_list[isort_ep], y=ab[isort_ep],
                          we=weights[isort_ep])
        myodr = ODR.ODR(mydata, func)
        myoutput = myodr.run()
        ep = myoutput.beta[1]
        err_ep = myoutput.sd_beta[1]
        del func, mydata, myodr, myoutput
        
        func = ODR.polynomial(1)
        mydata = ODR.Data(x=rw_list[isort_rw], y=ab[isort_rw],
                          we=weights[isort_rw])
        myodr = ODR.ODR(mydata, func)
        myoutput = myodr.run()
        rw = myoutput.beta[1]
        err_rw = myoutput.sd_beta[1]
        del func, mydata, myodr, myoutput
        
        m = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew)))
        err_abI = np.std(ab)
        
        abII = abund['ab']['FeII']
        wave = abund['wave']['FeII']
        ew = abund['EW']['FeII']
        err_ew = abund['err_EW']['FeII']
        m2 = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew)))
        dif = m.n-m2.n
        err_dif = m.s
        
    else:
        ab = vals_moog[0]
        abII = vals_moog[1]
        ep = vals_moog[2]
        rw = vals_moog[3]
        
        err_rw = 0.03
        err_ep = 0.008
        err_abI = np.std(ab)
        err_dif = 0.1
        dif = np.mean(ab)-np.mean(abII)
        

    # Initial values
    micro_i = ufloat(vt_moog, 0.1)
    T_i = ufloat(t_moog, 100.0)
    xmetal_i = ufloat(xmet_moog, 0.1)
    logg_i = ufloat(logg_moog, 0.1)
    
    # Error for the microturbulence
    if 'velocity' in hold:
        micro_i = ufloat(vt_moog, err_init_vals[3])
    if use_vt != 'no':
        micro_i = ufloat(vt_moog, 0.1)
    if (use_vt == 'no') and ('velocity' not in hold):
    
        z1 = unumpy.uarray([-3.25, 8.8089, -7.2177], [0.318, 0.79886, 0.494])
        z2 = unumpy.uarray([1.0726, -0.0987], [0.005, 0.006])
        rw_i = ufloat(rw, err_rw)
        c1 = np.polyval(z1, vt_moog)
        c2 = np.polyval(z2, vt_moog)
        micro_i = np.polyval([c1, c2], rw_i)
        logging.info('Error in micro is %.4f', micro_i.s)
        
    # Error for the temperature
    if (use_Tc != 'no') or ('temperature' in hold):
        T_i = ufloat(t_moog, err_init_vals[1])
    else:
        
        z1 = unumpy.uarray([-2592.84, -5028.971], [222.543, 29.47])
        z2 = unumpy.uarray([1.002, -6.914], [0.006, 27.305])
        ep_i = ufloat(ep, err_ep)
        c1 = np.polyval(z1, xmet_moog)
        c2 = np.polyval(z2, t_moog)
        T_i = np.polyval([c1, c2], ep_i)
        logging.info('Error in T is %.4f', T_i.s)
        
    if 'metallicity' in hold:
        xmetal_i = ufloat(xmet_moog, err_init_vals[0])
    else:
        xmetal_i = ufloat(xmetal_i.n, err_abI)
        logging.info('Error in met is %.4f', xmetal_i.s)
    
    if 'pressure' in hold:
        logg_i = ufloat(logg_moog, err_init_vals[2])
    else:
        z1 = unumpy.uarray([-2.379, 29.588, -137.375, 282.347, -219.058],\
                           [0.516, 6.374, 29.443, 60.222, 46.022])
        z2 = unumpy.uarray([0.978, 0.068], [0.009, 0.029])
        #dif_i = ufloat(dif, err_dif)
        dif_i = ufloat(0.0, err_dif)
        c1 = np.polyval(z1, logg_moog)
        c2 = np.polyval(z2, logg_moog)
        logg_i = np.polyval([c1, c2], dif_i)
        logging.info('Error in logg is %.4f', logg_i.s)
    
    return min(micro_i.s, 1.0), min(micro_i.s, 1.0), min(T_i.s, 300.), min(T_i.s, 300.),\
           min(xmetal_i.s, 1.0), min(logg_i.s, 1.0)
    
