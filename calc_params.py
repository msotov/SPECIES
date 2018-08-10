'''
Last modified in: 14/03/2018
'''

import os, glob, re, math, time, logging
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, Column
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from find_mass import find_mass_age as find_mass_age
import photometric_relations as pr
from calc_errors import obtain_errors
from interpol_function import interpol
from scipy.stats import sigmaclip
from atmos import atmos
import calc_broadening as cb
from astropy.stats import sigma_clip
from scipy import stats


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def multi_run_wrapper(args):
    try:
        return run_iteration(*args)
    except Exception as e:
        print e
        logging.error('Error in code, stopping the computation')
        logging.error(e)
        return (args[0], 0.0, 0.0, 0.0, 0.0, 2, 0, 0, 0.0, 0.0,\
                0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0,\
                0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0,\
                0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0,\
                0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0,\
                2, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 'no', 'no', 'no', 0.0, 0.0, 'None')

#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def nlines(starname):
    """
    Finds the number of lines used to compute the parameters.
    Input: star name
    Returns: nFeI  : number of Fe I lines
           nFeII : number of Fe II lines
    """
    linelist = np.genfromtxt('./Spectra/linelist.dat', dtype = None, skip_header = 2,\
                            usecols = (0,4), names = ('line','ion'))
    line = linelist['line']
    ion = linelist['ion']

    file_ew = np.genfromtxt('./EW/'+starname+'.ares', dtype = None,\
                            usecols = (0, 4, 5), names = ('line','ew', 'ew_err'))
    line_ew_b = file_ew['line']
    ew_b = file_ew['ew']
    ew_err_b = file_ew['ew_err']

    i = np.where((ew_b >= 10.) & (ew_b <= 150.) & (ew_err_b/ew_b <= 2.0))[0]
    line_ew_o = line_ew_b[i]
    ew_o = ew_b[i]
    ew_err_o = ew_err_b[i]

    div = ew_err_o/ew_o
    indices = np.where(div <=1.)[0]
    line_ew = line_ew_o[indices]
    ew = ew_o[indices]
    ew_err = ew_err_o[indices]


    nFeI = 0
    nFeII = 0

    for i in range(len(line_ew)):
        index = np.where(line == line_ew[i])[0]
        if len(index)>0:
            index = int(index)
            if ion[index] == 26.0:
                nFeI += 1
            else:
                nFeII += 1

    del linelist, line, ion, file_ew, line_ew_b, ew_b, ew_err_b, line_ew_o,\
        ew_o, ew_err_o, indices, line_ew, ew, ew_err

    return nFeI, nFeII


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def compute_average_abundance(starname, w = False):
    filemoog = open('./output/' + starname + '_out.test', 'r')
    flag = nfailed = 0
    ep = dif = rw = final_Fe = -99

    abund = {'ab':{'FeI':-99, 'FeII': -999},\
             'lines':{'FeI':[], 'FeII':[]}}

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

    if w:

        filename = './EW/' + starname + '.ares'
        filelines = ascii.read(filename, include_names = ('col1', 'col5', 'col6'))
        file_wave = filelines['col1']
        file_ew = filelines['col5']
        file_e_ew = filelines['col6']

        for p in names:
            wave = np.array([float(fe.split()[0]) for fe in abund['lines'][p]])
            ab = np.array([float(fe.split()[6]) for fe in abund['lines'][p]])

            w = np.array([1./file_e_ew[int(np.where(file_wave == wave[i])[0])] for i in range(len(wave))])

            if sum(w) != 0.:
                abund['ab'][p] = round(np.average(ab, weights = w), 3)
            else:
                abund['ab'][p] = np.mean(ab)

        del filename, filelines, file_wave, file_ew, file_e_ew, wave, ab, w

    dif = abund['ab']['FeI'] - abund['ab']['FeII']
    final_Fe = abund['ab']['FeI']

    filemoog.close()

    del abund, names, names_file

    return final_Fe, ep, dif, rw, nfailed


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def runMOOG(atmos_params):
    starname = atmos_params['star']
    T = atmos_params['temperature']['value']
    g = atmos_params['gravity']['value']
    m = atmos_params['metallicity']['value']
    vt = atmos_params['velocity']['value']
    interpol(starname, T, g, m, vt)
    cmd = 'MOOGSILENT > temp.log 2>&1 <<EOF\nMOOGFEB2017/ab_%s.par\n\nEOF' % starname
    os.system(cmd)
    ab, ep, dif, rw, nfailed = compute_average_abundance(starname, w = False)
    ab = ab - 7.50
    return ab, ep, dif, rw, nfailed


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def create_linelist(file_star, ab = False):
    if ab:
        file_linelist = 'Spectra/lines_ab.dat'
    else:
        file_linelist = 'Spectra/linelist.dat'

    linelist = np.genfromtxt(file_linelist, dtype = None, skip_header = 2,\
                             names = ('line', 'excit', 'loggf', 'num', 'ion'))
    line = linelist['line']
    if ab:
        line = np.array([round(l,2) for l in line])
    excit = linelist['excit']
    loggf = linelist['loggf']
    num = linelist['num']
    ion = linelist['ion']

    file_ew = np.genfromtxt('./EW/' + file_star, dtype = None, usecols = (0, 4, 5),\
                            names = ('line', 'ew', 'ew_err'))
    line_ew_b = file_ew['line']
    ew_b = file_ew['ew']
    ew_err_b = file_ew['ew_err']

    #if ab==False:
    #    file_ew = np.genfromtxt('./EW/new_ew_file.dat', dtype = None,\
    #                            names = ('line', 'ew', 'ew_err'))
    #    line_ew_b = file_ew['line']
    #    ew_b = file_ew['ew']
    #    ew_err_b = file_ew['ew_err']

    #Take only the lines that 10. <= EW <=150 and the error in the EW is lower than the EW
    ilines = np.where((ew_b >= 10.) & (ew_b <= 150.) & (ew_err_b/ew_b <= 2.0))[0]
    line_ew = line_ew_b[ilines]
    ew = ew_b[ilines]
    ew_err = ew_err_b[ilines]

    if ab == False:

        div = ew_err/ew
        indices = np.where(div <=1.)[0]
        line_ew = line_ew[indices]
        ew = ew[indices]
        ew_err = ew_err[indices]

        del indices

    output = open('MOOG_linelist/lines.' + file_star, 'w')
    output.writelines(' %s\n' % file_star)

    for i in range(len(line_ew)):
        index = np.where(line == line_ew[i])[0]
        if len(index) == 1:
            index = int(index[0])
            output.writelines('  %7.2f    %4.1f       %5.2f     %6.3f                       %7.4f\n' %\
                              (line[index], ion[index], excit[index], loggf[index], ew[i]))

    output.close()

    del linelist, line, excit, loggf, num, ion, file_ew, line_ew_b, ew_b, ew_err_b,\
        ilines, line_ew, ew, ew_err


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************



def Metallicity(atmos):
    """
    Runs MOOG with a model atmosphere in which the
    metallicity is the value that changes.
    It stops running when the abundances derived by MOOG
    are the same as the input value.

    Input: dictionary object with the atmospheric parameters.
           It can be created using the function
           'create_atmos_params' in atmos.py

    Returns: same dictionary as input, with the metallicity and
             moog output as the ones derived here, and the last
             output from MOOG.
    """

    c = 0
    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atmos.write_debug_message('Iteration in Metallicity was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atmos.new_change()
            break

        xmetal = atmos.atmos_params['metallicity']['value']
        xmetal_antes = xmetal
        ab, ep, dif, rw = atmos.atmos_params['moog']
        if abs(ab - xmetal) <= 0.02:
            atmos.new_change()
            break

        else:
            if c > 50:
                if abs(ep) <= 0.002:
                    atmos.new_change('temperature')
                else:
                    atmos.new_change()
                break

            atmos.change_metallicity(ab)
            xmetal = atmos.atmos_params['metallicity']['value']
            bound_min = atmos.atmos_params['metallicity']['boundaries'][0]
            bound_max = atmos.atmos_params['metallicity']['boundaries'][1]
            if xmetal < bound_min or xmetal > bound_max:
                atmos.write_debug_message('Not possible to compute parameters '\
                                          'for [Fe/H] < %.1f or [Fe/H] > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atmos.new_change()
                atmos.change_metallicity(xmetal_antes)
                break
            ab, ep, dif, rw, nfailed = runMOOG(atmos.atmos_params)
            atmos.moog_output([ab, ep, dif, rw], nfailed)
            atmos.write_debug_moog()
            del ab, ep, dif, rw, xmetal, xmetal_antes
            if nfailed > 0:
                atmos.new_change()
                break
            c += 1

    del c, nit_total, nfailed
    return atmos.atmos_params


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def Temperature(atmos, i):
    """
    Runs MOOG with a model atmosphere in which the
    temperature is the value that changes.
    It stops running when the correlation between Ab(FeI)
    and excitation potential is less than 0.002

    Input: - Dictionary object with the atmospheric parameters.
             It can be created using the function
             'create_atmos_params' in atmos.py
           - i: index, from 0 to 2.

    Returns: - Same dictionary as input, with the temperature and
               moog output as the ones derived here, and the last
               output from MOOG.
             - Updated value for i.
    """

    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atmos.write_debug_message('Iteration in Temperature was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atmos.new_change()
            break

        Tantes = atmos.atmos_params['temperature']['value']
        ab, ep, dif, rw = atmos.atmos_params['moog']

        if abs(ep) <= 0.002:
            if ab == atmos.atmos_params['metallicity']['value']:
                atmos.new_change()
            else:
                atmos.new_change('velocity')
            break

        else:
            atmos.change_parameter('temperature', 1, 250., 0)
            T = atmos.atmos_params['temperature']['value']
            bound_min = atmos.atmos_params['temperature']['boundaries'][0]
            bound_max = atmos.atmos_params['temperature']['boundaries'][1]
            if T > bound_max or T < bound_min:
                atmos.write_debug_message('Not possible to compute parameters '\
                                          'for T < %d or T > %d. '\
                                          'Check the boundaries of your parameter.' % \
                                          (int(bound_min), int(bound_max)))
                atmos.new_change()
                atmos.atmos_params['temperature']['value'] = Tantes #???????????????????
                if T < 3500.:
                    atmos.atmos_params['temperature']['value'] = 3500.
                break

            ab, ep, dif, rw, nfailed = runMOOG(atmos.atmos_params)
            atmos.moog_output([ab, ep, dif, rw], nfailed)

            atmos.write_debug_moog()
            atmos.atmos_params['metallicity']['ranges'][i] = ab
            del ab, ep, dif, rw, T, Tantes
            i += 1
            if i == 3: i = 0

            if nfailed > 0:
                atmos.new_change()
                break

            if atmos.check_met():
                atmos.new_change('velocity')
                break

    del nit_total, nfailed
    return atmos.atmos_params, i


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def Pressure(atmos, i):
    """
    Runs MOOG with a model atmosphere in which the
    surface gravity is the value that changes.
    It stops running when the difference between the abundances
    derived for FeI and FeII is less than 0.002.

    Input: - Dictionary object with the atmospheric parameters.
             It can be created using the function
             'create_atmos_params' in atmos.py
           - i: index, from 0 to 2.

    Returns: - Same dictionary as input, with the gravity and
               moog output as the ones derived here, and the last
               output from MOOG.
             - Updated value for i.
    """

    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atmos.write_debug_message('Iteration in Pressure was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atmos.new_change()
            break

        logg_antes = atmos.atmos_params['gravity']['value']
        ab, ep, dif, rw = atmos.atmos_params['moog']

        if abs(dif) <= 0.002:
            if ab == atmos.atmos_params['metallicity']['value']:
                atmos.new_change()
            else:
                atmos.new_change('velocity')
            break

        else:
            atmos.change_parameter('gravity', 2, 0.25, 5)
            logg = atmos.atmos_params['gravity']['value']
            bound_min = atmos.atmos_params['gravity']['boundaries'][0]
            bound_max = atmos.atmos_params['gravity']['boundaries'][1]
            if logg < bound_min or logg > bound_max:
                atmos.write_debug_message('Not possible to compute parameters '\
                                          'for log g < %.1f or log g > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atmos.new_change()
                atmos.atmos_params['gravity']['value'] = logg_antes
                break

            ab, ep, dif, rw, nfailed = runMOOG(atmos.atmos_params)
            atmos.moog_output([ab, ep, dif, rw], nfailed)
            atmos.write_debug_moog()
            atmos.atmos_params['metallicity']['ranges'][i] = ab
            del ab, ep, dif, rw, logg, logg_antes
            i += 1
            if i == 3: i = 0

            if nfailed > 0:
                atmos.new_change()
                break

            if atmos.check_met():
                atmos.new_change('velocity')
                break

    del nit_total, nfailed

    return atmos.atmos_params, i


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def Velocity(atmos, i):
    """
    Runs MOOG with a model atmosphere in which the
    microturbulence velocity is the value that changes.
    It stops running when the correlation between Ab(FeI)
    and reduced EW is less than 0.002

    Input: - Dictionary object with the atmospheric parameters.
             It can be created using the function
             'create_atmos_params' in atmos.py
           - i: index, from 0 to 2.

    Returns: - Same dictionary as input, with the microturbulence and
               moog output as the ones derived here, and the last
               output from MOOG.
             - Updated value for i.
    """

    nfailed = 0
    nit_total = 0
    v = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atmos.write_debug_message('Iteration in Velocity was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atmos.new_change()
            break

        micro_antes = atmos.atmos_params['velocity']['value']
        ab, ep, dif, rw = atmos.atmos_params['moog']
        v += 1

        if v == 50:
            atmos.new_change()
            break

        if abs(rw) <= 0.002:
            if ab == atmos.atmos_params['metallicity']['value']:
                atmos.new_change('metallicity')
            else:
                atmos.new_change('velocity')
            break

        else:
            atmos.change_parameter('velocity', 3, 0.25, 5)
            micro = atmos.atmos_params['velocity']['value']
            bound_min = atmos.atmos_params['velocity']['boundaries'][0]
            bound_max = atmos.atmos_params['velocity']['boundaries'][1]
            if micro < bound_min or micro > bound_max:
                atmos.write_debug_message('Not possible to compute parameters '\
                                          'for micro < %.1f or micro > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atmos.new_change()
                atmos.atmos_params['velocity']['value'] = micro_antes
                break

            ab, ep, dif, rw, nfailed = runMOOG(atmos.atmos_params)
            atmos.moog_output([ab, ep, dif, rw], nfailed)
            atmos.write_debug_moog()
            atmos.atmos_params['metallicity']['ranges'][i] = ab
            del ab, ep, dif, rw, micro_antes, micro
            i += 1
            if i == 3: i = 0

            if nfailed > 0:
                atmos.new_change()
                break

            if atmos.check_met():
                atmos.new_change('velocity')
                break

    del nit_total, nfailed, v

    return atmos.atmos_params, i


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_params(star, hold, init_vals, debug, log_f = None,\
                set_boundaries = False, vals_boundaries = None,\
                in_errors = False):
    """
    Computes the stellar parameters for a certain star.
    Uses the methods: Metallicity, Temperature, Pressure and Velocity.

    Input: star: name of the star.
           hold : array with the values that are not to be changed.
           When empty, all the values will be changed.
           init_vals : initial parameters.
                       When empty, default values will be used.
           debug : True/False,
                   if you want to turn on or off the debugging option.
           file_debug : File where the output from the debugging
                        will be stored.

    Returns: T: temperature of the model.
             logg: surface gravity.
             xmetal: metallicity.
             micro: microturbulence velocity.
             exception: 1 or 2. If 1, then it founds the correct parameters.
                        If 2, it encountered a problem and couldn't converge to a solution.
    """

    #########################################################################
    # Creates the atmosphere object, where all the data concerrning the
    # computation of the atmospheric parameters will be stored.
    #########################################################################

    a = atmos(star, hold, init_vals, debug, log_f, in_errors, \
              set_boundaries, vals_boundaries)

    h_array = a.show_hold()
    a.write_debug_message('hold_t=%s, hold_p=%s, hold_m=%s, hold_v=%s' % \
                          (h_array[0], h_array[1], h_array[2], h_array[3]))
    line_hold_param = 'Boundaries are: '
    names_params = ['temperature', 'gravity', 'metallicity', 'velocity']
    for name_p in names_params:
        line_hold_param += '%s = (%.2f, %.2f), ' % (name_p, \
                           a.atmos_params[name_p]['boundaries'][0],\
                           a.atmos_params[name_p]['boundaries'][1])
    a.write_debug_message(line_hold_param[:-2])
    del line_hold_param, names_params
    a.write_debug_message('Initial values are: T=%f, logg=%f, met=%f, vt=%f' % \
                          (init_vals[0], init_vals[1], init_vals[2], init_vals[3]))

    #########################################################################
    # Modifies the output from MOOG if hold metallicity is 'yes'
    #########################################################################
    xmetal_i = a.atmos_params['metallicity']['value']
    a.check_hold(xmetal_i)

    #########################################################################
    # First iteration with MOOG.
    #########################################################################

    ab, ep, dif, rw, nfailed = runMOOG(a.atmos_params)


    a.moog_output([ab, ep, dif, rw], nfailed)
    a.write_debug_message('ab=%.3f, ep=%.3f, dif=%.3f, rw=%.3f, nfailed=%d' % \
                          (ab, ep, dif, rw, nfailed))
    a.atmos_params['metallicity']['ranges'][2] = ab

    a.write_debug_message('change=%s' % (a.atmos_params['change']))
    i = 0
    while True:
        #########################################################################
        # Values that need to be reset each iteration
        #########################################################################
        a.new_iteration(xmetal_i)

        #########################################################################
        # Found the right values
        #########################################################################
        if a.check_correct_vals() == -1: break

        #########################################################################
        # If all the parameters are out of range,
        # break the calculation
        #########################################################################
        if a.check_nout() == -1: break

        #########################################################################
        # Parameters still in the permitted ranges
        #########################################################################
        change = a.atmos_params['change']

        if change == 'metallicity':
            dic_params = Metallicity(a)
            a.atmos_params = dic_params
            ab = a.atmos_params['moog'][0]
            a.atmos_params['metallicity']['ranges'] = [-999., -999., ab]
            a.atmos_params['temperature']['ranges'] = [-999., -999.]
            a.atmos_params['gravity']['ranges'] = [-999., -999.]
            a.atmos_params['velocity']['ranges'] = [-999., -999.]
            i = 0

            a.check_nfailed_it(change)

        elif change == 'temperature':
            dic_params, i = Temperature(a, i)
            a.atmos_params = dic_params

            a.check_nfailed_it(change)

        elif change == 'pressure':
            dic_params, i = Pressure(a, i)
            a.atmos_params = dic_params

            a.check_nfailed_it(change)

        else:
            dic_params, i = Velocity(a, i)
            a.atmos_params = dic_params

            a.check_nfailed_it(change)

        #########################################################################
        # If an iteration failed, change the input parameters
        # according to a normal distribution
        #########################################################################
        a.check_nfailed()

        #########################################################################
        # If the iteration has failed more than 5 times,
        # break the calculation
        #########################################################################
        if a.check_nbreak() == -1: break

        #########################################################################
        # If the parameters for an iteration are the same
        # as a previous one, save them
        #########################################################################
        a.check_params_rep()

        #########################################################################
        # If the parameters are repeated more than
        # 500 times, break the calculation
        #########################################################################
        if a.check_nrepeat() == -1: break

        #########################################################################
        # If mnore than 1 million iterations have been performed
        # and the values have not converge, stop the calculation
        #########################################################################
        if a.check_nit_total() == -1: break


    T, logg, xmetal, micro = a.values()
    exception = a.atmos_params['exception']
    a.write_debug_message('Final parameters for %s: T=%.3f, logg=%.3f, xmetal=%.3f, micro=%.3f' %\
                          (star, T, logg, xmetal, micro))

    del a, i, xmetal_i, h_array

    return T, logg, xmetal, micro, exception


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f, set_boundaries, vals_boundaries, \
                    try_with_vt_hold, check_with_photometry_T, relation_temp, photometry,\
                    inst, nFeI, T_c):

    use_vt = 'no'

    T, logg, xmetal, micro, exception = calc_params(star, hold, init_vals, \
                                                    debug, log_f, \
                                                    set_boundaries, \
                                                    vals_boundaries)

    if (exception == 2) and (try_with_vt_hold == True) and (('velocity' in hold) == False) \
            and (nFeI > 10):
        hold_c1 = hold[:]
        init_vals_c1 = init_vals[:]
        hold_c1.append('velocity')
        init_vals_c1[3] = 1.2
        if debug:
            log_f.debug('Trying again with vt = 1.2 km/s.')
        T, logg, xmetal, micro, exception = calc_params(star, \
                                                    hold_c1, init_vals_c1, \
                                                    debug, log_f, \
                                                    set_boundaries,\
                                                    vals_boundaries)
        use_vt = 'yes'

    #########################################################################
    # Computes the temperature using the relations in Casagrande et al. 2010.
    # If the difference in T is greater than 500 K, compute everything again
    # but fixing the temperature to the one from Casagrande et al. 2010.
    #########################################################################

    use_casagrande = 'no'

    if (check_with_photometry_T == True) and (('temperature' in hold) == False):
        # If Casagrande et al. or Mann et al. were used to for Tc,
        # recompute using the metallicity from SPECIES

        if (relation_temp == 'casagrande') or (relation_temp == 'mann') and (exception == 1):
            T_c2, err_T_c2, color_c2, relation2 = pr.check_relation(photometry, xmetal, exception, inst)
            if relation2 == 'casagrande':
                if T_c2 != 0.0:
                    T_c = T_c2#_n
                    if debug:
                        log_f.debug('New temperature using the Casagrande et al. 2010 formula is %f +- %f, using %s',\
                                    T_c2, err_T_c2, color_c2)

            if relation2 == 'mann':
                if T_c2 != 0.0:
                    if exception == 1:
                        T_c2_n = pr.correct_mann(T_c2, color_c2, inst, True)
                    else:
                        T_c2_n = pr.correct_mann(T_c2, color_c2, inst, False)

                    T_c = T_c2_n
                    if debug:
                        log_f.debug('New temperature using the Mann et al. 2015 formula is %f +- %f, using %s',\
                                    T_c2_n, err_T_c2, color_c2)

            err_init_vals[0] = err_T_c2

        if exception == 2:
            if T_c > 3000. and T_c < 10000.:
                use_casagrande = 'yes'

                hold_c = hold[:]
                hold_c.append('temperature')
                init_vals_c = init_vals[:]
                init_vals_c[0] = T_c
                if debug == True:
                    log_f.debug('Difference between T_c and T is larger than 200 K '\
                                'or exception is equal to two.')
                    log_f.debug('Computing the parameters again but using T_c.')

                T, logg, xmetal, micro, exception = calc_params(star, \
                                                    hold_c, init_vals_c, \
                                                    debug, log_f, \
                                                    set_boundaries,\
                                                    vals_boundaries)

                use_vt = 'no'

                if (exception == 2) and (try_with_vt_hold == True) and (('velocity' in hold_c) == False) \
                    and (nFeI > 10):
                    hold_c.append('velocity')
                    init_vals_c[3] = 1.2
                    if debug:
                        log_f.debug('Trying again with vt = 1.2 km/s.')
                        T, logg, xmetal, micro, exception = calc_params(star, \
                                                    hold_c, init_vals_c, \
                                                    debug, log_f, \
                                                    set_boundaries,\
                                                    vals_boundaries)
                        use_vt = 'yes'


            else:
                if debug:
                    log_f.debug('Temperature from photometry outside of the permitted ranges. '\
                                'Stopping the calculation.')

        else:
            if debug:
                log_f.debug('There is agreement between the obtained temperature and from photometry.')

    return T, logg, xmetal, micro, exception, use_vt, use_casagrande, err_init_vals, log_f


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************



def run_iteration(starlist,\
                  debug = True,\
                  hold_params = False,\
                  ab = True,\
                  err = True,\
                  file_parameters = None,\
                  save_coefs_error = False,\
                  make_plots_coefs_error = False,\
                  check_with_photometry_T = False,\
                  colors_from_file = False,\
                  name_file_colors = None,\
                  hold_mass = False,\
                  file_hold_mass = None,\
                  use_coords = False,\
                  file_with_coords = None,\
                  set_boundaries = False,\
                  file_with_boundaries = None,
                  mass_from_file = False,
                  try_with_vt_hold = False,
                  recompute_with_logg_p = True,
                  age_limits = None,
                  PATH_SPECTRA = './Spectra',
                  is_giant = False):
    """
    Compute all the parameters for one star.

    Input: -starlist :
                name of the star.
           -debug :
                Bool. If you want to turn on or off the debugging option.
                Default is True.
           -hold_params:
                Bool. If you want to read from a file the initial parameters and hold options.
                If False, it will read from a random file that shouldn't exist.
           -ab :
                Bool. If you want to compute the abundances for certain elements.
                Default is True.
           -err :
                Bool. If you want to compute the errors in T, logg, metallicity and vt.
                Default is True.
           -file_parameters :
                Name of the file with the initial and hold values. Used only if hold_params = True.
                Default is None.
           -save_coefs_error :
                Bool. If you want to save the coefficients of the fits made when computing the errors.
                Used only if err = True. Default if False.
           -make_plots_coefs_error :
                Bool. If you want to save plots of the correlations between the parameters,
                when computing the errors. Used only if err = True. Default is False.
           -check_with_photometry_T :
                Bool. Check the obtained temperature with the one from using the phototric relations,
                and set it to the photometric value if they don't agree. Default is False.
           -colors_from_file :
                Bool. Retrieve the photometric magnitudes from a file instead of from catalogues in Vizier.
                Default is False.
           -name_file_colors:
                Name of the file with the magnitudes. Used only if colors_from_file is True.
           -hold_mass :
                Bool. If you want to set the mass, age, and radius of the star to a certain value.
                Default is False.
           -file_hold_mass:
                Name of the file with the values of mass, age, and radius.
                Used only if hold_mass is True.
           -use_coords :
                Bool. Use the star's coordinates instead of its name to search for the
                photometric magnitudes. Default is False.
           -file_with_coords:
                Name of the file with the star's coordinates, in deg. Used only if use_coords is True.
           -set_boundaries :
                Bool. Set the boundaries of the atmospheric parameters. Default is False.
           -file_with_boundaries :
                Name of the file with the boundaries for the atmospheric parameters.
           -mass_from_file :
                Bool. Read the information for the mass, age, and radius from a previously created file.
                The file must have been created using the SPECIES, from a previous computation.
                Default is False.
           -try_with_vt_hold :
                Bool. Set the microturbulence to 1.2 km/s when the atmospheric parameters do not converge.
                Default is False.
           -recompute_with_logg_p :
                Bool. Check the logg with the one from the isochrones models, and recompute the atmospheric
                parameters if they disagree (diference between both measurements > error in logg from isochrones).
                Default is True.
    """
    #########################################################################
    # Begins the calculation for each star.
    #########################################################################

    star = starlist
    star_name = star
    inst = 'nofound'
    indices = [indx for indx, c in enumerate(star) if '_' == c]
    if len(indices) > 0:
        inst = star[indices[-1]+1:]
        if inst == 'o':
            inst = star[indices[-2]+1:]
            star_name = star[:indices[-2]]
        else:
            star_name = star[:indices[-1]]

    #########################################################################
    # Create a symbolic copy of the restframe spectra into the default spectra
    # directory if they're not the same
    #########################################################################

    if os.path.samefile(PATH_SPECTRA, './Spectra') == False:
        if os.path.isfile(os.path.join('./Spectra', star + '_res.fits')) == False:
                os.symlink(os.path.join(PATH_SPECTRA, star + '_res.fits'), \
                           os.path.join('./Spectra', star + '_res.fits'))


    #########################################################################
    # Creates the debug file, only if debug == True
    #########################################################################

    if debug:
        log_f = logging.getLogger()
        log_f.setLevel(logging.DEBUG)
        fh = logging.FileHandler(filename = 'debugging_files/%s.log' % \
                                star, mode = 'w')
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
                        fmt='%(asctime)s %(funcName)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S'
                        )
        fh.setFormatter(formatter)
        log_f.addHandler(fh)

        log_f.debug(starlist)
        if hasattr(os, 'getppid'):
            log_f.debug('parent process: %d', os.getppid())
        log_f.debug('process id: %d', os.getpid())

    else:
        log_f = None


    #########################################################################
    # Obtain photometric information
    #########################################################################

    photometry = {}
    in_file = False

    if use_coords:
        if file_with_coords == None:
            file_with_coords = 'any_name.txt'
        if os.path.isfile(file_with_coords):
            star_coords = ascii.read(file_with_coords)
            if star_name in star_coords['Starname']:
                in_file = True
                icoords = int(np.where(star_coords['Starname'] == star_name)[0])
                star_ra = star_coords['RA'][icoords]
                star_dec = star_coords['DEC'][icoords]
                photometry = pr.vizier_params(star_name, use_coords = True, \
                                           RA = star_ra, DEC = star_dec)

            else:
                h = fits.getheader('Spectra/' + star + '.fits', 0)
                try:
                    star_ra = h['RA']
                    star_dec = h['DEC']
                    in_file = True
                    photometry = pr.vizier_params(star_name, use_coords = True, \
                                               RA = star_ra, DEC = star_dec)
                except KeyError:
                    pass
                del h

            del star_coords

        else:
            h = fits.getheader('Spectra/' + star + '.fits', 0)
            try:
                star_ra = h['RA']
                star_dec = h['DEC']
                in_file = True
                photometry = pr.vizier_params(star_name, use_coords = True, \
                                           RA = star_ra, DEC = star_dec)
            except KeyError:
                pass

            del h


    if (use_coords == False) or (in_file == False):

        if debug:
            log_f.debug('Attempting to find colors for ' + star)
        photometry = pr.vizier_params(star_name)
        photo_keys = photometry.keys()
        if len(photo_keys) == 1 or (len(photo_keys) == 3 and ('RA' in photo_keys) and ('DEC' in photo_keys)):
            try:
                new_name = fits.getval('./Spectra/' + star + '.fits', 'HIP_ID', 0)
                if new_name == '':
                    new_name = fits.getval('./Spectra/' + star + '.fits', 'OBJECT', 0)
                if debug:
                    log_f.debug('Attempting to find colors for ' + new_name)
                photometry = pr.vizier_params(new_name)
            except Exception as e:
                pass



    if colors_from_file == True:
        is_file = False
        in_file = False

        if name_file_colors == None:
            name_file_colors = 'any_name.txt'

        if os.path.isfile(name_file_colors):
            vals_magnitudes = np.genfromtxt(name_file_colors, dtype = None,
                                            names = True)

            color_names = vals_magnitudes['Starname']
            do_extinction = vals_magnitudes['Extinction']

            if star_name in color_names:
                in_file = True
                imag = int(np.where(color_names == star_name)[0])
                filter_names = ['CTIO B', 'CTIO V', 'CTIO R', 'CTIO I', '2MASS J', '2MASS H',\
                                '2MASS Ks', 'HIPPARCOS BT', 'HIPPARCOS VT', 'Stromgren b', 'Stromgren y']

                for index_mag,mag in enumerate(['B', 'V', 'R', 'I', 'J', 'H', 'K', 'Bt', 'Vt', 'b', 'y']):
                    try:
                        value_mag = vals_magnitudes[mag][imag]
                    except IndexError:
                        value_mag = vals_magnitudes[mag]
                    if value_mag != 'no':
                        photometry[mag] = [float(value_mag), 0.01, 'Given by user', filter_names[index_mag], 0.0]

                try:
                    if do_extinction[imag] == 'yes':
                        if 'RA' and 'DEC' in vals_magnitudes.keys():
                            photometry['RA'] = [float(vals_magnitudes['RA'][imag]), 'Given by user']
                            photometry['DEC'] = [float(vals_magnitudes['DEC'][imag]), 'Given by user']
                            if 'Parallax' in vals_magnitudes.keys():
                                if vals_magnitudes['Parallax'][imag] != 'no':
                                    photometry['parallax'] = [float(vals_magnitudes['Parallax'][imag]), 'Given by user']

                            photometry = pr.correct_extinction(photometry)
                except IndexError:
                    if do_extinction == 'yes':
                        if 'RA' and 'DEC' in vals_magnitudes.keys():
                            photometry['RA'] = [float(vals_magnitudes['RA']), 'Given by user']
                            photometry['DEC'] = [float(vals_magnitudes['DEC']), 'Given by user']
                            if 'Parallax' in vals_magnitudes.keys():
                                if vals_magnitudes['Parallax'] != 'no':
                                    photometry['parallax'] = [float(vals_magnitudes['Parallax']), 'Given by user']

                            photometry = pr.correct_extinction(photometry)

                if debug:
                    log_f.debug('Obtaining colors from file.')

            del vals_magnitudes, color_names, do_extinction, filter_names

    else:
        if debug:
            log_f.debug('Using colors only from Vizier.')


    if debug:
        log_f.debug('Photometry information is: %s', photometry)


    #########################################################################
    # Derive initial conditions from photometry
    #########################################################################

    relation_temp = 'None'

    # Stellar class:
    sp_class = pr.stellar_class(photometry)
    if debug:
        log_f.debug('Luminosity class is %s, from photometry.', sp_class)

    # Initial metallicity:
    ini_met = pr.ini_met(photometry)
    if debug:
        log_f.debug('Initial metallicity is %f', ini_met)

    # Initial temperature:
    if sp_class == 'dwarf':
        T_c, err_T_c, color_c, relation = pr.check_relation(photometry, ini_met, 1, inst)
        if T_c != 0.0:
            if relation == 'casagrande':
                relation_temp = 'casagrande'
                if debug:
                    log_f.debug('Temperature using the Casagrande et al. 2010 formula is %f +- %f, using %s',\
                                T_c, err_T_c, color_c)
            else:
                relation_temp = 'mann'
                if debug:
                    log_f.debug('Temperature using the Mann et al. 2015 formula is %f +- %f, using %s',\
                                T_c, err_T_c, color_c)

        else:
            T_c, err_T_c = pr.mamajek(photometry, inst)
            if T_c != 0.0:
                relation_temp = 'mamajek'
                if debug:
                    log_f.debug('Temperature using the Mamajek table is %f +- %f', T_c, err_T_c)

    else:
        T_c, err_T_c = pr.gonzalez_hernandez(photometry, ini_met)
        if T_c != 0.0:
            relation_temp = 'gonzalez_hernandez'
            if debug:
                log_f.debug('Temperature using Gonzalez-Hernandez and Bonifacio 2009 is %f +- %f', T_c, err_T_c)

    valid_ini_logg = False
    if T_c != 0.0:

        # Initial logg
        ini_logg = pr.ini_logg(T_c, sp_class)
        if ini_logg < 4.8:
            valid_ini_logg = True
        ini_logg = min(ini_logg, 4.8)
        if debug:
            log_f.debug('Initial logg is %f', ini_logg)


        init_vals = [T_c, ini_logg, ini_met, 1.23]
        err_init_vals = [err_T_c, 0.1, 0.1, 0.1]

    else:
        if debug:
            log_f.debug('No valid temperature from photometry was computed. '\
                        'Check the photometric information of your star')

        init_vals = [5500., 4.36, 0.0, 1.23]
        err_init_vals = [0.0, 0.0, 0.0, 0.0]


    #########################################################################
    # Finds the hold parameters and initial values.
    # (used only if hold_params = True)
    #########################################################################

    hold_broad = 0

    if hold_params == True:
        file_params = file_parameters
        try:
            arch_params = ascii.read(file_params)
            stars_params = np.array(arch_params['Starname']).astype('str')
            hold_par = np.array(arch_params['hold'])
            T_params = np.array(arch_params['T'])
            logg_params = np.array(arch_params['logg'])
            met_params = np.array(arch_params['met'])
            micro_params = np.array(arch_params['micro'])
            vsini_params = arch_params['vsini']
            vmac_params = arch_params['vmac']
            del arch_params


            if star_name in stars_params:
                try:
                    i = int(np.where(stars_params == star_name)[0])
                except TypeError:
                    i = int(np.where(stars_params == star_name)[0][0])


                if hold_par[i] == 'no':
                    hold = []
                else:
                    hold = hold_par[i].split(',')

                if str(T_params[i]) != 'no':
                    init_vals[0] = float(str(T_params[i]).split(',')[0])
                    if ',' in str(T_params[i]):
                        err_init_vals[0] = float(str(T_params[i]).split(',')[1])
                if str(logg_params[i]) != 'no':
                    init_vals[1] = float(str(logg_params[i]).split(',')[0])
                    if ',' in str(logg_params[i]):
                        err_init_vals[1] = float(str(logg_params[i]).split(',')[1])
                if str(met_params[i]) != 'no':
                    init_vals[2] = float(str(met_params[i]).split(',')[0])
                    if ',' in str(met_params[i]):
                        err_init_vals[2] = float(str(met_params[i]).split(',')[1])
                if str(micro_params[i]) != 'no':
                    init_vals[3] = float(str(micro_params[i]).split(',')[0])
                    if ',' in str(micro_params[i]):
                        err_init_vals[3] = float(str(micro_params[i]).split(',')[1])

                if (str(vsini_params[i]) != 'no') and (str(vmac_params[i]) != 'no'):
                    vs = float(str(vsini_params[i]).split(',')[0])
                    if ',' in str(vsini_params[i]):
                        err_vs = float(str(vsini_params[i]).split(',')[1])
                    else:
                        err_vs = 0.5
                    vm = float(str(vmac_params[i]).split(',')[0])
                    if ',' in str(vmac_params[i]):
                        err_vm = float(str(vmac_params[i]).split(',')[1])
                    else:
                        err_vm = 0.5
                    hold_broad = 1


            else:
                hold = []

            del T_params, logg_params, met_params, micro_params, hold_par,\
                stars_params, vsini_params, vmac_params


        except IOError:
            hold = []

    else:
        file_params = 'random_file.txt'
        hold = []



    #########################################################################
    # Set the boundaries, in case they are given by the user
    #########################################################################

    if set_boundaries == True and file_with_boundaries != None:
        try:
            arch_bound = ascii.read(file_with_boundaries)
            stars_bound = arch_bound['Starname']
            temperature_bound = arch_bound['temperature']
            metallicity_bound = arch_bound['metallicity']
            gravity_bound = arch_bound['gravity']
            velocity_bound = arch_bound['velocity']
            del arch_bound

            if (star_name in stars_bound) == False:
                set_boundaries = False
                vals_boundaries = {}

            else:
                i = int(np.where(stars_bound == star_name)[0])
                vals_boundaries = {}
                if (temperature_bound[i] != 'no'):
                    v = temperature_bound[i].split(',')
                    vals_boundaries['temperature'] = (float(v[0]), float(v[1]))
                if metallicity_bound[i] != 'no':
                    v = metallicity_bound[i].split(',')
                    vals_boundaries['metallicity'] = (float(v[0]), float(v[1]))
                if gravity_bound[i] != 'no':
                    v = gravity_bound[i].split(',')
                    vals_boundaries['gravity'] = (float(v[0]), float(v[1]))
                if velocity_bound[i] != 'no':
                    v = velocity_bound[i].split(',')
                    vals_boundaries['velocity'] = (float(v[0]), float(v[1]))

            del stars_bound, temperature_bound, metallicity_bound,\
                gravity_bound, velocity_bound

        except IOError:
            set_boundaries = False
            vals_boundaries = {}

    else:
        set_boundaries = False
        vals_boundaries = {}

    if (set_boundaries == False):
        vals_boundaries = {}
        if (T_c != 0.0) and (check_with_photometry_T == True):
            vals_boundaries['temperature'] = (T_c - 300.0, T_c + 300.0)
            if is_giant:
                vals_boundaries['temperature'] = (T_c - 1000.0, T_c + 1000.0)

        set_boundaries = True


    #########################################################################
    # Searches for all the stars with EW values.
    #########################################################################

    path = './MOOG_linelist/'

    #########################################################################
    # hold mass part (if hold_mass == True)
    # THIS OPTION HAS BEEN DEPRECATED.
    #########################################################################

    use_hold_mass = 0
    file_mass = 'random_file.txt'
    '''

    if hold_mass == True:
        file_mass = file_hold_mass
        try:
            arch_mass_params = ascii.read(file_mass)
            stars_hold_mass = arch_mass_params['Starname']
            mass_params = arch_mass_params['mass']
            err_mass_params = arch_mass_params['err_mass']
            age_params = arch_mass_params['age']
            err_age_params = arch_mass_params['err_age']
            logg_p_params = arch_mass_params['logg_p']
            err_logg_p_params = arch_mass_params['err_logg_p']
            radius_params = arch_mass_params['radius']
            err_radius_params = arch_mass_params['err_radius']

            if (star in stars_hold_mass) == True:
                i_mass = int(np.where(stars_hold_mass == star)[0])
                mass = mass_params[i_mass]
                err_mass = err_mass_params[i_mass]
                age = age_params[i_mass]
                err_age = err_age_params[i_mass]
                s_logg = logg_p_params[i_mass]
                err_s_logg = err_logg_p_params[i_mass]
                radius = radius_params[i_mass]
                err_radius = err_radius_params[i_mass]

                use_hold_mass = 1

            del stars_hold_mass, mass_params, err_mass_params, age_params, \
                err_age_params, logg_p_params, err_logg_p_params,\
                arch_mass_params, radius_params, err_radius_params
        except IOError:
            pass
    else:
        file_mass = 'random_file.txt'
    '''


    #########################################################################
    # Creates the EW file needed by MOOG.
    #########################################################################

    create_linelist(star + '.ares')
    nFeI, nFeII = nlines(star)

    #########################################################################
    # Modify abfind.par
    #########################################################################


    cmd = 'cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/ab_%s.par' % (star)
    os.system(cmd)

    par_out = open('./MOOGFEB2017/ab_%s_2.par' % (star), 'w')
    par_out.writelines("abfind\n")
    par_out.writelines("terminal    null\n")
    par_out.writelines("atmosphere  1\n")
    par_out.writelines("molecules   1\n")
    par_out.writelines("lines       1\n")
    par_out.writelines("flux/int    0\n")
    par_out.writelines("plot        1\n")
    par_out.writelines("damping     2\n")
    par_out.writelines("units       0\n")
    par_out.writelines("standard_out  './output/%s.test'\n" % (star))
    par_out.writelines("summary_out   './output/%s_out.test'\n" % (star))
    par_out.writelines("model_in      './atm_models/%s.atm'\n" % (star))
    par_out.writelines("lines_in      './MOOG_linelist/lines.%s.ares'\n" % (star))

    par_out.close()

    cmd = 'cp ./MOOGFEB2017/ab_%s_2.par ./MOOGFEB2017/ab_%s.par' % (star, star)
    os.system(cmd)


    #########################################################################
    # Computes the temperature, logg, metallicity and microturbulence
    #########################################################################

    T, logg, xmetal, micro, exception, use_vt, use_casagrande, err_init_vals,\
        log_f = calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f, \
                                set_boundaries, vals_boundaries, try_with_vt_hold, \
                                check_with_photometry_T,relation_temp, photometry, inst, nFeI, T_c)


    #########################################################################
    # Compute uncertainties for the atmospheric parameters
    #########################################################################

    if err == True:

        err_vt, err_vt2, err_T, err_T2,\
        err_met, err_logg = obtain_errors(star, T, xmetal, logg, micro,\
                                      debug, log_f,\
                                      save_coefs_error,\
                                      make_plots_coefs_error,\
                                      err_init_vals, hold,\
                                      use_casagrande, use_vt)

    #########################################################################
    # Compute mass, radius, age and spectroscopic logg,
    # and recomputing the atmospheric parameters if the logg measurements
    # do not agree
    #########################################################################

    print '\t\tComputing mass...'

    use_logg_p = 'no'

    if use_hold_mass == 0:

        if err == False:
            err_T, err_logg, err_met = 100., 0.25, 0.25

        if mass_from_file:
            if os.path.isfile('./isochrones/' + star + '_samples.h5') == False:
                mass_from_file = False

        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, err_logg_iso2,\
                        radius, err_radius1, err_radius2, \
                        logL, err_logL1, err_logL2 = find_mass_age(star, T, logg, xmetal,\
                                                    err_T, err_logg, err_met, exception, photometry,\
                                                    debug, log_f, mass_from_file = mass_from_file,\
                                                    age_limits = age_limits)

        if err_logg_iso1 == 0.0:
            err_logg_iso1 = 0.5
        if err_logg_iso2 == 0.0:
            err_logg_iso2 = 0.5


        if recompute_with_logg_p and (exception == 1) and (logg_iso < 5.0):
            if np.abs(logg-logg_iso) > 0.22:
                use_logg_p = 'yes'
                if debug:
                    log_f.debug('Difference between the spectroscopic and photometric logg '\
                                'are too large, recomputing the parameters again but using logg_iso.')
                print '\t\tRecomputing with logg = logg_iso'

                hold.append('pressure')
                init_vals[1] = logg_iso
                err_init_vals[1] = max(err_logg_iso1, err_logg_iso2)

                T, logg, xmetal, micro, exception, use_vt, use_casagrande, err_init_vals,\
                    log_f = calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f, \
                                            set_boundaries, vals_boundaries, try_with_vt_hold, \
                                            check_with_photometry_T, relation_temp, photometry, inst, nFeI, T_c)

                if err == True:

                    err_vt, err_vt2, err_T, err_T2,\
                    err_met, err_logg = obtain_errors(star, T, xmetal, logg, micro,\
                                                debug, log_f,\
                                                save_coefs_error,\
                                                make_plots_coefs_error,\
                                                err_init_vals, hold,\
                                                use_casagrande, use_vt)

                if err == False:
                    err_T, err_logg, err_met = 100., 0.25, 0.25

                if err_logg == 0.0:
                    err_logg = 0.5

                mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, err_logg_iso2,\
                                radius, err_radius1, err_radius2,\
                                logL, err_logL1, err_logL2 = find_mass_age(star, T, logg, xmetal,\
                                                            err_T, err_logg, err_met, exception, photometry,\
                                                            debug, log_f, mass_from_file = mass_from_file,\
                                                            age_limits = age_limits)


    #########################################################################
    # Computes the abundances of elements
    #########################################################################


    ab_FeI, ep, dif, rw, nfailed = compute_average_abundance(star, w = True)

    ab_FeI = ab_FeI - 7.50
    ab_FeII = ab_FeI - dif


    if ab == True:
        ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
        ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
        ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI,\
        ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI,\
        exception_Fe, exception_Ti = calc_ab(star, T, logg, xmetal, micro)



    #########################################################################
    # Compute broadening and vsini
    #########################################################################

    print '\t\tComputing broadening'

    if hold_broad == 0:

        broadening, error_broadening = 0.0, 0.0

        if ab == False:
            ab_NiI = xmetal
            dev_NiI = 0.05
        if dev_NiI == 0.0:
            dev_NiI = 0.05
        if np.isnan(ab_NiI):
            ab_NiI = xmetal
        if np.isnan(dev_NiI):
            dev_NiI = 0.05

        if os.path.isfile('./EW/' + star + '_vsini.ares'):

            if debug:
                log_f.debug('Values used for finding vsini and vmac are: '\
                            '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f',\
                            star, T, xmetal, logg, micro, ab_NiI,\
                            err_T, err_logg, err_met, dev_NiI)

            if exception == 1:

                vs, err_vs, vm, err_vm = cb.calc_broadening(star, T, xmetal, logg, micro, ab_NiI,\
                                                            err_T, err_logg, err_met, dev_NiI)
                if debug:
                    log_f.debug('Finished computing vsini and vmac.')

            else:
                vs, err_vs, vm, err_vm = 0.0, 0.0, 0.0, 0.0
                if debug:
                    log_f.debug('Stopping the calculation because exception = 2.')

        else:
            if debug:
                log_f.debug('File with EW for lines used in the computation of '\
                            'vsini is not present.')

            broadening, error_broadening, \
                vs, err_vs, vm, err_vm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if debug == True:
        log_f.debug('vsini = %f +/- %f, vmac = %f +/- %f', vs, err_vs, vm, err_vm)



    if debug == True:
        log_f.debug('Finished with star %s', star)

    del photometry

    #########################################################################
    # Deleted all the files that won't be used anymore.
    #########################################################################

    cmd = 'rm -f ./MOOGFEB2017/ab_%s.par ./MOOGFEB2017/ab_%s_2.par' % (star, star)
    os.system(cmd)
    cmd = 'rm -f ./atm_models/%s.atm' % (star)
    os.system(cmd)
    os.system('cp ./output/%s.test ./output/%s_out.test '\
              './output/MOOG_output_files/' % (star, star))
    cmd = 'rm -f ./output/%s_out.test ./output/%s.test' % (star, star)
    os.system(cmd)
    cmd = 'rm -f ./MOOG_linelist/lines.%s.ares' % (star)
    os.system(cmd)

    if ab == True:
        cmd = 'rm -f ./MOOGFEB2017/ab_%s_ab.par ./MOOGFEB2017/ab_%s_ab_2.par' % (star, star)
        os.system(cmd)
        cmd = 'rm -f ./atm_models/%s_ab.atm' % (star)
        os.system(cmd)
        os.system('cp ./output/%s_ab.test ./output/%s_ab_out.test '\
                  './output/MOOG_output_files/' % (star, star))
        cmd = 'rm -f ./output/%s_ab_out.test ./output/%s_ab.test' % (star, star)
        os.system(cmd)
        cmd = 'rm -f ./MOOG_linelist/lines.%s_ab.ares' % (star)
        os.system(cmd)

    if debug == True:
        fh.close()
        log_f.removeHandler(fh)
        del log_f, fh


    #########################################################################
    # Final values to return
    #########################################################################

    if ab == False and err == True:
        ab_NaI = dev_NaI = ab_MgI = dev_MgI = ab_AlI = dev_AlI = ab_SiI = dev_SiI =\
                ab_CaI = dev_CaI = ab_TiI = dev_TiI = ab_TiII = dev_TiII =\
                ab_CrI = dev_CrI = ab_MnI = dev_MnI = ab_NiI = dev_NiI =\
                ab_CuI = dev_CuI = ab_ZnI = dev_ZnI = -99
        nNaI = nMgI = nAlI = nSiI = nCaI = nTiI = nTiII = nCrI = nMnI = nNiI =\
                nCuI = nZnI = exception_Fe = exception_Ti = 0
        exception_Fe = exception_Ti = -99

    elif ab == False and err == False:
        ab_NaI = dev_NaI = ab_MgI = dev_MgI = ab_AlI = dev_AlI = ab_SiI = dev_SiI =\
                ab_CaI = dev_CaI = ab_TiI = dev_TiI = ab_TiII = dev_TiII =\
                ab_CrI = dev_CrI = ab_MnI = dev_MnI = ab_NiI = dev_NiI =\
                ab_CuI = dev_CuI = ab_ZnI = dev_ZnI = -99
        nNaI = nMgI = nAlI = nSiI = nCaI = nTiI = nTiII = nCrI = nMnI = nNiI =\
                nCuI = nZnI = exception_Fe = exception_Ti = 0
        exception_Fe = exception_Ti = -99
        err_T = err_T2 = err_logg = err_met = err_vt = err_vt2 = -99

    elif ab == True and err == False:
        err_T = err_T2 = err_logg = err_met = err_vt = err_vt2 = -99

    print '\t\tFinished with star ' + star

    return (star, T, logg, xmetal, micro, exception, nFeI, nFeII, ab_FeI, ab_FeII,\
            ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
            ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
            ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI,\
            ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI,\
            exception_Fe, exception_Ti, err_T, err_T2, err_logg, err_met, err_vt, err_vt2,\
            vs, err_vs, vm, err_vm, mass, err_mass1, err_mass2, age, err_age1, err_age2,\
            logg_iso, err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2,\
            logL, err_logL1, err_logL2, use_casagrande, use_vt, use_logg_p, \
            T_c, err_T_c, relation_temp)




#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************



def calc_ab(star, T, logg, xmetal, micro):
    """
    Input: star : Name of the star
           T : Tempertaure of the atmosphere model.
           logg : logg of the atmosphere model.
           xmetal : Metallicity of the atmosphere model.
           micro : Microturbulence velocity of the
                 atmosphere model.
    Output: ab : abundance of each element.
            dev : deviation from the mean of the abundance
                for each element.
            n : Number of lines used to compute the
                abundance for each element.
    """

    path = './MOOG_linelist/'

    #########################################################################
    # Creates the EW file needed by MOOG.
    #########################################################################

    create_linelist(star + '_ab.ares', ab = True)

    #########################################################################
    # Modify abfind.par
    #########################################################################


    cmd = 'cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/ab_%s_ab.par' % (star)
    os.system(cmd)

    par_out = open('./MOOGFEB2017/ab_%s_ab_2.par' % (star), 'w')
    par_out.writelines("abfind\n")
    par_out.writelines("terminal    null\n")
    par_out.writelines("atmosphere  1\n")
    par_out.writelines("molecules   1\n")
    par_out.writelines("lines       1\n")
    par_out.writelines("flux/int    0\n")
    par_out.writelines("plot        1\n")
    par_out.writelines("damping     2\n")
    par_out.writelines("units       0\n")
    par_out.writelines("standard_out  './output/%s_ab.test'\n" % (star))
    par_out.writelines("summary_out   './output/%s_ab_out.test'\n" % (star))
    par_out.writelines("model_in      './atm_models/%s_ab.atm'\n" % (star))
    par_out.writelines("lines_in      './MOOG_linelist/lines.%s_ab.ares'\n" % (star))

    par_out.close()

    cmd = 'cp ./MOOGFEB2017/ab_%s_ab_2.par ./MOOGFEB2017/ab_%s_ab.par' % (star, star)
    os.system(cmd)

    #########################################################################
    # Run MOOG with the parameters needed for the model atmosphere
    #########################################################################

    ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
    ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
    ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI,\
    ab_FeI, dev_FeI, nFeI, ab_FeII, dev_FeII, nFeII, ab_NiI, dev_NiI, nNiI,\
    ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI = runMOOG_ab(star + '_ab', T, \
                                                              logg, xmetal, micro)

    exception_Ti = 1
    exception_Fe = 1
    if abs(ab_TiI - ab_TiII) > 1.:
        exception_Ti = 2
    if abs(ab_FeI - ab_FeII) > 1 or abs(xmetal - ab_FeI) > 1:
        exception_Fe = 2

    return ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
           ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
           ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI,\
           ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI,\
           exception_Fe, exception_Ti



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

def runMOOG_ab(starname, T, g, m, vt, use_w = False):
    use_w = True
    interpol(starname, T, g, m, vt)
    cmd = 'MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/ab_%s.par\n\n\n\n\n\n\n\n\n\n\n\n\n\nEOF' % starname
    os.system(cmd)
    output = open('./output/' + starname + '_out.test')

    abund = {}

    lines = {}
    lines['NaI'] = lines['MgI'] = lines['AlI'] = lines['SiI'] = lines['CaI'] = \
                   lines['TiI'] = lines['TiII'] = lines['CrI'] = lines['MnI'] = \
                   lines['FeI'] = lines['FeII'] = lines['NiI'] = lines['CuI'] = \
                   lines['ZnI'] = np.array([])


    wave = {}
    wave['NaI'] = wave['MgI'] = wave['AlI'] = wave['SiI'] = wave['CaI'] = \
                  wave['TiI'] = wave['TiII'] = wave['CrI'] = wave['MnI'] = \
                  wave['FeI'] = wave['FeII'] = wave['NiI'] = wave['CuI'] = \
                  wave['ZnI'] = np.array([])


    ab = {}
    ab['NaI'] = ab['MgI'] = ab['AlI'] = ab['SiI'] = ab['CaI'] = \
                ab['TiI'] = ab['CrI'] = ab['MnI'] = \
                ab['FeI'] = ab['NiI'] = ab['CuI'] = \
                ab['ZnI'] = -99
    ab['TiII'] = ab['FeII'] = -999


    dev = {}
    dev['NaI'] = dev['MgI'] = dev['AlI'] = dev['SiI'] = dev['CaI'] = \
                 dev['TiI'] = dev['CrI'] = dev['MnI'] = \
                 dev['FeI'] = dev['NiI'] = dev['CuI'] = \
                 dev['ZnI'] = -99
    dev['TiII'] = dev['FeII'] = -999


    abund['lines'] = lines
    abund['wave'] = wave
    abund['ab'] = ab
    abund['dev'] = dev


    names = ['NaI', 'MgI','AlI', 'SiI', 'CaI', 'TiI', 'TiII', 'CrI', 'MnI', \
             'FeI', 'FeII', 'NiI', 'CuI', 'ZnI']

    names_file = ['Na I', 'Mg I','Al I', 'Si I', 'Ca I', 'Ti I', 'Ti II', 'Cr I', 'Mn I', \
                  'Fe I', 'Fe II', 'Ni I', 'Cu I', 'Zn I']


    flagab = 0
    nfailed = 0
    for linea in output:
        linea = linea.strip()
        #m = re.search(r'OH NO! ANOTHER FAILED ITERATION!', linea)
        #if m:
        #  nfailed += 1
        for p in range(len(names)):
            m = re.search(r'Abundance Results for Species ' + names_file[p] + '\.*', linea)
            if m:
                flagab = p+1

        m = re.search(r'[a-z]', linea)
        if m == None:
            m = re.search(r'[\d]', linea)
            if m:
                linea = linea.split()
                wave_linea = float(linea[0])
                ab_linea = float(linea[6])
                dev_linea = float(linea[7])
                if abs(dev_linea) < 10.0:
                    abund['lines'][names[flagab-1]] = np.append(abund['lines'][names[flagab-1]], ab_linea)
                    abund['wave'][names[flagab-1]] = np.append(abund['wave'][names[flagab-1]], wave_linea)


    del names_file


    for p in names:
        c, low, upp = sigmaclip(abund['lines'][p], 1.5, 1.5)
        index = np.where(np.in1d(abund['lines'][p], c))[0]
        if index.size > 0:
            abund['lines'][p] = abund['lines'][p][index]
            abund['wave'][p] = abund['wave'][p][index]
        del c, low, upp, index



    if use_w == True:
        filename = './EW/' + starname + '.ares'
        filelines = ascii.read(filename, include_names = ('col1', 'col5', 'col6'))
        file_wave = filelines['col1']
        file_ew = filelines['col5']
        file_e_ew = filelines['col6']

        w = {}
        for p in names:
            w[p] = np.array([1./file_e_ew[int(np.where(file_wave == abund['wave'][p][i])[0])] for i in range(len(abund['lines'][p]))])

            if sum(w[p]) != 0:
                abund['ab'][p] = round(np.average(abund['lines'][p], weights = w[p]), 3)
                abund['dev'][p] = round(np.sqrt(np.average((abund['lines'][p] - np.mean(abund['lines'][p]))**2., \
                                         weights = w[p])), 3)
            else:
                abund['ab'][p] = np.mean(abund['lines'][p])
                abund['dev'][p] = np.std(abund['lines'][p])

        del w, filelines, file_wave, file_ew, file_e_ew



    else:

        for p in names:
            if len(abund['lines'][p]) != 0:
                abund['ab'][p] = np.mean(abund['lines'][p])
                abund['dev'][p] = np.std(abund['lines'][p])


    nNaI = len(abund['lines']['NaI'])
    nMgI = len(abund['lines']['MgI'])
    nAlI = len(abund['lines']['AlI'])
    nSiI = len(abund['lines']['SiI'])
    nCaI = len(abund['lines']['CaI'])
    nTiI = len(abund['lines']['TiI'])
    nTiII = len(abund['lines']['TiII'])
    nCrI = len(abund['lines']['CrI'])
    nMnI = len(abund['lines']['MnI'])
    nFeI = len(abund['lines']['FeI'])
    nFeII = len(abund['lines']['FeII'])
    nNiI = len(abund['lines']['NiI'])
    nCuI = len(abund['lines']['CuI'])
    nZnI = len(abund['lines']['ZnI'])

    ab_NaI = abund['ab']['NaI'] - 6.24
    ab_MgI = abund['ab']['MgI'] - 7.60
    ab_AlI = abund['ab']['AlI'] - 6.45
    ab_SiI = abund['ab']['SiI'] - 7.51
    ab_CaI = abund['ab']['CaI']- 6.34
    ab_TiI = abund['ab']['TiI']- 4.95
    ab_TiII = abund['ab']['TiII']- 4.95
    ab_CrI = abund['ab']['CrI'] - 5.64
    ab_MnI = abund['ab']['MnI'] - 5.43
    ab_FeI = abund['ab']['FeI'] - 7.50
    ab_FeII = abund['ab']['FeII'] - 7.50
    ab_NiI = abund['ab']['NiI'] - 6.22
    ab_CuI = abund['ab']['CuI'] - 4.19
    ab_ZnI = abund['ab']['ZnI'] - 4.56

    dev_NaI = abund['dev']['NaI']
    dev_MgI = abund['dev']['MgI']
    dev_AlI = abund['dev']['AlI']
    dev_SiI = abund['dev']['SiI']
    dev_CaI = abund['dev']['CaI']
    dev_TiI = abund['dev']['TiI']
    dev_TiII = abund['dev']['TiII']
    dev_CrI = abund['dev']['CrI']
    dev_MnI = abund['dev']['MnI']
    dev_FeI = abund['dev']['FeI']
    dev_FeII = abund['dev']['FeII']
    dev_NiI = abund['dev']['NiI']
    dev_CuI = abund['dev']['CuI']
    dev_ZnI = abund['dev']['ZnI']


    output.close()

    del abund, names

    return ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
           ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
           ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI, \
           ab_FeI, dev_FeI, nFeI, ab_FeII, dev_FeII, nFeII, ab_NiI, dev_NiI, nNiI, \
           ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI
