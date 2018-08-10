import numpy as np
import re, math, os
from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit
from astropy.io import ascii
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip

def func_gauss(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2./(2*sigma**2.))


def compute_average_abundance_error(starname, w = False):
    filemoog = open('./output/' + starname + '_out.test', 'r')

    flagfe = 0
    ab_FeI = -99
    ep = -99
    dif = -99
    rw = -99
    ab_FeII = -999
    nflailed = 0
    linesfe1 = []
    linesfe2 = []
    ab = -99
    ab2 = -999
    final_Fe = -99
    nfailed = 0

    for line in filemoog:
        line = line.strip()
        m = re.search(r'OH NO! ANOTHER FAILED ITERATION!', line)
        if m:
            nfailed += 1
        m = re.search(r'CANNOT DECIDE ON A LINE WAVELENGTH STEP FOR', line)
        if m:
            nfailed += 1
        m = re.search(r'Abundance Results for Species (Fe I\s)\.*', line)
        if m:
            flagfe = 1
        m = re.search(r'Abundance Results for Species (Fe II)\.*', line)
        if m:
            flagfe = 2

        m = re.search(r'[a-z]', line)
        if m == None:
            m = re.search(r'[\d]', line)
            if m:
                if flagfe == 1:
                    linesfe1.append(line)
                elif flagfe == 2:
                    linesfe2.append(line)

        m = re.search(r'average abundance', line)
        if m:
            if flagfe == 1:
                lines = line.split()
                ab = float(lines[3])
                del lines
            elif flagfe == 2:
                lines = line.split()
                ab2 = float(lines[3])
                del lines

    filemoog.close()

    if w == False:
        final_Fe = ab
        dif = ab - ab2

        del flagfe, ab_FeI, linesfe1, linesfe2, ab, ab2

        return final_Fe, ep, dif, rw, nfailed

    else:
        wave_feI = np.array([float(fe1.split()[0]) for fe1 in linesfe1])
        abund_feI = np.array([float(fe1.split()[6]) for fe1 in linesfe1])
        e_abund_feI = np.array([float(fe1.split()[7]) for fe1 in linesfe1])
        wave_feII = np.array([float(fe2.split()[0]) for fe2 in linesfe2])
        abund_feII = np.array([float(fe2.split()[6]) for fe2 in linesfe2])
        e_abund_feII = np.array([float(fe2.split()[7]) for fe2 in linesfe2])
        ep_feI = np.array([float(fe1.split()[2]) for fe1 in linesfe1])

        filename = './EW/' + starname + '.ares'
        filelines = ascii.read(filename, include_names = ('col1', 'col5', 'col6'))
        file_wave = filelines['col1']
        file_ew = filelines['col5']
        file_e_ew = filelines['col6']

        w_feI = np.zeros(wave_feI.size)
        w_feII = np.zeros(wave_feII.size)
        e_w_feI = np.zeros(wave_feI.size)
        e_w_feII = np.zeros(wave_feII.size)

        for i in range(len(wave_feI)):
            index = int(np.where(file_wave == wave_feI[i])[0])
            e_w_feI[i] = file_e_ew[index]
            w_feI[i] = file_ew[index]
        for i in range(len(wave_feII)):
            index = int(np.where(file_wave == wave_feII[i])[0])
            e_w_feII[i] = file_e_ew[index]
            w_feII[i] = file_ew[index]

        del wave_feII, e_abund_feII, filename, filelines, file_wave, file_ew,\
            file_e_ew, w_feII, e_w_feII, index, flagfe, ab_FeI, \
            linesfe1, linesfe2, ab, ab2

        return np.array(wave_feI), np.array(abund_feI), np.array(e_abund_feI),\
               np.array(w_feI), np.array(e_w_feI), np.array(ep_feI),\
               np.array(abund_feII)


#******************************************************************************
#******************************************************************************


def func_lineal(x, a, b):
    return x*a + b

def func_cuadratic(x, a, b, c):
    return x*x*a + b*x + c

def piecewise_linear(p, x):
    a, b, c, x0 = p
    return np.piecewise(x, [x < x0, x >= x0], [lambda x:a*x + b, \
                        lambda x:c*(x-x0) + a*x0 + b])

def piecewise_quadratic(p, x):
    a, b, c, d, x0 = p
    return np.piecewise(x, [x < x0, x >= x0], \
                        [lambda x:a*x + b, lambda x:c*(x**2-x0**2) + \
                        d*(x-x0) + a*x0 + b])


#******************************************************************************
#******************************************************************************

def error_for_vt(starname, t_moog, xmet_moog, logg_mogg, vt_moog,\
                debug = False, file_debug = None,\
                save_coefs_error = False, file_coefs_error = None,\
                make_plots_coefs_error = False):

    #########################################################################
    # Write in the debug file
    #########################################################################
    if debug == True:
        file_debug.debug('Computing error in vt for %s', starname)

    if save_coefs_error == True:
        file_coefs_error.writelines('Error for vt\n')

    #########################################################################
    # Computes the dependency of vt with the slope of the Ab vs RW plot
    #########################################################################

    wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
        ep_feI, abund_feII = compute_average_abundance_error(starname, w = True)

    if len(w_feI)==0 or len(abund_feI)==0:
        if debug == True:
            file_debug.debug('No information found in output file. '\
                                  'Please check for next iteration. '\
                                  'Using default value for err_vt = 0.1')
        if save_coefs_error == True:
            file_coefs_error.writelines('No information found in output file. '\
                                        'Please check for next iteration. '\
                                        'Using default value for err_vt = 0.1\n')
        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
            ep_feI, abund_feII
        return 0.1, 0.1

    else:
        popt, pcov = curve_fit(func_lineal, w_feI, abund_feI)

        s = popt[0]
        err_rw0 = np.sqrt(np.diag(pcov))[0]


        c0 = 0.99014168*vt_moog + 0.4654291
        c1 = 0.02960763*vt_moog**2. + 0.8103779*vt_moog + 0.30579742
        c2 = piecewise_quadratic([0.48657998, 0.46092457,\
                                  0.19139902, 0.2605826, 1.03615874], vt_moog)
        c3 = piecewise_linear([0.07113249, 0.09986789, 0.97995149, 0.62953598],\
                               vt_moog)

        ini = -0.0025
        fin = 0.003

        knots = np.array([ini, ini, ini, ini, fin, fin, fin, fin])
        coefs = np.array([c0, c1, c2, c3, 0., 0., 0., 0.])
        tck = (knots, coefs, 3)
        ycero = splev(s, tck, der = 1)

        err1 = np.abs(ycero*err_rw0)

        if save_coefs_error == True:
            file_coefs_error.writelines('The coefficients for the spline '\
                                        'approximation for vt vs rw are\n')
            file_coefs_error.writelines('c0 = %f, c1 = %f, c2 = %f, c3 = %f\n' %\
                                        (c0, c1, c2, c3))
            file_coefs_error.writelines('Obtained dev(vt)/dev(slope) = %f\n' % ycero)
            file_coefs_error.writelines('Error in the slope is %f\n' % err_rw0)
            file_coefs_error.writelines('The error 1 for vt is %f\n' % err1)

        #########################################################################
        #########################################################################

        orden_wfeI = np.argsort(w_feI)
        w_feI_s = w_feI[orden_wfeI]
        e_abund_feI_s = e_abund_feI[orden_wfeI]
        e_w_feI_s = e_w_feI[orden_wfeI]

        c = np.abs(e_abund_feI_s)/e_w_feI_s

        # Looks for the range in EW in which the slope of the c vs EW is closest to
        # being flat, and use the value of c average over that range

        pares = []
        for i in range(len(w_feI)-5):
            for j in range(i+5, len(w_feI)):
                pares.append([i,j])

        pares = np.array(pares)

        s = []
        const = []

        for pi in pares:
            x_data = w_feI_s[pi[0]:pi[1]+1]
            y_data = c[pi[0]:pi[1]+1]
            popt,_ = curve_fit(func_lineal, x_data, y_data)
            s.append(popt[0])
            const.append(popt[1])

            del x_data, y_data

        s = np.abs(np.array(s))
        const = np.array(const)

        if len(s) == 0:
            err_vt_final = err1
            if save_coefs_error == True:
                file_coefs_error.writelines('Could not compute error. '\
                                            'It will be set to have the same '\
                                            'value as the first error for vt '\
                                            'computed previously.\n')

        else:

            c_final = const[np.argmin(s)]

            err_vt_final = c_final*(np.var(e_w_feI)/np.var(w_feI))*np.abs(ycero)

        if debug == True:
            file_debug.debug('The error for vt is %f' % (err_vt_final))
        if save_coefs_error == True:
            file_coefs_error.writelines('Variance in the errors of the EW '\
                                        'of the FeI lines is %f, and the '\
                                        'variance of the EW is %f\n' % \
                                        (np.var(e_w_feI), np.var(w_feI)))
            file_coefs_error.writelines('The error for vt is %f\n' % \
                                        (err_vt_final))

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
            ep_feI, abund_feII
        del popt, pcov, err_rw0, c0, c1, c2, c3, ini, fin, knots, coefs,\
            tck, ycero, orden_wfeI, w_feI_s, e_abund_feI_s, e_w_feI_s,\
            c, pares, const

        return err1, err_vt_final


#******************************************************************************
#******************************************************************************


def error_for_T(starname, t_moog, xmet_moog, logg_moog, vt_moog, err_vt,\
                debug = False, file_debug = None, save_coefs_error = False,\
                file_coefs_error = None, make_plots_coefs_error = False):

    #########################################################################
    # Write in the debug file
    #########################################################################

    if debug == True:
        file_debug.debug('Computing error in temperature for %s' %\
                              (starname))
    if save_coefs_error == True:
        file_coefs_error.writelines('Error in T\n')


    #########################################################################
    # Error in T because of the error
    # in the microturbulence
    #########################################################################

    c0 = -8.68977910E-03*t_moog + 8.17457805E01
    c1 = -5.39616573E-02*t_moog + 5.80470240E02
    c2 = 0.8825294*t_moog + 405.78169057

    sigma_T_vt2 = ((2.*c0*vt_moog + c1)**2.)*(err_vt**2.)

    if save_coefs_error == True:
        file_coefs_error.writelines('T = %f*vt^2 + %f*vt + %f\n' % (c0, c1, c2))
    if debug == True:
        file_debug.debug('Found T = %f*vt^2 + %f*vt + %f' % (c0, c1, c2))
        file_debug.debug('Found sigma^2 = %f' % (sigma_T_vt2))

    #########################################################################
    # Error in T because of the uncertainty
    # in the slope of the Ab vs EP plot
    #########################################################################

    if debug == True:
        file_debug.debug('Error in T because of the uncertainty in the '\
                              'slope of the Ab vs EP plot.')

    c0 = 7.42285856*t_moog - 4.00838616E04
    c1 = -8.76908285E-04*t_moog**2. + 8.81138009*t_moog - 2.73931146E04
    c2 = 1.00659297*t_moog - 41.00080888

    wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
        ep_feI, abund_feII = compute_average_abundance_error(starname, w = True)

    if len(ep_feI)==0 or len(abund_feI)==0:
        if debug == True:
            file_debug.debug('No information found in output file. '\
                                  'Please check for next iteration. '\
                                  'Using default value for err_T = 100.0')
        if save_coefs_error == True:
            file_coefs_error.writelines('No information found in output file. '\
                                        'Please check for next iteration. '\
                                        'Using default value for err_T = 100.0\n')

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
            ep_feI, abund_feII, c0, c1, c2

        return 100., 100.

    else:
        popt, pcov = curve_fit(func_lineal, ep_feI, abund_feI)
        err_ep0 = np.sqrt(np.diag(pcov))[0]
        s = popt[0]

        sigma_T_ep2_1 = (2.*c0*s + c1)**2.*(err_ep0**2.)
        sigma_T_ep2_2 = (s*t_moog**2./5040.)**2.*err_ep0**2.

        if save_coefs_error == True:
            file_coefs_error.writelines('Error in the slope of the Ab vs EP '\
                                        'plot is %f\n' % err_ep0)
            file_coefs_error.writelines('T = %f*ep^2 + %f*ep + %f\n' %\
                                        (c0, c1, c2))
        if debug == True:
            file_debug.debug('Found T = %f*ep^2 + %f*ep + %f' % \
                                  (c0, c1, c2))

        sigma_T1 = sigma_T_vt2 + sigma_T_ep2_1
        sigma_T2 = sigma_T_vt2 + sigma_T_ep2_2

        err_T_final1 = math.sqrt(sigma_T1)
        err_T_final2 = math.sqrt(sigma_T2)

        if debug == True:
            file_debug.debug('Found sigma_T1^2=%f, sigma_T2^2=%f' % \
                                  (sigma_T1, sigma_T2))
        if save_coefs_error == True:
            file_coefs_error.writelines('Found sigma_T1^2=%f, sigma_T2^2=%f\n' % \
                                        (sigma_T1, sigma_T2))


        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI,\
            ep_feI, abund_feII, c0, c1, c2, popt, pcov,\
            err_ep0, s, sigma_T_ep2_1, sigma_T_ep2_2, sigma_T1, sigma_T2

        return err_T_final1, err_T_final2


#******************************************************************************
#******************************************************************************


def error_for_met(starname, t_moog, xmet_moog, logg_moog, vt_moog, \
                  err_vt, err_T, debug = False, file_debug = None, \
                  save_coefs_error = False, file_coefs_error = None, \
                  make_plots_coefs_error = False):

    if debug == True:
        file_debug.debug('Computing error in metallicity for %s' % (starname))
    if save_coefs_error == True:
        file_coefs_error.writelines('Error in metallicity\n')

    #########################################################################
    # Error in metallicity because of the error in
    # microturbulence and temperature
    #########################################################################

    if debug == True:
        file_debug.debug('Error in metallicity because of the error '\
                              'in microturbulence and temperature.')

    c0 = 2.02478847E-05*xmet_moog - 1.82105651E-05
    c1 = -9.56658802E-05*xmet_moog + 7.02697837E-04
    c2 = -0.19272574*xmet_moog - 0.04037161
    c3 = 1.48739581*xmet_moog - 4.01159216

    sigma_vt_t2 = (c0*vt_moog + c1)**2.*err_T**2. + \
                  (c0*t_moog + c2)**2.*err_vt**2.

    if save_coefs_error == True:
        file_coefs_error.writelines('Found [Fe/H]=(%f*vt + %f)*T + '\
                                    '(%f*vt + %f)\n' % (c0, c1, c2, c3))
    if debug == True:
        file_debug.debug('Found [Fe/H]=(%f*vt + %f)*T + '\
                              '(%f*vt + %f)' % (c0, c1, c2, c3))


    #########################################################################
    # Error in metallicity because of the scatter in the line abundances
    #########################################################################

    wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
        ep_feI, abund_feII = compute_average_abundance_error(starname, w = True)

    if len(abund_feI)==0:
        if debug == True:
            file_debug.debug('No information found in output file. '\
                                  'Please check for next iteration. '\
                                  'Using default value for err_met = 0.2')
        if save_coefs_error == True:
            file_coefs_error.writelines('No information found in output file. '\
                                        'Please check for next iteration. '\
                                        'Using default value for err_met = 0.2\n')

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
            ep_feI, abund_feII, c0, c1, c2, c3, sigma_vt_t2

        return 0.2

    else:
        #sigma_xmet2 = np.var(abund_feI)
        #print e_abund_feI
        #weights_ab = 1./(e_abund_feI**2.)
        #print weights_ab
        #inan = np.where((np.isnan(weights_ab) == False) & (np.isinf(weights_ab) == False))[0]
        #print np.sum(weights_ab[inan])
        #sigma_xmet2 = (np.sum(weights_ab[inan]*abund_feI[inan]**2.)/np.sum(weights_ab[inan]) - np.average(abund_feI[inan], weights = weights_ab[inan])**2.)*len(abund_feI[inan])/(len(abund_feI[inan])-1)

        #nh, binsh = np.histogram(abund_feI[inan], bins=20, weights = weights_ab[inan])
        #bins2 = np.array([binsh[i] + (binsh[i+1] - binsh[i])/2. for i in range(len(binsh)-1)])
        #popt, pcov = curve_fit(func_gauss, bins2, nh, p0 = [max(nh), np.nanmean(abund_feI[inan]), np.std(abund_feI[inan])])
        #print popt
        #sigma_xmet2 = popt[2]**2.
        #print nh

        #fig_h, ax_h = plt.subplots()
        #ax_h.hist(abund_feI[inan], bins=20, alpha=0.5)
        #ax_h.plot(bins2, nh, marker = 'o', ls = 'None')
        #ax_h.plot(bins2, func_gauss(bins2, *popt))
        #fig_h.savefig('./output/coeffs_error_files/' + starname + '_scatter_feI.pdf')
        #plt.close(fig_h)

        imask = sigma_clip(abund_feI)
        sigma_xmet2 = np.var(abund_feI[~imask.mask])

        if debug == True:
            file_debug.debug('Error in metallicity because of the '\
                                  'scatter in the line abundances.')
        if save_coefs_error == True:
            file_coefs_error.writelines('Scatter in the FeI lines is '\
                                        'sigma_FeI^2 = %f\n' % sigma_xmet2)

        error_xmet = math.sqrt(sigma_xmet2 + sigma_vt_t2)

        if debug == True:
            file_debug.debug('Found error_xmet = %f' % (error_xmet))
        if save_coefs_error == True:
            file_coefs_error.writelines('Found error_xmet = %f\n' % (error_xmet))

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
            ep_feI, abund_feII, c0, c1, c2, c3, sigma_vt_t2,\
            sigma_xmet2

        return error_xmet


#******************************************************************************
#******************************************************************************


def error_for_logg(starname, t_moog, xmet_moog, logg_moog, vt_moog, \
                  err_T, debug = False, file_debug = None, \
                  save_coefs_error = False, file_coefs_error = None, \
                  make_plots_coefs_error = False):

    if debug == True:
        file_debug.debug('Computing error in logg for %s' % (starname))
    if save_coefs_error == True:
        file_coefs_error.writelines('Error in logg\n')

    #########################################################################
    # Error in logg because of the scatter of FeII lines
    #########################################################################

    if debug == True:
        file_debug.debug('Error in logg because of the scatter of FeII lines.')

    wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
        ep_feI, abund_feII = compute_average_abundance_error(starname, w = True)

    if len(abund_feII)==0:
        if debug == True:
            file_debug.debug('No information found in output file. '\
                                  'Please check for next iteration. '\
                                  'Using default value for err_logg = 0.5')
        if save_coefs_error == True:
            file_coefs_error.writelines('No information found in output file. '\
                                        'Please check for next iteration. '\
                                        'Using default value for err_logg = 0.5\n')

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
            ep_feI, abund_feII

        return 0.5

    else:
        var_feII = np.var(abund_feII)
        feII = np.mean(abund_feII)

        c0 = 8.31781934E-05*t_moog + 2.08059358
        c1 = -5.36075667E-04*t_moog - 1.18307730E01

        sigma_logg2 = (c0**2.)*var_feII

        if save_coefs_error == True:
            file_coefs_error.writelines('Found logg = %f*AbFeII + %f\n' % (c0, c1))
            file_coefs_error.writelines('Variance of the FeII lines is '\
                                        'sigma_FeII^2 = %f\n' % var_feII)
        if debug == True:
            file_debug.debug('Found logg = %f*AbFeII + %f' % (c0, c1))
            file_debug.debug('Obtained sigma_logg^2 = %f' % (sigma_logg2))

        #########################################################################
        # Error in logg because of the error in T
        #########################################################################

        if debug == True:
            file_debug.debug('Error in logg because of the error in T.')

        c2 = 9.20326399E-10*t_moog - 5.83146385E-06
        c3 = -2.22982964E-06*t_moog + 2.08682974E-02
        c4 = 2.60909720E-06*t_moog**2. - 2.70864935E-02*t_moog + 4.43795251E01

        sigma_logg_T2 = (2.*c2*t_moog + c3)**2.*err_T**2.

        #########################################################################
        # Sum both errors
        #########################################################################

        sigma_logg = sigma_logg2 + sigma_logg_T2
        err_logg = math.sqrt(sigma_logg)

        if save_coefs_error == True:
            file_coefs_error.writelines('Found logg = %f*T^2 + %f*T + %f\n' % \
                                        (c2, c3, c4))
            file_coefs_error.writelines('Obtained err_logg = %f\n' % (err_logg))
        if debug == True:
            file_debug.debug('Found logg = %f*T^2 + %f*T + %f' % \
                                  (c2, c3, c4))
            file_debug.debug('Obtained err_logg = %f' % (err_logg))

        del wave_feI, abund_feI, e_abund_feI, w_feI, e_w_feI, \
            ep_feI, abund_feII, var_feII, feII, c0, c1, sigma_logg2,\
            c2, c3, c4, sigma_logg_T2, sigma_logg

        return err_logg


#******************************************************************************
#******************************************************************************


def obtain_errors(starname, t_moog, xmet_moog, logg_moog, vt_moog, \
                  debug = False, file_debug = None, save_coefs_error = False, \
                  make_plots_coefs_error = False, err_init_vals = None, \
                  hold = [], use_casagrande = 'no', use_vt = 'no'):


    if debug == True:
        file_debug.debug('Computation of errors.')

    if save_coefs_error == True:
        file_coefs_error = open('./output/coeffs_error_files/%s_coefs_error.dat' % starname, 'w')
        file_coefs_error.writelines('%s\tT=%f, logg=%f, [Fe/H]=%f, vt=%f\n' %\
                                    (starname, t_moog, logg_moog, xmet_moog, vt_moog))

    if make_plots_coefs_error == True:
        if not os.path.exists('./output/plots_err'):
            os.makedirs('./output/plots_err')

    print '\t\tError for vt...'
    if (use_vt == 'no') and (('velocity' in hold) == False):
        err_vt, err_vt2 = error_for_vt(starname, t_moog, xmet_moog, logg_moog, \
                                   vt_moog, debug, file_debug, \
                                   save_coefs_error, file_coefs_error, \
                                   make_plots_coefs_error)
    else:
        if ('velocity' in hold):
            err_vt, err_vt2 = err_init_vals[3], err_init_vals[3]
        else:
            err_vt, err_vt2 = 0.1, 0.1

    print '\t\tError for T...'
    if (use_casagrande == 'no') and (('temperature' in hold) == False):
        err_T, err_T2 = error_for_T(starname, t_moog, xmet_moog, logg_moog, \
                                vt_moog, err_vt, debug, file_debug, \
                                save_coefs_error, file_coefs_error, \
                                make_plots_coefs_error)
    else:
        err_T, err_T2 = err_init_vals[0], err_init_vals[0]


    print '\t\tError for metallicity...'
    if ('metallicity' in hold) == False:
        err_met = error_for_met(starname, t_moog, xmet_moog, logg_moog, \
                            vt_moog, err_vt, err_T, debug, file_debug, \
                            save_coefs_error, file_coefs_error, \
                            make_plots_coefs_error)
    else:
        err_met = err_init_vals[2]

    print '\t\tError for logg...'
    if ('pressure' in hold) == False:
        err_logg = error_for_logg(starname, t_moog, xmet_moog, logg_moog, \
                              vt_moog, err_T, debug, file_debug, \
                              save_coefs_error, file_coefs_error, \
                              make_plots_coefs_error)
    else:
        err_logg = err_init_vals[1]

    if save_coefs_error:
        file_coefs_error.close()

    return min(err_vt,10.), min(err_vt2,10.), min(err_T,300.), min(err_T2,300.),\
           min(err_met,5.), min(err_logg,5.)
