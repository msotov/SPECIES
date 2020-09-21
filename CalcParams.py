"""
Last modified in: 07/11/2019
"""

from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
import os
import logging
import string
import random
import sys
from past.utils import old_div
import numpy as np
from astropy.io import ascii, fits
import astropy.units as u
from PyAstronomy import pyasl
from FindIsochrones import find_mass_age
import PhotometricRelations_class as pr
#from CalcErrors import obtain_errors
from CalcErrors_new import obtain_errors
from Atmos import calc_params, compute_average_abundance, runMOOG_ab, plot_output_file
import CalcBroadening as cb

#from EWComputation_pyspeckit import EW_calc
from EWComputation import EW_calc
from CalcErrorsGiant_new import obtain_errors as obtain_errors_giant



########################################################

def multi_run_wrapper(args):
    try:
        return run_iteration(*args)
    except Exception as e:
        print('%s: ERROR IN CODE, STOPPING THE COMPUTATION' % args[0])
        print(e)
        logging.error('%s: ERROR IN CODE, STOPPING THE COMPUTATION', args[0])
        exc_type, exc_obj, exc_tb = sys.exc_info()
        logging.error('line %d: %s', exc_tb.tb_lineno, e)
        logging.error(exc_type)
        logging.error(exc_obj)
        logging.error(exc_tb)
        return (args[0], 0.0, 0.0, 0.0, 0.0, 2, 0, 0, 0.0, 0.0,\
                {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 'no', 'no', 'no', 0.0, 0.0, 'None', 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


#@profile
def run_iteration(starlist,\
                  debug=True,\
                  hold_params=False,\
                  ab=True,\
                  err=True,\
                  file_params='any_file.txt',\
                  save_coefs_error=False,\
                  make_plots_coefs_error=False,\
                  check_with_photometry_T=False,\
                  colors_from_file=False,\
                  name_file_colors='any_file.txt',\
                  use_coords=False,\
                  file_with_coords='any_file.txt',\
                  set_boundaries=False,\
                  file_with_boundaries='any_file.txt',
                  mass_from_file=False,
                  try_with_vt_hold=False,
                  recompute_with_logg_p=True,
                  age_limits=None,
                  PATH_SPECTRA='./Spectra',
                  is_giant=False,
                  PATH_ISO='./isochrones',
                  skip_mass=False,
                  PATH_EW='./EW/plots_EW',
                  skip_vsini=False,
                  only_EW=False,
                  linelist='linelist.dat',
                  EW_recomp=True,
                  refit=False,
                  plot_mass=False,
                  minimization='per_parameter',
                  read_mode='linearregression',
                  nlist=(0, 0),
                  linelist_ab='lines_ab.dat',
                  EW_recomp_ab=True,
                  makeplotEW=False):
    """
    Compute all the parameters for one star.
    """
    #########################################################################
    # Begins the calculation for each star.
    #########################################################################

    star = starlist
    star_name = star
    inst = 'nofound'
    indices = [indx for indx, c in enumerate(star) if c == '_']
    if indices:
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

    if os.path.samefile(PATH_SPECTRA, './Spectra') is False:
        if os.path.isfile(os.path.join('./Spectra', star + '_res.fits')) is False:
            os.symlink(os.path.join(PATH_SPECTRA, star + '_res.fits'), \
                       os.path.join('./Spectra', star + '_res.fits'))


    #########################################################################
    # Creates the debug file, only if debug == True
    #########################################################################

    if debug:
        log_f = logging.getLogger()
        log_f.setLevel(logging.INFO)
        fh = logging.FileHandler(filename='debugging_files/%s.log' % \
                                star, mode='w')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter(\
                        fmt='%(asctime)s %(levelname)s %(funcName)s: %(message)s',\
                        datefmt='%Y-%m-%d %H:%M:%S')
        fh.setFormatter(formatter)
        log_f.addHandler(fh)

        logging.info(starlist)
        if hasattr(os, 'getppid'):
            logging.info('parent process: %d', os.getppid())
        logging.info('process id: %d', os.getpid())

    else:
        log_f = None


    print('\t\tStarting computation for %s (%d/%d)' % (star, nlist[0], nlist[1]))

    #########################################################################
    # Creates the EW file needed by MOOG.
    #########################################################################

    hdu = fits.open('./Spectra/%s_res.fits' % star)
    data = hdu[0].data
    header = hdu[0].header
    snr = None
    if header['NAXIS'] == 1:
        wave, flux = pyasl.read1dFitsSpec('./Spectra/%s_res.fits' % star)
    else:
        wave = data[0]
        flux = data[1]
    if 'SNR' in list(header.keys()):
        snr = header['SNR']
    if 'SNR0' in list(header.keys()) and (wave.ndim > 1):
        iorders = len(wave)
        snr = np.array([header['SNR%d' % o] for o in range(iorders)])
    if snr is None:
        raise ValueError('No valid SNR in spectrum header. Stopping calculation')
    hdu.close()

    recomp = False
    if linelist != 'linelist.dat' and EW_recomp:
        print('\t\tRecomputing EW with lines from %s' % linelist)
        logging.info('Recomputing EW with lines from %s' % linelist)
        logging.info('To disable the recomputation of the EWs when '\
                   'the linelist is different than the default one, '\
                   'use option -no_EW_recomp')
        recomp = True

    if (os.path.isfile('./EW/%s.txt' % star) is False) or recomp or refit:
        EW_calc(star, wave, flux, snr=snr, path_plots=PATH_EW, linelist=linelist,\
                makeplot=makeplotEW)
    logging.info('Created EW file.')


    if only_EW:
        logging.info('only_EW key was used, stopping the computation.')
        logging.info('Finished with star %s', star)
        if debug:
            fh.close()
            log_f.removeHandler(fh)
            del log_f, fh
        print('\t\tFinished with star %s' % star)
        return (star, 0.0, 0.0, 0.0, 0.0, 2, 0, 0, 0.0, 0.0,\
                {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                0.0, 0.0, 0.0, 'no', 'no', 'no', 0.0, 0.0, 'None', 0.0,\
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


    alias = randomString()
    #alias = star_name+'_v2'
    print('\t\tAlias for %s is %s' % (star, alias))
    logging.info('Alias for %s is %s', star, alias)

    nFeI, nFeII = create_linelist(star, linelist=linelist, alias=alias)

    #########################################################################
    # Obtain photometric information
    #########################################################################
            
    photometry = {}
    in_file = False

    if use_coords:
        star_ra, star_dec, in_file = get_coords(star, star_name, file_with_coords, PATH_SPECTRA)
        if in_file:
            photometry = pr.vizier_params(star_name, use_coords=True, \
                                          RA=star_ra, DEC=star_dec)

    if (use_coords is False) or (in_file is False):

        logging.info('Attempting to find colors for %s', star)
        photometry = pr.vizier_params(star_name)
        photo_keys = list(photometry.keys())
        if len(photo_keys) == 1 or (len(photo_keys) == 3 and ('RA' in photo_keys) and \
           ('DEC' in photo_keys)):
            try:
                new_name = fits.getval(os.path.join(PATH_SPECTRA, star + '.fits'), 'HIP_ID', 0)
                if new_name == '':
                    new_name = fits.getval(os.path.join(PATH_SPECTRA, star + '.fits'), 'OBJECT', 0)
                logging.info('Attempting to find colors for %s', new_name)
                photometry = pr.vizier_params(new_name)
            except (IOError, KeyError):
                try:
                    new_name = fits.getval(os.path.join(PATH_SPECTRA, star + '.fits'), 'OBJECT', 0)
                    logging.info('Attempting to find colors for %s', new_name)
                    photometry = pr.vizier_params(new_name)
                except (IOError, KeyError):
                    pass
                    

    if colors_from_file:
        photometry, in_file = read_colors(star_name, photometry, name_file_colors)

    else:
        logging.info('Using colors only from Vizier.')


    logging.info('Photometry information is:')
    str_ = "\n%s" % (' '*35)
    logging.info(str_.join("{} : {}".format(k, v) for k, v in photometry.items()))

    Av = pr.get_Av(photometry)


    #########################################################################
    # Derive initial conditions from photometry
    #########################################################################

    relation_temp = 'None'

    # Stellar class:
    sp_class = 'dwarf'
    sp_class1 = pr.stellar_class_pm(photometry)
    sp_class2 = pr.stellar_class(photometry)

    if (sp_class1 == 'giant') or (sp_class2 == 'giant') or is_giant:
        sp_class = 'giant'
    logging.info('Luminosity class is %s, from photometry.', sp_class)

    # Initial metallicity:
    ini_met = pr.ini_met(photometry)
    logging.info('Initial metallicity is %f', ini_met)

    # Initial temperature:
    if sp_class == 'dwarf':
        T_c, err_T_c, color_c, relation = pr.check_relation(photometry, ini_met, 1, inst)
        if T_c != 0.0:
            if relation == 'casagrande':
                relation_temp = 'casagrande'
                logging.info('Temperature using the Casagrande et al. 2010 formula is '\
                                '%f +- %f, using %s',\
                                T_c, err_T_c, color_c)
            else:
                relation_temp = 'mann'
                logging.info('Temperature using the Mann et al. 2015 formula is %f +- %f, '\
                                'using %s',\
                                T_c, err_T_c, color_c)

        else:
            T_c, err_T_c = pr.mamajek(photometry, inst)
            if T_c != 0.0:
                relation_temp = 'mamajek'
                logging.info('Temperature using the Mamajek table is %f +- %f', T_c, err_T_c)

    else:
        T_c, err_T_c = pr.alonso1999(photometry, ini_met)
        if T_c != 0.0:
            relation_temp = 'Alonso1999'
            logging.info('Temperature using the Alonso et al. 1999 relations is %f +- %f', \
                             T_c, err_T_c)

    if T_c != 0.0 and check_with_photometry_T:

        # Initial logg
        ini_logg = pr.ini_logg(T_c, sp_class)
        if sp_class == 'dwarf':
            ini_logg = min(ini_logg, 4.8)
        else:
            ini_logg = min(ini_logg, 3.5)
        logging.info('Initial logg is %f', ini_logg)

        if (3.6 < ini_logg < 4.8) and (4400 < T_c < 6400) and (-0.8 < ini_met < 0.4):
            ini_micro = 6.932*1e-4*T_c - 0.348*ini_logg - 1.437
            if ini_micro <= 0.0:
                ini_micro = 0.5
        else:
            ini_micro = 1.23

        init_vals = [ini_met, T_c, ini_logg, ini_micro]
        err_init_vals = [0.1, err_T_c, 0.1, 0.1]
    elif not check_with_photometry_T:
        logging.warning('Usage of the photometric temperature is disabled.')
        if is_giant:
            init_vals = [0.0, 4500., 3.5, 1.23]
        else:
            init_vals = [0.0, 5500., 4.36, 1.23]
        err_init_vals = [0.0, 0.0, 0.0, 0.0]
    else:
        logging.warning('No valid temperature from photometry was computed. '\
                        'Check the photometric information of your star')

        if is_giant:
            init_vals = [0.0, 4500., 3.5, 1.23]
        else:
            init_vals = [0.0, 5500., 4.36, 1.23]
        err_init_vals = [0.0, 0.0, 0.0, 0.0]


    #########################################################################
    # Finds the hold parameters and initial values.
    # (used only if hold_params = True)
    #########################################################################

    hold_broad = 0
    hold = []

    if hold_params:
        if os.path.isfile(file_params):
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

                if hold_par[i] != 'no':
                    hold = hold_par[i].split(',')

                if str(T_params[i]) != 'no':
                    init_vals[1] = float(str(T_params[i]).split(',')[0])
                    if ',' in str(T_params[i]):
                        err_init_vals[1] = float(str(T_params[i]).split(',')[1])
                if str(logg_params[i]) != 'no':
                    init_vals[2] = float(str(logg_params[i]).split(',')[0])
                    if ',' in str(logg_params[i]):
                        err_init_vals[2] = float(str(logg_params[i]).split(',')[1])
                if str(met_params[i]) != 'no':
                    init_vals[0] = float(str(met_params[i]).split(',')[0])
                    if ',' in str(met_params[i]):
                        err_init_vals[0] = float(str(met_params[i]).split(',')[1])
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

            del T_params, logg_params, met_params, micro_params, hold_par,\
                stars_params, vsini_params, vmac_params

    #########################################################################
    # Set the boundaries, in case they are given by the user
    #########################################################################
    
    set_boundaries = False
    vals_boundaries = {}

    if set_boundaries and (file_with_boundaries is not None):
        if os.path.isfile(file_with_boundaries):
            arch_bound = ascii.read(file_with_boundaries)
            stars_bound = arch_bound['Starname']
            temperature_bound = arch_bound['temperature']
            metallicity_bound = arch_bound['metallicity']
            gravity_bound = arch_bound['gravity']
            velocity_bound = arch_bound['velocity']
            del arch_bound

            if (star_name in stars_bound) is False:
                set_boundaries = False
                vals_boundaries = {}

            else:
                i = int(np.where(stars_bound == star_name)[0])
                vals_boundaries = {}
                if temperature_bound[i] != 'no':
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

    if set_boundaries is False:
        vals_boundaries = {}
        if (T_c != 0.0) and (check_with_photometry_T):
            vals_boundaries['temperature'] = (T_c - 1000.0, T_c + 1000.0)
            #vals_boundaries['temperature'] = (3500., 8500.)
            if is_giant:
                vals_boundaries['temperature'] = (3500., 5500.)
                if check_with_photometry_T and (T_c > 5500 or T_c < 3500.):
                    check_with_photometry_T = False
                    init_vals = [0.0, 4500., 2.5, 1.23]
        if is_giant:
            vals_boundaries['temperature'] = (3500., 5500.)
            vals_boundaries['gravity'] = (0.5, 3.8)

        set_boundaries = True

    #########################################################################
    # hold mass part (if hold_mass == True)
    # THIS OPTION HAS BEEN DEPRECATED.
    #########################################################################

    use_hold_mass = 0

    #########################################################################
    # Computes the temperature, logg, metallicity and microturbulence
    #########################################################################
    T, logg, xmetal, micro, exception, use_vt, use_Tc, err_init_vals,\
        log_f = calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f,\
                                set_boundaries, vals_boundaries, try_with_vt_hold,\
                                check_with_photometry_T, relation_temp, photometry,\
                                inst, nFeI, nFeII,\
                                T_c, alias, minimization=minimization,\
                                read_mode=read_mode)


    #########################################################################
    # Compute uncertainties for the atmospheric parameters
    #########################################################################

    err_vt, err_vt2, err_T, err_T2,\
            err_met, err_logg = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    if err and exception == 1:
        err_vt, err_vt2, err_T, err_T2,\
            err_met, err_logg = obtain_errors(star, T, xmetal, logg, micro,\
                                              err_init_vals, hold,\
                                              use_Tc, use_vt, alias=alias, read_mode=read_mode)
        

    #########################################################################
    # Compute mass, radius, age and spectroscopic logg,
    # and recomputing the atmospheric parameters if the logg measurements
    # do not agree
    #########################################################################

    print('\t\tComputing mass...')

    use_logg_p = 'no'

    if use_hold_mass == 0:

        if err is False:
            err_T, err_logg, err_met = 100., 0.25, 0.25
        err_T = max(err_T, 0.01)
        err_logg = max(err_logg, 0.01)
        err_met = max(err_met, 0.01)

        if mass_from_file:
            if os.path.isfile(os.path.join(os.path.relpath(PATH_ISO, '.'),\
                              '%s_samples.h5' % star)) is False:
                mass_from_file = False

        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, \
                        err_logg_iso2,\
                        radius, err_radius1, err_radius2, \
                        logL, err_logL1, err_logL2, AV, \
                        eep, err_eep1, err_eep2, \
                        p_prems, p_ms, p_rgb, p_hb, p_posthb = find_mass_age(star, T, logg, xmetal,\
                                                    err_T, err_logg, err_met, exception, \
                                                    photometry, mass_from_file=mass_from_file,\
                                                    age_limits=age_limits, is_giant=is_giant,\
                                                    PATH_ISO=PATH_ISO, skip_mass=skip_mass,\
                                                    plot_mass=plot_mass, Av=Av)

        if err_logg_iso1 == 0.0:
            err_logg_iso1 = 0.5
        if err_logg_iso2 == 0.0:
            err_logg_iso2 = 0.5

        if skip_mass:
            recompute_with_logg_p = False

        if recompute_with_logg_p and (exception == 1) and (logg_iso < 5.0):
            if np.abs(logg-logg_iso) > 0.22:
                use_logg_p = 'yes'
                logging.warning('Difference between the spectroscopic and photometric logg '\
                                'are too large, recomputing the parameters again but '\
                                'using logg_iso.')
                print('\t\tRecomputing with logg = logg_iso')

                hold.append('gravity')
                init_vals[2] = logg_iso
                err_init_vals[2] = max(err_logg_iso1, err_logg_iso2)

                T, logg, xmetal, micro, exception, use_vt, use_Tc, err_init_vals,\
                    log_f = calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f,\
                                            set_boundaries, vals_boundaries, try_with_vt_hold,\
                                            check_with_photometry_T, relation_temp, photometry,\
                                            inst, nFeI, nFeII, T_c, alias,\
                                            minimization=minimization,\
                                            read_mode=read_mode)

                if err:
                    err_vt, err_vt2, err_T, err_T2,\
                        err_met, err_logg = obtain_errors(star, T, xmetal, logg, micro,\
                                                          err_init_vals, hold,\
                                                          use_Tc, use_vt, alias=alias,\
                                                          read_mode=read_mode)

                else:
                    err_T, err_logg, err_met = 100., 0.25, 0.25

                if err_logg == 0.0:
                    err_logg = 0.5

                mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1,\
                err_logg_iso2,\
                radius, err_radius1, err_radius2,\
                logL, err_logL1, err_logL2, AV,\
                eep, err_eep1, err_eep2,\
                p_prems, p_ms, p_rgb, p_hb, p_posthb = find_mass_age(star, T, logg, xmetal,\
                                                               err_T, err_logg, err_met,
                                                               exception, photometry,\
                                                               mass_from_file=mass_from_file,\
                                                               age_limits=age_limits, \
                                                               is_giant=is_giant,\
                                                               PATH_ISO=PATH_ISO,\
                                                               skip_mass=skip_mass,\
                                                               plot_mass=plot_mass, Av=Av)

    #########################################################################
    # Computes the abundances of elements
    #########################################################################

    ab_FeI, ab_FeII = 0.0, 0.0
    if exception == 1:
        ab_FeI, _, dif, _, _ = compute_average_abundance(star, w=True, alias=alias, mode=read_mode)

        ab_FeI = ab_FeI - 7.50
        ab_FeII = ab_FeI - dif

    abund_dict = {}
    if ab and exception == 1:
        recomp = False
        if linelist_ab != 'lines_ab.dat' and EW_recomp_ab:
            print('\t\tRecomputing EW with lines from %s' % linelist_ab)
            logging.info('Recomputing EW with lines from %s' % linelist_ab)
            logging.info('To disable the recomputation of the EWs when '\
                       'the linelist is different than the default one, '\
                       'use option -no_EW_recomp_ab')
            recomp = True
        print('\t\tCreating abundance EW file for %s' % star)
        if (os.path.isfile('EW/%s_ab.txt' % star) is False) or refit or recomp:
            EW_calc('%s_ab' % star, wave, flux, snr=snr, linelist=linelist_ab,\
                    makeplot=makeplotEW)
        logging.info('Created abundance EW file.')
        
        # Check number of ions and their names
        nameions = np.unique(ascii.read('Spectra/%s' % linelist_ab, comment='-')['ele'])
        
        abund_dict = calc_ab(star, T, logg, xmetal, micro, alias,
                             nions=len(nameions), nameions=nameions, linelist=linelist_ab)
                             


    #########################################################################
    # Compute vsini and vmac
    #########################################################################

    if skip_vsini:
        logging.warning('Skipping the computation of vsini and vmac')
        vs, err_vs, vm, err_vm = 0.0, 0.0, 0.0, 0.0
    print('\t\tComputing vsini and vmac')

    dev_NiI = 0.1
    ab_NiI = xmetal
    if hold_broad == 0 and not skip_vsini:
        if ab is False:
            ab_NiI = xmetal
            dev_NiI = 0.2
        elif 'NiI' in abund_dict:
            ab_NiI = abund_dict['NiI'][0]
            dev_NiI = abund_dict['NiI'][1]
        if dev_NiI == 0.0:
            dev_NiI = 0.1
        if np.isnan(ab_NiI):
            ab_NiI = xmetal
        if np.isnan(dev_NiI):
            dev_NiI = 0.2

        if (os.path.isfile('./EW/%s_vsini.txt' % star) is False) or refit:
            print('\t\tCreating vsini EW file for ' + star)
            EW_calc('%s_vsini' % star, wave, flux, snr=snr,\
                    linelist='linelist_vsini.dat', makeplot=makeplotEW)
            logging.info('Created vsini EW file.')


        if os.path.isfile('./EW/%s_vsini.txt' % star):
            if hasattr(snr, "__len__"):
                snr = np.mean(snr)
            logging.info('Values used for finding vsini and vmac are: '\
                            '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s',\
                            star, T, xmetal, logg, micro, ab_NiI,\
                            err_T, err_logg, err_met, dev_NiI, snr, alias)

            if exception == 1:
                vs, err_vs, vm, err_vm = cb.calc_broadening(star, T, xmetal, logg, micro, ab_NiI,\
                                                            err_T, err_logg, err_met, dev_NiI,\
                                                            snr, alias)
                logging.info('Finished computing vsini and vmac.')

            else:
                vs, err_vs, vm, err_vm = 0.0, 0.0, 0.0, 0.0
                logging.warning('Stopping the calculation because exception = 2.')

        else:
            logging.warning('File with EW for lines used in the computation of '\
                            'vsini is not present.')
            vs, err_vs, vm, err_vm = 0.0, 0.0, 0.0, 0.0

    logging.info('vsini = %f +/- %f, vmac = %f +/- %f', vs, err_vs, vm, err_vm)
    logging.info('Finished with star %s', star)

    del photometry

    #########################################################################
    # Deleted all the files that won't be used anymore.
    #########################################################################

    del wave, flux, snr, hdu, data, header

    cmd = 'rm -f ./MOOGFEB2017/ab_%s.par ./MOOGFEB2017/ab_%s_2.par' % (alias, alias)
    os.system(cmd)
    cmd = 'rm -f ./atm_models/%s.atm' % (alias)
    os.system(cmd)
    if os.path.isfile('./output/%s_out.test' % alias):
        os.system('cp ./output/%s_out.test '\
                  './output/MOOG_output_files/%s_out.test' % (alias, star))
        plot_output_file(star)
    cmd = 'rm -f ./output/%s_out.test ./output/%s.test' % (alias, alias)
    os.system(cmd)
    cmd = 'rm -f ./MOOG_linelist/lines.%s.txt' % (alias)
    os.system(cmd)

    if ab:
        cmd = 'rm -f ./MOOGFEB2017/ab_%s_ab.par ./MOOGFEB2017/ab_%s_ab_2.par' % (alias, alias)
        os.system(cmd)
        cmd = 'rm -f ./atm_models/%s_ab.atm' % (alias)
        os.system(cmd)
        os.system('cp ./output/%s_ab_out.test '\
                  './output/MOOG_output_files/%s_ab_out.test' % (alias, star))
        cmd = 'rm -f ./output/%s_ab_out.test ./output/%s_ab.test' % (alias, alias)
        os.system(cmd)
        cmd = 'rm -f ./MOOG_linelist/lines.%s_ab.txt' % (alias)
        os.system(cmd)

    if debug:
        fh.close()
        log_f.removeHandler(fh)
        del log_f, fh


    #########################################################################
    # Final values to return
    #########################################################################

    print('\t\tFinished with star ' + star)

    return (star, T, logg, xmetal, micro, exception, nFeI, nFeII, ab_FeI, ab_FeII,\
            abund_dict, err_T, err_T2, err_logg, err_met, err_vt, err_vt2,\
            vs, err_vs, vm, err_vm, mass, err_mass1, err_mass2, age, err_age1, err_age2,\
            logg_iso, err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2,\
            logL, err_logL1, err_logL2, use_Tc, use_vt, use_logg_p, \
            T_c, err_T_c, relation_temp, AV, eep, err_eep1, err_eep2, \
            p_prems, p_ms, p_rgb, p_hb, p_posthb)


#******************************************************************************


def create_linelist(file_star, ab=False, linelist='linelist.dat', alias='test'):
    #if ab:
    #    file_linelist = 'Spectra/lines_ab.dat'
    #else:
    #    file_linelist = 'Spectra/%s' % linelist
    file_linelist = 'Spectra/%s' % linelist

    linelist = np.genfromtxt(file_linelist, dtype=None, skip_header=2,\
                             names=('line', 'excit', 'loggf', 'num', 'ion'))
    line = linelist['line']

    excit = linelist['excit']
    loggf = linelist['loggf']
    num = linelist['num']
    ion = linelist['ion']

    file_ew = np.genfromtxt('./EW/%s.txt' % file_star, dtype=None,\
                            names=('line', 'ew', 'ew_e', 'ew_err1', 'ew_err2'))

    line_ew_b = file_ew['line']
    ew_b = file_ew['ew']
    ew_err_b = np.maximum(file_ew['ew_err1'], file_ew['ew_err2'])

    #Take only the lines that 10. <= EW <=150 and the error in the EW is lower than the EW
    ilines = np.where((ew_b >= 10.) & (ew_b <= 150.) &\
                      (old_div(ew_err_b, ew_b) <= 0.5) & (ew_err_b > 0.0))[0]
    line_ew = line_ew_b[ilines]
    ew = ew_b[ilines]
    ew_err = ew_err_b[ilines]

    if ab is False:
        div = old_div(ew_err, ew)
        indices = np.where(div <= 1.)[0]
        line_ew = line_ew[indices]
        ew = ew[indices]
        ew_err = ew_err[indices]
        del indices

    nFeI = 0
    nFeII = 0
    
    with open('MOOG_linelist/lines.%s.txt' % alias, 'w') as output:
        output.write(' %s\n' % file_star)
        for i, l in enumerate(line_ew):
            index = np.where(line == l)[0]
            if len(index) == 1:
                index = int(index[0])
                if ion[index] == 26.0:
                    nFeI += 1
                    output.write('%s%7.2f%s%4.1f%s%5.2f%s%6.3f%s%7.4f\n' %\
                                      (' '*2, line[index], ' '*4, ion[index], ' '*7, excit[index],\
                                       ' '*5, loggf[index], ' '*23, ew[i]))
                elif ion[index] == 26.1:
                    nFeII += 1
                    output.write('%s%7.2f%s%4.1f%s%5.2f%s%6.3f%s%7.4f\n' %\
                                      (' '*2, line[index], ' '*4, ion[index], ' '*7, excit[index],\
                                       ' '*5, loggf[index], ' '*23, ew[i]))
                else:
                    output.write('%s%7.2f%s%4.1f%s%5.2f%s%6.3f%s%7.4f\n' %\
                                      (' '*2, line[index], ' '*4, ion[index], ' '*7, excit[index],\
                                       ' '*5, loggf[index], ' '*23, ew[i]))
    del output

    del linelist, line, excit, loggf, num, ion, file_ew, line_ew_b, ew_b, ew_err_b,\
        ilines, line_ew, ew, ew_err

    return nFeI, nFeII


#******************************************************************************

def calc_atm_params(star, hold, init_vals, err_init_vals, debug, log_f, set_boundaries,
                    vals_boundaries,\
                    try_with_vt_hold, check_with_photometry_T, relation_temp, photometry,\
                    inst, nFeI, nFeII, T_c, alias, minimization='per_parameter',\
                    read_mode='linearregression'):

    if nFeII == 0:
        logging.error('No FeII lines were detected. Stopping the computation')
        return 0.0, 0.0, 0.0, 0.0, 2, 'no', 'no', err_init_vals, log_f

    os.system('cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/ab_%s.par' % alias)

    with open('./MOOGFEB2017/ab_%s_2.par' % alias, 'w') as par_out:
        par_out.write("abfind\n")
        par_out.write("terminal    null\n")
        par_out.write("atmosphere  1\n")
        par_out.write("molecules   1\n")
        par_out.write("lines       1\n")
        par_out.write("flux/int    0\n")
        par_out.write("plot        1\n")
        par_out.write("damping     2\n")
        par_out.write("units       0\n")
        par_out.write("standard_out  './output/%s.test'\n" % alias)
        par_out.write("summary_out   './output/%s_out.test'\n" % alias)
        par_out.write("model_in      './atm_models/%s.atm'\n" % alias)
        par_out.write("lines_in      './MOOG_linelist/lines.%s.txt'\n" % alias)

    os.system('cp ./MOOGFEB2017/ab_{0}_2.par ./MOOGFEB2017/ab_{0}.par'.format(alias))

    use_vt = 'no'

    T, logg, xmetal, micro, exception = calc_params(star, hold, init_vals, \
                                                    debug, log_f, \
                                                    set_boundaries, \
                                                    vals_boundaries, alias=alias,
                                                    minimization=minimization,
                                                    read_mode=read_mode)
                                                    
    if (exception == 2) and (minimization == 'downhill_simplex'):
        logging.warning('Trying again, but only requiring convergence just once.')
        T, logg, xmetal, micro, exception = calc_params(star, hold, init_vals, \
                                                    debug, log_f, \
                                                    set_boundaries, \
                                                    vals_boundaries, alias=alias,
                                                    minimization=minimization,
                                                    one_round=True,
                                                    read_mode=read_mode)

    if (exception == 2) and (try_with_vt_hold) and ('velocity' not in hold) and (nFeI > 10):
        hold_c1 = hold[:]
        init_vals_c1 = init_vals[:]
        hold_c1.append('velocity')
        logging.info('Trying again with vt = %.2f km/s.', init_vals_c1[3])
        T, logg, xmetal, micro, exception = calc_params(star, \
                                                        hold_c1, init_vals_c1, \
                                                        debug, log_f, \
                                                        set_boundaries,\
                                                        vals_boundaries, alias=alias,
                                                        minimization=minimization,
                                                        read_mode=read_mode)
        use_vt = 'yes'

    #########################################################################
    # Computes the temperature using the relations in Casagrande et al. 2010.
    # If the difference in T is greater than 500 K, compute everything again
    # but fixing the temperature to the one from Casagrande et al. 2010.
    #########################################################################

    use_Tc = 'no'

    if (check_with_photometry_T) and (('temperature' in hold) is False):
        # If Casagrande et al. or Mann et al. were used to for Tc,
        # recompute using the metallicity from SPECIES

        if (relation_temp == 'casagrande') or (relation_temp == 'mann') and (exception == 1):
            T_c2, err_T_c2, color_c2, relation2 = pr.check_relation(photometry, xmetal, \
                                                                    exception, inst)
            if relation2 == 'casagrande':
                if T_c2 != 0.0:
                    T_c = T_c2#_n
                    logging.info('New temperature using the Casagrande et al. 2010 '\
                               'formula is %f +- %f, using %s',\
                               T_c2, err_T_c2, color_c2)

            if relation2 == 'mann':
                if T_c2 != 0.0:
                    if exception == 1:
                        T_c2_n = pr.correct_mann(T_c2, color_c2, inst, True)
                    else:
                        T_c2_n = pr.correct_mann(T_c2, color_c2, inst, False)

                    T_c = T_c2_n
                    logging.info('New temperature using the Mann et al. 2015 formula '\
                               'is %f +- %f, using %s',\
                               T_c2_n, err_T_c2, color_c2)

            err_init_vals[1] = err_T_c2

        if exception == 2:
            if 3000. < T_c < 10000.:
                use_Tc = 'yes'

                hold_c = hold[:]
                hold_c.append('temperature')
                init_vals_c = init_vals[:]
                init_vals_c[1] = T_c
                logging.warning('Difference between T_c and T is larger than 200 K '\
                           'or exception is equal to two.')
                logging.info('Computing the parameters again but using T_c.')

                T, logg, xmetal, micro, exception = calc_params(star, \
                                                                hold_c, init_vals_c, \
                                                                debug, log_f, \
                                                                set_boundaries,\
                                                                vals_boundaries, alias=alias,
                                                                minimization=minimization,
                                                                read_mode=read_mode)

                use_vt = 'no'

                if (exception == 2) and (try_with_vt_hold) and \
                   (('velocity' in hold_c) is False) and (nFeI > 10):
                    hold_c.append('velocity')
                    logging.info('Trying again with vt = %.2f km/s.', init_vals_c[3])
                    T, logg, xmetal, micro, exception = calc_params(star, \
                                                                    hold_c, init_vals_c, \
                                                                    debug, log_f, \
                                                                    set_boundaries,\
                                                                    vals_boundaries, alias=alias,
                                                                    minimization=minimization,
                                                                    read_mode=read_mode)
                    use_vt = 'yes'

            else:
                logging.warning('Temperature from photometry outside of the permitted ranges. '\
                           'Stopping the calculation.')

        else:
            logging.info('There is agreement between the obtained temperature'\
                       ' and from photometry.')

    return T, logg, xmetal, micro, exception, use_vt, use_Tc, err_init_vals, log_f


#******************************************************************************


def randomString(stringLength=6):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))



def get_coords(star, star_name, file_with_coords, PATH_SPECTRA):
    in_file = False
    def get_from_header():
        h = fits.getheader(os.path.join(PATH_SPECTRA, '%s.fits' % star), 0)
        if 'RA' in h and 'DEC' in h:
            star_ra = h['RA']
            star_dec = h['DEC']
            in_file = True
            del h
            return star_ra, star_dec, in_file
        del h
        return None, None, False

    if file_with_coords is None:
        star_ra, star_dec, in_file = get_from_header()
    else:
        if os.path.isfile(file_with_coords):
            star_coords = ascii.read(file_with_coords)
            if star_name in star_coords['Starname']:
                in_file = True
                icoords = int(np.where(star_coords['Starname'] == star_name)[0])
                star_ra = star_coords['RA'][icoords]
                star_dec = star_coords['DEC'][icoords]
            else:
                star_ra, star_dec, in_file = get_from_header()
        else:
            star_ra, star_dec, in_file = get_from_header()
    return star_ra, star_dec, in_file
    
    
def read_colors(star_name, photometry, name_file_colors):
    in_file = False
    if os.path.isfile(str(name_file_colors)):
        vals_magnitudes = ascii.read(name_file_colors, delimiter=' ')
        color_names = vals_magnitudes['Starname']
        if star_name in color_names:
            in_file = True
            logging.info('Obtaining colors from file.')
            imag = int(np.where(color_names == star_name)[0])
            filter_names = ['CTIO B', 'CTIO V', 'CTIO R', 'CTIO I', '2MASS J', '2MASS H',\
                            '2MASS Ks', 'HIPPARCOS BT', 'HIPPARCOS VT', \
                            'Stromgren b', 'Stromgren y', 'Gaia']
            magnitudes = ['B', 'V', 'R', 'I', 'J', 'H', 'K', 'Bt', 'Vt', 'b', 'y', 'G']
                                
            for index_mag, mag in enumerate(magnitudes):
                try:
                    value_mag = vals_magnitudes[mag][imag]
                except IndexError:
                    value_mag = vals_magnitudes[mag]
                if value_mag != 'no':
                    photometry[mag] = [float(value_mag), 0.01, 'Given by user', \
                                       filter_names[index_mag], 0.0]
                
            for mag in ['pmRA', 'pmDEC', 'RA', 'DEC', 'Parallax']:
                if mag in vals_magnitudes.colnames:
                    try:
                        value_mag = vals_magnitudes[mag][imag]
                    except IndexError:
                        value_mag = vals_magnitudes[mag]
                    if value_mag != 'no':
                        if isinstance(value_mag, str) and (',' in value_mag):
                            photometry[mag] = [float(value_mag.split(',')[0]), \
                                               float(value_mag.split(',')[1]), 'Given by user']
                        else:
                            photometry[mag] = [float(value_mag), 0.1, 'Given by user']
                        if mag in ['RA', 'DEC']:
                            photometry[mag] = [float(value_mag), 'Given by user']
                        elif mag == 'Parallax':
                            photometry['parallax'] = [float(value_mag.split(',')[0])*u.mas,\
                                               float(value_mag.split(',')[1])*u.mas,\
                                               'Given by user']
            do_extinction = vals_magnitudes['Extinction']
            try:
                extinction = do_extinction[imag]
            except IndexError:
                extinction = do_extinction
            if extinction == 'yes':
                photometry = pr.correct_extinction(photometry)
            del do_extinction
        del vals_magnitudes, color_names
    return photometry, in_file


#******************************************************************************



def calc_ab(star, T, logg, xmetal, micro, alias='test', nions=14,\
            nameions=['NaI', 'MgI', 'AlI', 'SiI', 'CaI', 'TiI', 'TiII', 'CrI', 'MnI', \
                      'FeI', 'FeII', 'NiI', 'CuI', 'ZnI'], linelist='lines_ab.dat'):

    _, _ = create_linelist(star + '_ab', ab=True, alias=alias + '_ab', linelist=linelist)


    os.system('cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/ab_%s_ab.par' % alias)

    par_out = open('./MOOGFEB2017/ab_%s_ab_2.par' % alias, 'w')
    par_out.write("abfind\n")
    par_out.write("terminal    null\n")
    par_out.write("atmosphere  1\n")
    par_out.write("molecules   1\n")
    par_out.write("lines       1\n")
    par_out.write("flux/int    0\n")
    par_out.write("plot        1\n")
    par_out.write("damping     2\n")
    par_out.write("units       0\n")
    par_out.write("standard_out  './output/%s_ab.test'\n" % (alias))
    par_out.write("summary_out   './output/%s_ab_out.test'\n" % (alias))
    par_out.write("model_in      './atm_models/%s_ab.atm'\n" % (alias))
    par_out.write("lines_in      './MOOG_linelist/lines.%s_ab.txt'\n" % (alias))
    par_out.close()

    cmd = 'cp ./MOOGFEB2017/ab_%s_ab_2.par ./MOOGFEB2017/ab_%s_ab.par' % (alias, alias)
    os.system(cmd)

    #########################################################################
    # Run MOOG with the parameters needed for the model atmosphere
    #########################################################################
    
    ab_dict = runMOOG_ab(star + '_ab', T, logg, xmetal, micro, alias=alias+'_ab',
                         nions=nions, nameions=nameions)

    return ab_dict

