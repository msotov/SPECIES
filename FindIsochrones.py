from __future__ import print_function
from __future__ import division
from builtins import str
import os
import sys
import logging
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['seaborn-muted'])
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from isochrones import SingleStarModel, get_ichrone
from isochrones.mist import MIST_Isochrone
from isochrones.priors import FlatPrior


#@profile
def calc_values(starname, T, err_T, logg, err_logg, met, err_met, photometry,\
                makeplot=True, mass_from_file=False,\
                is_giant=False, age_limits=None,\
                PATH_ISO='./isochrones', Av=[]):
    print('\t\tCreating isochrones')

    mist = get_ichrone('mist')
    props = {'Teff': (T, min(err_T, 250)),
             'logg': (logg, min(err_logg, 0.5)),
             'feh': (met, err_met)}
    props.update(photometry)
    
    #props.update({'density': (2.72, 0.16)})

    if os.path.isfile(os.path.join(os.path.relpath(PATH_ISO, '.'), '%s_samples.h5' % starname))\
        and mass_from_file:
        model = SingleStarModel.load_hdf(os.path.join(os.path.relpath(PATH_ISO, '.'),
                                                      '%s_samples.h5' % starname))
        print('\t\tModel read from file')
        makeplot = False

    else:
        model = SingleStarModel(mist, name=starname, **props)
        model._priors['eep'].bounds = (200, 1710)

        if is_giant:
            model._priors['eep'].bounds = (353, 1710)
            
        logging.info('EEP prior boundaries set to (%d, %d)', *model._priors['eep'].bounds)
        if Av:
            Avmax = min(max(np.max(Av), 0.1), 1.0)
            #Avmax = 1.0
            model.set_prior(AV=FlatPrior((0.0, Avmax)))
            logging.info('Prior for Av changed to flat prior with boundaries: (%.1f, %.1f)',
                         0, Avmax)
        logging.info(props)

        model.fit(overwrite=True, basename=starname+'_chains_',\
                  verbose=False, refit=True, n_live_points=1000)
        model.save_hdf(os.path.join(os.path.relpath(PATH_ISO, '.'), '%s_samples.h5' % starname),
                       overwrite=True)

    derived = model.derived_samples
    mass_s = derived['mass']
    age_s = derived['age']
    logg_s = derived['logg']
    radius_s = derived['radius']
    logL_s = derived['logL']
    AV = derived['AV']
    EEP = derived['eep']

    age_s = 10**(age_s)/(1E9)

    mass_p = np.percentile(mass_s, [16, 50, 84])
    age_p = np.percentile(age_s, [16, 50, 84])
    logg_p = np.percentile(logg_s, [16, 50, 84])
    radius_p = np.percentile(radius_s, [16, 50, 84])
    logL_p = np.percentile(logL_s, [16, 50, 84])
    AV_p = np.median(AV)
    EEP_p = np.percentile(EEP, [16, 50, 84])

    P_preMS, P_MS, P_RGB, P_HB, P_postHB = Pevol(EEP)

    if makeplot:
        try:
            fig = model.corner_physical()
            fig.savefig(os.path.join(os.path.relpath(PATH_ISO, '.'),
                                     'plots', '%s_corner_physical.pdf' % starname))
            plt.close('all')
            del fig

            fig = model.corner_observed()
            fig.savefig(os.path.join(os.path.relpath(PATH_ISO, '.'),
                                     'plots', '%s_corner_observed.pdf' % starname))
            plt.close('all')
            del fig

            fig = model.corner_derived(['Teff', 'logg', 'age', 'mass',
                                        'radius', 'age', 'distance', 'eep', 'logL'])
            fig.savefig(os.path.join(os.path.relpath(PATH_ISO, '.'),
                                     'plots', '%s_corner_derived.pdf' % starname))
            plt.close('all')
            del fig

            plot_isochrones(starname, derived, props, N=10, PATH_ISO=PATH_ISO)
        except RuntimeError:
            logging.error('Error while plotting. Skipping.')

    os.system('rm -f chains/' + starname + '_chains_*')

    del model, derived, mass_s, age_s, logg_s, radius_s, logL_s, AV, EEP


    return mass_p[1], mass_p[2]-mass_p[1], mass_p[1]-mass_p[0],\
           age_p[1], age_p[2]-age_p[1], age_p[1]-age_p[0],\
           logg_p[1], logg_p[2]-logg_p[1], logg_p[1]-logg_p[0],\
           radius_p[1], radius_p[2]-radius_p[1], radius_p[1]-radius_p[0],\
           logL_p[1], logL_p[2]-logL_p[1], logL_p[1]-logL_p[0], AV_p,\
           EEP_p[1], EEP_p[2]-EEP_p[1], EEP_p[1]-EEP_p[0],\
           P_preMS, P_MS, P_RGB, P_HB, P_postHB


def Pevol(eep):
    i = np.where(~np.isnan(eep))[0]
    N = i.size
    if N == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    N_preMS = np.where(eep[i] < 202)[0].size
    N_MS = np.where((eep[i] >= 202) & (eep[i] <= 454))[0].size
    N_RGB = np.where((eep[i] > 454) & (eep[i] < 631))[0].size
    N_HB = np.where((eep[i] >= 631) & (eep[i] <= 707))[0].size
    N_postHB = np.where(eep[i] > 707)[0].size

    return N_preMS/N, N_MS/N, N_RGB/N, N_HB/N, N_postHB/N


def plot_isochrones(starname, derived, props, N=10, PATH_ISO='.'):
    try:
        mist = MIST_Isochrone()

        age_samples = np.random.choice(derived['age'], size=N)
        feh_samples = np.random.choice(derived['feh'], size=N)

        fig, ax = plt.subplots()
        for i in range(10):
            grid = mist.isochrone(age_samples[i], feh_samples[i])
            ax.plot(grid['logTeff'], grid['logL'],
                    color='gray', alpha=0.4, lw=0.8, label='_nolegend_')
            del grid

        age_mean = np.median(derived['age'])
        feh_mean = np.median(derived['feh'])
        grid = mist.isochrone(age_mean, feh_mean)
        ll = ax.scatter(grid['logTeff'], grid['logL'], c=grid['mass'],
                        cmap='viridis', s=8, label='_nolegend_', edgecolors='face')
        cbar = fig.colorbar(ll, ax=ax)
        cbar.set_label(r'$M_{\star}$', fontsize='x-large')
        del grid

        t16, t50, t84 = np.percentile(derived['logTeff'], [16, 50, 84])
        l16, l50, l84 = np.percentile(derived['logL'], [16, 50, 84])

        ax.errorbar(t50, l50, xerr=[[t50-t16], [t84-t50]], yerr=[[l50-l16], [l84-l50]],
                    marker='o', color='red', markeredgecolor='red', markersize=10,
                    label=r'$T$ = %.0f K'
                          '\n'
                          r'[Fe/H] = %.2f'
                          '\n'
                          r'$M$ = %.2f $M_{\odot}$'
                          '\n'
                          r'$R$ = %.2f $R_{\odot}$' % (props['Teff'][0], props['feh'][0],\
                                                       np.median(derived['mass']),\
                                                       np.median(derived['radius'])))

        ax.invert_xaxis()
        ax.legend(numpoints=1, loc='lower left')

        ax.set_xlabel('logTeff', fontsize='x-large')
        ax.set_ylabel('logL', fontsize='x-large')

        fig.savefig(os.path.join(os.path.relpath(PATH_ISO, '.'),
                                 'plots', '%s_isochrones.pdf' % starname))

        plt.close('all')

        del mist, age_samples, feh_samples, fig

    except Exception as e:
        logging.error('ERROR IN CODE, STOPPING THE COMPUTATION')
        _, _, exc_tb = sys.exc_info()
        logging.error('line %d: %s', exc_tb.tb_lineno, e)

    return 0


def find_mass_age(starname, T, logg, met, err_T, err_logg, err_met, exception, \
                  photometry, make_plot=True, mass_from_file=False,\
                  age_limits=None, is_giant=False, PATH_ISO='./isochrones',\
                  skip_mass=False, plot_mass=False, Av=0.0):

    if skip_mass:
        logging.warning('Skipping mass computation.')
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    err_T = max(err_T, 0.01)
    err_logg = max(err_logg, 0.01)
    err_met = max(err_met, 0.01)
    AV = 0.0

    use_mags = True
    logging.info('Computing the mass, age, and photometric logg.')

    # Check that PATH_ISO contains a directory called '/plots'
    # If not, create it.

    if not os.path.exists(os.path.join(os.path.relpath(PATH_ISO, '.'), 'plots')):
        os.makedirs(os.path.join(os.path.relpath(PATH_ISO, '.'), 'plots'))

    if exception == 1:
        photometry_mass = {}
        mags = ['J', 'H', 'K', 'B', 'V', 'G']#, 'W1', 'W2', 'W3']
        #mags = ['B', 'V', 'G']#, 'W1', 'W2', 'W3']
        mags_Av = ['B', 'V', 'Vt', 'Bt']
        Av = [Av]
        if is_giant:
            use_mags = False
        for m in mags_Av:
            if m in photometry:
                if photometry[m][4] > 1.0:
                    use_mags = False
        if not use_mags:
            mags = ['J', 'H', 'K', 'G']
        for m in mags:
            if m in photometry:
                #if 'HIPPARCOS' not in photometry[m][3]: 
                c = photometry[m][0]# - photometry[m][4]
                if photometry[m][1] <= 0.3:
                    photometry_mass[m] = (c, max(photometry[m][1], 0.01))
                    if photometry[m][4] > 0.0 and photometry[m][4] < 1.0:
                        Av.append(photometry[m][4])
                    elif photometry[m][4] <= 0.0:
                        Av.append(0.1)
        if 'parallax' in photometry:
            photometry_mass['parallax'] = (photometry['parallax'][0].value, \
                                           photometry['parallax'][1].value)

        logging.info('Values used for finding mass and age are: '\
                        '%s, %f, %f, %f, %f, %f, %f',\
                        starname, T, logg, met, err_T, err_logg, err_met)
        logging.info('Photometry used is: %s', photometry_mass)

        if (-2.5 <= met <= 0.5) is False:
            print('\t\t[Fe/H] out of bounds, no extrapolation curve created')
            met = 0.5 if met > 0.5 else -2.5
            print('\t\tIsochrones created for [Fe/H] = %s' % str(met))


        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, \
                err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2, \
                logL, err_logL1, err_logL2, AV, eep, err_eep1, err_eep2, \
                p_prems, p_ms, p_rgb, p_hb, p_posthb =\
                calc_values(starname, T, err_T, logg, err_logg, met, err_met,
                            photometry_mass, makeplot=make_plot,
                            mass_from_file=mass_from_file, PATH_ISO=PATH_ISO,
                            is_giant=is_giant, age_limits=age_limits, Av=Av)

        del photometry_mass

    else:
        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, \
                        err_logg_iso2, radius, err_radius1, err_radius2, \
                        logL, err_logL1, err_logL2, AV, eep, err_eep1, err_eep2, \
                        p_prems, p_ms, p_rgb, p_hb, p_posthb = \
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        logging.warning('Stopping the calculation because exception = 2.')

    logging.info('Obtained M=%f + %f - %f, Age=%f + %f - %f, '
                 'logg=%f + %f - %f, radius=%f + %f - %f',\
                 mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, \
                 err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2)


    return mass, err_mass1, err_mass2, age, err_age1, err_age2,\
           logg_iso, err_logg_iso1, err_logg_iso2,\
           radius, err_radius1, err_radius2, logL, err_logL1, err_logL2, AV,\
           eep, err_eep1, err_eep2, p_prems, p_ms, p_rgb, p_hb, p_posthb
