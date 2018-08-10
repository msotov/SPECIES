import numpy as np
import math, os, glob, time, corner, sys, re, logging
import matplotlib
import matplotlib.pyplot as plt
from isochrones import StarModel
from isochrones.mist import MIST_Isochrone
from astropy.io import ascii
from astropy.stats import sigma_clip

def calc_values(starname, T, err_T, logg, err_logg, met, err_met, photometry, make_plot = True, mass_from_file = False, RA = None, DEC = None, use_emcee = False, is_giant = False, age_limits = None):
    print '\t\tCreating isochrones'
    dar = MIST_Isochrone()

    Teff = (T, max(err_T,100.))
    logg_s = (logg, err_logg)
    feh = (met, err_met)

    if os.path.isfile('./isochrones/' + starname + '_samples.h5') and (mass_from_file == True):
        model = StarModel.load_hdf('./isochrones/' + starname + '_samples.h5')
        print '\t\tModel read from file'

    else:
        print '\t\tBuilding model'
        model  = StarModel(dar, Teff=Teff, logg=logg_s, feh=feh, RA = RA, DEC = DEC, **photometry)
        #model.fit_multinest(overwrite=True, basename=starname+'_chains_',  verbose=False, refit = True)
        if 'parallax' in photometry.keys():
            photometry['max_distance'] = 1000./(photometry['parallax'][0]-3.*photometry['parallax'][1])


        if age_limits != None:
            model.set_bounds(age=age_limits)
        else:
            model.obs.add_limit(feh=(feh[0]-2.*feh[1], feh[0]+2.*feh[1]))

        if use_emcee:

            model.use_emcee = True
            p0 = model.sample_from_prior(1)[0]
            p0[2] = feh[0]
            model.fit(nburn=200, niter=1500, nwalkers=500, p0=p0)

        else:
            model.fit_multinest(overwrite=True, basename=starname+'_chains_',  verbose=False, refit = True, n_live_points=1000)

        model.save_hdf('./isochrones/' + starname + '_samples.h5', overwrite = True)

        print '\t\tFinished with model'

    mass_s = model.samples['mass_0_0']
    age_s = model.samples['age_0']
    logg_new = model.samples['logg_0_0']
    radius_s = model.samples['radius_0_0']
    logL_s = model.samples['logL_0_0']

    age_s = 10.**(age_s)/1E9

    mass_p = np.percentile(mass_s, [16, 50, 84])
    age_p = np.percentile(age_s, [16, 50, 84])
    logg_p = np.percentile(logg_new, [16, 50, 84])
    radius_p = np.percentile(radius_s, [16, 50, 84])
    logL_p = np.percentile(logL_s, [16, 50, 84])

    if make_plot:
        make_cmd_dist_plots(starname, mass_s, age_s, logg_new, T, logg, radius_s, met, err_T, err_logg, model)
        data = np.vstack([mass_s, age_s, logg_new, radius_s]).T
        figure = corner.corner(data, labels = [r"Mass", r"Age", r"logg", r"Radius"],\
                              truths = [mass_p[1], age_p[1], logg_p[1], radius_p[1]])
        figure.savefig('./isochrones/plots/' + starname + '_mass_age.pdf')
        plt.close('all')

        plot_keys = [u'B_mag_0_0', u'H_mag_0_0', u'J_mag_0_0', u'K_mag_0_0',
                     u'Teff_0_0', u'V_mag_0_0', u'age_0_0',
                     u'logL_0_0', u'logg_0_0', u'mass_0_0', u'radius_0_0',
                     u'feh_0', u'distance_0', u'AV_0']
        plot_names = ['B', 'H', 'J', 'K', 'Teff', 'V', 'age', 'logL', 'logg', 'mass',
                      'radius', 'feh', 'distance', 'AV']

        truth_vals = [None for i in plot_keys]
        for i,k in enumerate(plot_keys):
            if plot_names[i] in photometry.keys():
                truth_vals[i] = photometry[plot_names[i]][0]
            elif (plot_names[i] == 'distance') and ('parallax' in photometry.keys()):
                truth_vals[i] = 1000./photometry['parallax'][0]
            elif plot_names[i] == 'Teff':
                truth_vals[i] = Teff[0]
            elif plot_names[i] == 'logg':
                truth_vals[i] = logg_s[0]
            elif plot_names[i] == 'feh':
                truth_vals[i] = feh[0]

        fig = corner_plot(model, plot_keys, truths = truth_vals, show_titles=True, quantiles=[0.16,0.5,0.84],\
                          labels=plot_names, sigma_clip=False)
        fig.savefig('./isochrones/plots/' + starname + '_corner_plot.pdf')

        plt.close('all')

        fig = corner_plot(model, plot_keys, truths = truth_vals, show_titles=True, quantiles=[0.16,0.5,0.84],\
                          labels=plot_names, sigma_clip=True)
        fig.savefig('./isochrones/plots/' + starname + '_corner_plot_w_sigma_clip.pdf')

        plt.close('all')

        del data, plot_keys, plot_names


    del dar, Teff, logg_s, model, mass_s, age_s, logg_new, radius_s

    os.system('rm -f chains/' + starname + '_chains_*')


    return mass_p[1], mass_p[2]-mass_p[1], mass_p[1]-mass_p[0],\
           age_p[1], age_p[2]-age_p[1], age_p[1]-age_p[0],\
           logg_p[1], logg_p[2]-logg_p[1], logg_p[1]-logg_p[0],\
           radius_p[1], radius_p[2]-radius_p[1], radius_p[1]-radius_p[0],\
           logL_p[1], logL_p[2]-logL_p[1], logL_p[1]-logL_p[0]


def find_mass_age(starname, T, logg, met, err_T, err_logg, err_met, exception, \
                  photometry, debug, log_f, make_plot = True, mass_from_file = False,\
                  age_limits = None):

    if debug:
        log_f.debug('Computing the mass, age, and photometric logg.')

    if exception == 1:
        photometry_mass = {}
        mags = ['B', 'V', 'R', 'I', 'J', 'H', 'K', 'b', 'y', 'Bt', 'Vt']
        photo_keys = photometry.keys()
        Av = []
        for m in mags:
            if m in photo_keys:
                c = photometry[m][0] - photometry[m][4]
                photometry_mass[m] = (c, max(photometry[m][1], 0.01))
                Av.append(photometry[m][4])

        if 'parallax' in photo_keys:
            photometry_mass['parallax'] = (photometry['parallax'][0].value, \
                                           photometry['parallax'][1].value)

        if np.mean(np.array(Av)) != 0.0:
            if 'parallax' in photo_keys:
                photometry_mass['maxAV'] = 1.0
            else:
                photometry_mass['maxAV'] = 1.0

        if 'RA' in photo_keys:
            RA = photometry['RA'][0]
        else:
            RA = None
        if 'DEC' in photo_keys:
            DEC = photometry['DEC'][0]
        else:
            DEC = None

        if debug:
            log_f.debug('Values used for finding mass and age are: '\
                        '%s, %f, %f, %f, %f, %f, %f',\
                        starname, T, logg, met, err_T, err_logg, err_met)
            log_f.debug('Photometry used is: %s', photometry_mass)

        if (-2.5 <= met <= 0.5) == False:
            print '\t\t[Fe/H] out of bounds, no extrapolation curve created'
            if met > 0.5: met = 0.5
            else: met = -2.5
            print '\t\tIsochrones created for [Fe/H] = ' + str(met)

        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, \
                err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2, \
                logL, err_logL1, err_logL2 =\
                calc_values(starname, T, err_T, logg, err_logg, met, err_met,
                            photometry_mass, make_plot = make_plot,
                            mass_from_file = mass_from_file, RA = RA, DEC = DEC,\
                            age_limits = age_limits)

        del photometry_mass

    else:
        mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, \
                        err_logg_iso2, radius, err_radius1, err_radius2, \
                        logL, err_logL1, err_logL2 = \
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                        0.0, 0.0, 0.0
        if debug:
            log_f.debug('Stopping the calculation because exception = 2.')

    if debug == True:
        log_f.debug('Obtained M=%f + %f - %f, Age=%f + %f - %f, logg=%f + %f - %f, radius=%f + %f - %f',\
                    mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, \
                            err_logg_iso1, err_logg_iso2, radius, err_radius1, err_radius2)


    return mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1, err_logg_iso2,\
           radius, err_radius1, err_radius2, logL, err_logL1, err_logL2


def make_cmd_dist_plots(starname, mass, age, logg, T, logg_s, radius, feh, err_T, err_logg_s, model):

    imass = sigma_clip(mass)
    iage = sigma_clip(age)
    iradius = sigma_clip(radius)

    new_mass = mass[~imass.mask]
    new_age = age[~iage.mask]
    new_radius = radius[~iradius.mask]

    mass_p = np.percentile(new_mass, [16, 50, 84])
    age_p = np.percentile(new_age, [16, 50, 84])
    logg_p = np.percentile(logg, [16, 50, 84])
    radius_p = np.percentile(new_radius, [16, 50, 84])

    feh_p = np.percentile(model.samples['feh_0'], [16, 50, 84])
    logL_p = np.percentile(model.samples['logL_0_0'], [16, 50, 84])

    dar = MIST_Isochrone()
    ev = dar.evtrack(mass_p[1], feh = feh_p[1])
    masses = np.arange(mass_p[1] - 0.6, mass_p[1] + 0.6, 0.2)
    ev_a = []
    for m in masses:
        ev_i = dar.evtrack(m, feh_p[1])
        ev_a.append(ev_i)
        del ev_i

    plt.style.use(['classic'])
    fontname = 'Courier New'
    matplotlib.rcParams.update({'font.family': fontname, 'font.weight': 'medium'})
    ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=10, weight='medium', stretch='normal')

    fig = plt.figure(figsize = (7,10))
    ax1 = plt.subplot2grid((9,1), (0, 0), rowspan = 3)
    ax2 = plt.subplot2grid((9,1), (3, 0), rowspan = 2)
    ax3 = plt.subplot2grid((9,1), (5, 0), rowspan = 2)
    ax4 = plt.subplot2grid((9,1), (7, 0), rowspan = 2)

    ax2.hist(new_mass, histtype = 'step', color = 'black', bins = 20, normed = True)
    ax3.hist(new_age, histtype = 'step', color = 'black', bins = 20, normed = True)
    ax4.hist(new_radius, histtype = 'step', color = 'black', bins = 20, normed = True)

    ax2.axvline(mass_p[0], color = 'red', ls = '--')
    ax2.axvline(mass_p[1], color = 'red', ls = '-')
    ax2.axvline(mass_p[2], color = 'red', ls = '--')

    ax3.axvline(age_p[0], color = 'red', ls = '--')
    ax3.axvline(age_p[1], color = 'red', ls = '-')
    ax3.axvline(age_p[2], color = 'red', ls = '--')

    ax4.axvline(radius_p[0], color = 'red', ls = '--')
    ax4.axvline(radius_p[1], color = 'red', ls = '-')
    ax4.axvline(radius_p[2], color = 'red', ls = '--')


    ax2.set_xlabel(r'Mass ($M_{\odot}$)', fontname = fontname)
    ax3.set_xlabel('Age (Gyr)', fontname = fontname)
    ax4.set_xlabel(r'Radius ($R_{\odot}$)', fontname = fontname)

    #dar = MIST_Isochrone()
    #masses = np.arange(mass_p[1] - 0.6, mass_p[1] + 0.6, 0.2)
    #for m in masses:
    #    ev = dar.evtrack(m, feh = feh)
    for ev_i in ev_a:
        ax1.plot(ev_i['Teff'],ev_i['logg'], color = 'grey', ls = '-',
                alpha = 0.5, marker = '.',
                markersize = 3)

    #ev = dar.evtrack(mass_p[1], feh = feh)
    ax1.plot(ev['Teff'],ev['logg'], color = 'black', marker = '.')

    ax1.errorbar(T, logg_s, xerr = err_T, yerr = err_logg_s, marker = 'o', color = 'red')
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    #ax1.set_ylim(logg_s + 5.*err_logg_s, logg_s - 5.*err_logg_s)
    #ax1.set_xlim(T + 5.*err_T, T - 5.*err_T)
    ax1.set_xlabel(r'$T_{eff}$ (K)', fontname = fontname)
    ax1.set_ylabel(r'logg (cm s$^{-2}$)', fontname = fontname)

    for ax in (ax1, ax2, ax3, ax4):
        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)


    plt.subplots_adjust(wspace=0.7, hspace=1.0, left = 0.095, right = 0.96,
                        top = 0.98, bottom = 0.07)

    fig.savefig('./isochrones/plots/%s_CMD_dist.pdf' % starname)

    plt.close('all')

    fig, ax = plt.subplots(1,2,figsize=(12,5))
    ax[0].plot(ev['Teff'], ev['logg'], color = 'black')
    ax[0].errorbar(T, np.median(model.samples['logg_0_0']), xerr = err_T, yerr = np.std(model.samples['logg_0_0']), marker = 'o', color = 'red')
    ts, te = ax[0].get_xlim()
    ls, le = ax[0].get_ylim()

    for ev_i in ev_a:
        ax[0].plot(ev_i['Teff'], ev_i['logg'], color = 'black', alpha=0.2)

    ax[0].set_xlim(ts,te)
    ax[0].set_ylim(ls,le)

    ax[0].invert_xaxis()
    ax[0].invert_yaxis()

    xerr_i = [np.array([err_T]), np.array([err_T])]
    yerr_i = [np.array([logL_p[1]-logL_p[0]]), np.array([logL_p[2]- logL_p[1]])]


    ax[1].plot(ev['Teff'], ev['logL'], color = 'black')
    ax[1].errorbar(np.array([T]), np.array([logL_p[1]]), xerr = xerr_i, yerr = yerr_i, marker = 'o', color = 'red')
    ts, te = ax[1].get_xlim()
    ls, le = ax[1].get_ylim()

    for ev_i in ev_a:
        ax[1].plot(ev_i['Teff'], ev_i['logL'], color = 'black', alpha=0.2)

    ax[1].set_xlim(ts,te)
    ax[1].set_ylim(ls,le)

    ax[1].invert_xaxis()
    #ax[1].invert_yaxis()

    ax[0].set_xlabel(r'$T_{eff}$ (K)', fontname = fontname)
    ax[0].set_ylabel(r'logg (cm s$^{-2}$)', fontname = fontname)

    ax[1].set_xlabel(r'$T_{eff}$ (K)', fontname = fontname)
    ax[1].set_ylabel(r'logL ($L_{\odot}$)', fontname = fontname)

    fig.savefig('./isochrones/plots/%s_CMD.pdf' % starname)

    plt.close('all')

    del dar, ev, masses, fig, ax1, ax2, ax3, ax4, mass_p, age_p, logg_p, radius_p, imass, iage, iradius, new_mass, new_age, new_radius, ev_a, model


def corner_plot(model, params, query=None, **kwargs):

    df = model.random_samples(5000)
    if query is not None:
        df = df.query(query)

    priors = []
    ranges_plot = []
    for p in params:
        isigma = sigma_clip(df[p])
        ranges_plot.append((min(df[p][~isigma.mask]), max(df[p][~isigma.mask])))
        del isigma
        if re.match('mass', p):
            priors.append(lambda x: model.prior('mass', x, bounds=model.bounds('mass')))
        elif re.match('age', p):
            priors.append(lambda x: model.prior('age', x, bounds=model.bounds('age')))
        elif re.match('feh', p):
            priors.append(lambda x: model.prior('feh', x, bounds=model.bounds('feh')))
        elif re.match('distance', p):
            priors.append(lambda x: model.prior('distance', x, bounds=model.bounds('distance')))
        elif re.match('AV', p):
            priors.append(lambda x: model.prior('AV', x, bounds=model.bounds('AV')))
        else:
            priors.append(None)

    if 'sigma_clip' in kwargs.keys():
        if kwargs['sigma_clip'] == False:
            ranges_plot=None
    try:
        fig = corner.corner(df[params], priors=priors, range = ranges_plot, **kwargs)
    except:
        logging.warning("Use Tim's version of corner to plot priors.")
        fig = corner.corner(df[params], range = ranges_plot, **kwargs)
    fig.suptitle(model.name, fontsize=22)
    return fig
