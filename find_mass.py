import numpy as np
import math, os, glob, time, corner, sys
import matplotlib
import matplotlib.pyplot as plt
from isochrones import StarModel
from isochrones.mist import MIST_Isochrone
from astropy.io import ascii


def calc_values(starname, T, err_T, logg, err_logg, met, err_met, photometry, make_plot = True, mass_from_file = False):
    print '\t\tCreating isochrones'
    dar = MIST_Isochrone()

    Teff = (T, err_T)
    logg_s = (logg, err_logg)
    feh = (met, err_met)

    if os.path.isfile('./dotter_isochrones/' + starname + '_samples.h5') and (mass_from_file == True):
        model = StarModel.load_hdf('./dotter_isochrones/' + starname + '_samples.h5')
        print '\t\tModel read from file'

    else:
        print '\t\tBuilding model'

        model  = StarModel(dar, Teff=Teff, logg=logg_s, feh=feh, **photometry)
        model.fit_multinest(overwrite=True, basename=starname+'_chains_',  verbose=False, refit = True)

        model.save_hdf('./dotter_isochrones/' + starname + '_samples.h5', overwrite = True)

        print '\t\tFinished with model'

    mass_s = model.samples['mass_0_0']
    age_s = model.samples['age_0']
    logg_new = model.samples['logg_0_0']
    radius_s = model.samples['radius_0_0']

    age_s = 10.**(age_s)/1E9

    # Apply offsets with the sun values
    age_s = age_s - 0.9
    logg_new = logg_new + 0.02
    radius_s = radius_s - 0.01


    mass_p = np.percentile(mass_s, [16, 50, 84])
    age_p = np.percentile(age_s, [16, 50, 84])
    logg_p = np.percentile(logg_new, [16, 50, 84])
    radius_p = np.percentile(radius_s, [16, 50, 84])

    if make_plot:
        make_cmd_dist_plots(starname, mass_s, age_s, logg_new, T, logg, radius_s, met, err_T, err_logg)
        data = np.vstack([mass_s, age_s, logg_new, radius_s]).T
        figure = corner.corner(data, labels = [r"Mass", r"Age", r"logg", r"Radius"],\
                              truths = [mass_p[1], age_p[1], logg_p[1], radius_p[1]])
        figure.savefig('./dotter_isochrones/plots/' + starname + '_mass_age.pdf')
        plt.close('all')

        del data

    mass = mass_p[1]
    err_mass = max(abs(mass_p[0] - mass_p[1]), abs(mass_p[1] - mass_p[2]))
    age = age_p[1]
    err_age = max(abs(age_p[0] - age_p[1]), abs(age_p[1] - age_p[2]))
    logg_photo = logg_p[1]
    err_logg_photo = max(abs(logg_p[0] - logg_p[1]), abs(logg_p[1] - logg_p[2]))
    radius = radius_p[1]
    err_radius = max(abs(radius_p[0] - radius_p[1]), abs(radius_p[1] - radius_p[2]))

    del dar, Teff, logg_s, model, mass_s, age_s, logg_new, mass_p, age_p, logg_p

    os.system('rm -f chains/' + starname + '_chains_*')

    return mass, err_mass, age, err_age, logg_photo, err_logg_photo, radius, err_radius



def find_mass_age(starname, T, logg, met, err_T, err_logg, err_met, photometry, make_plot = True, mass_from_file = False):

    if (-2.5 <= met <= 0.5) == False:
        print '\t\t[Fe/H] out of bounds, no extrapolation curve created'
        if met > 0.5: met = 0.5
        else: met = -2.5
        print '\t\tIsochrones created for [Fe/H] = ' + str(met)

    if err_T>500.: err_T = 500.

    mass, err_mass, age, err_age, logg_photo, err_logg_photo, radius, err_radius =\
                calc_values(starname, T, err_T, logg, err_logg, met, err_met,
                            photometry, make_plot = make_plot,
                            mass_from_file = mass_from_file)


    return mass, err_mass, age, err_age, logg_photo, err_logg_photo, radius, err_radius


def make_cmd_dist_plots(starname, mass, age, logg, T, logg_s, radius, feh, err_T, err_logg_s):

    plt.style.use(['classic'])
    fontname = 'Courier New'
    matplotlib.rcParams.update({'font.family': fontname, 'font.weight': 'medium'})
    ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=10, weight='medium', stretch='normal')

    mass_p = np.percentile(mass, [16, 50, 84])
    age_p = np.percentile(age, [16, 50, 84])
    logg_p = np.percentile(logg, [16, 50, 84])
    radius_p = np.percentile(radius, [16, 50, 84])

    fig = plt.figure(figsize = (7,10))
    ax1 = plt.subplot2grid((9,1), (0, 0), rowspan = 3)
    ax2 = plt.subplot2grid((9,1), (3, 0), rowspan = 2)
    ax3 = plt.subplot2grid((9,1), (5, 0), rowspan = 2)
    ax4 = plt.subplot2grid((9,1), (7, 0), rowspan = 2)

    ax2.hist(mass, histtype = 'step', color = 'black', bins = 20, normed = True)
    ax3.hist(age, histtype = 'step', color = 'black', bins = 20, normed = True)
    #ax4.hist(logg, histtype = 'step', color = 'black', bins = 20, normed = True)
    ax4.hist(radius, histtype = 'step', color = 'black', bins = 20, normed = True)

    ax2.axvline(mass_p[0], color = 'red', ls = '--')
    ax2.axvline(mass_p[1], color = 'red', ls = '-')
    ax2.axvline(mass_p[2], color = 'red', ls = '--')

    ax3.axvline(age_p[0], color = 'red', ls = '--')
    ax3.axvline(age_p[1], color = 'red', ls = '-')
    ax3.axvline(age_p[2], color = 'red', ls = '--')

    #ax4.axvline(logg_p[0], color = 'red', ls = '--')
    #ax4.axvline(logg_p[1], color = 'red', ls = '-')
    #ax4.axvline(logg_p[2], color = 'red', ls = '--')
    ax4.axvline(radius_p[0], color = 'red', ls = '--')
    ax4.axvline(radius_p[1], color = 'red', ls = '-')
    ax4.axvline(radius_p[2], color = 'red', ls = '--')


    ax2.set_xlabel(r'Mass ($M_{\odot}$)', fontname = fontname)
    ax3.set_xlabel('Age (Gyr)', fontname = fontname)
    #ax4.set_xlabel(r'logg (cm s$^{-2}$)', fontname = fontname)
    ax4.set_xlabel(r'Radius ($R_{\odot}$)', fontname = fontname)

    dar = MIST_Isochrone()
    masses = np.arange(mass_p[1] - 0.5, mass_p[1] + 0.5, 0.01)
    for m in masses:
        ev = dar.evtrack(m, feh = feh, dage = 0.02)
        ax1.plot(ev['Teff'],ev['logg'], color = 'grey', ls = '-',
                alpha = 0.5, marker = '.',
                markersize = 3)

        del ev

    ev = dar.evtrack(mass_p[1], feh = feh, dage = 0.02)
    ax1.plot(ev['Teff'],ev['logg'], color = 'black', marker = '.')

    ax1.errorbar(T, logg_s, xerr = err_T, yerr = err_logg_s, marker = 'o', color = 'red')
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylim(logg_s + 2.*err_logg_s, logg_s - 2.*err_logg_s)
    ax1.set_xlim(T + 2.*err_T, T - 2.*err_T)
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

    fig.savefig('./dotter_isochrones/plots/%s_CMD_dist.pdf' % starname)
    #fig.savefig('./dotter_isochrones/plots/%s_CMD_dist.ps' % starname)

    del dar, ev, masses, fig, ax1, ax2, ax3, ax4, mass_p, age_p, logg_p, radius_p

    plt.close('all')


#if __name__ == '__main__':
    #start = time.time()
    #mass, err_mass, age, err_age, logg, err_logg, n = find_mass_age('sun01_harps', 5750., 4.75, -0.01, 40., 0.005)
    #print mass, err_mass
    #print age, err_age
    #print logg, err_logg
    #print n
    #end = time.time()

    #print (end-start)/60.
