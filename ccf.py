"""
Computes the CCF between the spectra with a template
(binary mask) to correct to restframe.

Created: 24/08/2015
"""
from __future__ import print_function
from __future__ import division

from builtins import range
import os
import collections
from past.utils import old_div
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.io.fits import getheader
from scipy.optimize import curve_fit
from PyAstronomy import pyasl

plt.style.use(['seaborn-muted'])
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def compute_snr(x, data):
    """
    Computes the S/N for spectra using a set of ranges where there shouldn't be
    any lines, only continuum.
    """
    #ranges = [[5840.22, 5844.01], [5979.25, 5983.13], [6066.49, 6075.55],\
    #          [6195.95, 6198.74], [6386.36, 6392.14], [5257.83, 5258.58],\
    #          [5256.03, 5256.7]]
    ranges = [[5979.25, 5983.13], [6066.49, 6075.55],\
              [6195.95, 6198.74], [6386.36, 6392.14], [5257.83, 5258.58],\
              [5256.03, 5256.7]]
    sn = []
    beq = pyasl.BSEqSamp()

    #fig, ax = plt.subplots(2, 3)
    #ax_flat = ax.flatten()

    for ii, r in enumerate(ranges):
        data_range = data[np.where((x >= r[0]) & (x <= r[1]))[0]]
        if data_range.size > 4:
            if np.min(data_range) < 0.0:
                ymin = np.min(data_range)
                data_range -= ymin
            m = np.median(data_range)
            s1, _ = beq.betaSigma(data_range, 1, 1, returnMAD=True)
            n = len(data_range)
            s2 = 0.6052697*np.median(np.abs(2.0*data_range[2:n-2]-data_range[0:n-4]-data_range[4:n]))
            s = np.std(data_range)
            #print(m/s1, m/s2, m/s)
            sn.append(old_div(m, s1))

            #axi = ax_flat[ii]
            #axi.plot(x[np.where((x >= r[0]) & (x <= r[1]))[0]], data_range)
            #axi.set_title('m = %.3f, s = %.3f' % (m, s1))

            del m, s
        del data_range

    #plt.tight_layout()
    #fig.savefig('test_snr.pdf')
    #plt.close('all')

    del ranges
    inonan = np.where(~np.isnan(sn))[0]
    return min(np.median(np.array(sn)[inonan]), 300)

########################################################

def fit_func(t, p0, p1, p2, p3):
    return p3 + p0*np.exp(old_div(-(t-p1)**2., (2.*p2**2.)))

########################################################

def continuum_det(x, y): # Substracts the continuum and normalizes
    rejt = 0.98
    p = np.poly1d(np.polyfit(x, y, 2))
    yfit = p(x)

    for _ in range(3):
        dif = np.hstack((np.abs((y[:-1]-y[1:])/y[:-1]), [1.0]))
        i = np.where(((y-yfit*rejt) > 0) & (dif < 0.1))[0]
        vecx = x[i]
        vecy = y[i]
        p = np.poly1d(np.polyfit(vecx, vecy, 2))
        yfit = p(x)

    ynorm = y/p(x) - 1.0
    del yfit, p
    return ynorm

########################################################

def ccf(nombre, x_range, data_range2):
    make_plot = True
    rvrange = 200.0

    c = 299792.458

    x1 = np.arange(x_range[0]-rvrange, x_range[0], x_range[1]-x_range[0])
    x2 = np.arange(x_range[-1], x_range[-1]+rvrange, x_range[-1]-x_range[-2])
    xtem = np.hstack([x1, x_range, x2])

    lines1, lines2, flux_l = np.genfromtxt('./binary_masks/G2.mas', dtype=None, unpack=True)
    ilines = np.where((lines1 > x_range[0]) & (lines2 < x_range[-1]))[0]
    lines1_new = lines1[ilines]
    lines2_new = lines2[ilines]
    flux_l_new = flux_l[ilines]

    ftem = np.zeros(xtem.size)

    for i, f in enumerate(flux_l_new):
        indices = np.where((xtem >= lines1_new[i]) & (xtem <= lines2_new[i]))[0]
        if indices.size > 0:
            ftem[indices] = f
        del indices

    rv_temp, cc = pyasl.crosscorrRV(x_range, data_range2, xtem, ftem, -rvrange, rvrange, 0.05)

    del lines1, lines2, flux_l, ilines, lines1_new, lines2_new, flux_l_new

    index_min = np.argmin(cc)
    x_min = rv_temp[index_min]
    y_min = cc[index_min]
    popt, _ = curve_fit(fit_func, rv_temp, cc, p0=(y_min, x_min, 1., 0.))
    final_rv = popt[1]
    y_fit = fit_func(rv_temp, *popt)
    #final_rv = -3.614

    if make_plot:
        ticks_font = matplotlib.font_manager.FontProperties(style='normal',\
                                                            size=10, weight='medium',
                                                            stretch='normal')

        if not os.path.exists('./Spectra/plots_ccf'):
            os.makedirs('./Spectra/plots_ccf')
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        lines = [6021.8, 6024.06, 6027.06]

        i_range_plot = np.where((x_range >= 6020) & (x_range <= 6030))[0]
        x_range_plot = x_range[i_range_plot]
        data_range_plot = data_range2[i_range_plot]

        ax1.plot(x_range_plot, data_range_plot, color='darkgrey')
        for i in lines:
            ax1.axvline(i, color='red')

        xnew = old_div(x_range_plot, (1. + old_div(final_rv, c)))
        ax1.plot(xnew, data_range_plot, color='black')

        ax2.plot(rv_temp, cc, color='black')
        ax2.plot(rv_temp, y_fit, color='red')
        ax2.set_xlim(popt[1] - 5.*popt[2], popt[1] + 5.*popt[2])
        sy, ey = ax2.get_ylim()
        ax2.set_ylim(sy, 1.0)
        sx, ex = ax2.get_xlim()
        delta_x = (ex-sx)
        sy, ey = ax2.get_ylim()
        delta_y = (ey-sy)
        ax2.text(sx + 2.*delta_x/5., ey - delta_y/5., 'RV = %.2f  km/s' % final_rv)
        ax1.set_xlabel(r'$\lambda$ ($\AA$)', fontsize='medium')
        ax2.set_xlabel(r'Velocity (km/s)', fontsize='medium')

        for ax in (ax1, ax2):
            _ = [i.set_linewidth(0.5) for i in ax.spines.values()]
            for label in ax.get_yticklabels():
                label.set_fontproperties(ticks_font)
            for label in ax.get_xticklabels():
                label.set_fontproperties(ticks_font)

        repeats = collections.Counter(nombre)['/']
        for _ in range(repeats):
            indice = nombre.index('/')
            new_name = nombre[indice + 1:]
            nombre = new_name
        plt.tight_layout()
        fig.savefig('./Spectra/plots_ccf/%s_ccf.pdf' % nombre)
        plt.close('all')

        del fig, ax1, ax2, sy, ey, delta_y, repeats, \
            i_range_plot, x_range_plot, data_range_plot, sx, ex, delta_x

    del x1, x2, xtem, ftem, rv_temp, cc, y_fit, popt

    return final_rv

########################################################

def plot_lines(x, data, lines, name_lines, nombre, savefigure=True, prev_fig=None, figclear=False):
    space = 5.
    ticks_font = matplotlib.font_manager.FontProperties(style='normal', size=10,
                                                        weight='medium', stretch='normal')

    if prev_fig is None:
        fig = plt.figure()
    else:
        fig = prev_fig
    if figclear:
        fig.clear()
    for l, line in enumerate(lines):
        i_lineas = np.where((x > (line - space)) & (x < (line + space)))[0]
        x_lineas = x[i_lineas]
        data_lineas = data[i_lineas]
        nplot = l + 1
        if nplot == 1:
            ax = fig.add_subplot(221)
        elif nplot == 2:
            ax = fig.add_subplot(222)
        elif nplot == 3:
            ax = fig.add_subplot(223)
        else:
            ax = fig.add_subplot(224)

        ax.plot(x_lineas, data_lineas, color='black')
        ax.axvline(line, color='red')
        ax.set_title(name_lines[l])

        _ = [i.set_linewidth(0.5) for i in ax.spines.values()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)

        del i_lineas, x_lineas, data_lineas

    plt.draw()
    plt.tight_layout()

    if savefigure:
        if not os.path.exists('./Spectra/plots_fit'):
            os.makedirs('./Spectra/plots_fit')
        repeats = collections.Counter(nombre)['/']
        for _ in range(repeats):
            indice = nombre.index('/')
            new_name = nombre[indice + 1:]
            nombre = new_name
        fig.savefig('./Spectra/plots_fit/%s_fit_lines.pdf' % nombre)
    else:
        plt.show()
    plt.close('all')
    del fig

########################################################

def restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=True, make_1d=False):
    if header0 is not None:
        if 'CUNIT1' in header0:
            if header0['CUNIT1'] == 'NM':
                wave = wave*10.
                if wave_1d is not None:
                    wave_1d = wave_1d*10.
                header0['CUNIT1'] = 'Angstrom'

        elif np.mean(wave) < 1000.0:
            wave = wave*10.
            if wave_1d is not None:
                wave_1d = wave_1d*10.
            header0['CUNIT1'] = 'Angstrom'

    if wave_1d is None:
        x = wave
        data = flux
    else:
        x = wave_1d
        data = flux_1d
        

    lineas = [6021.8, 6024.06, 6027.06, 6562.81]
    name_lineas = ['6021.8', '6024.06', '6027.06', 'H-alpha']
    
    i_range = np.where((x >= 5000) & (x <= 6250))[0]
    x_range = x[i_range]
    data_range = data[i_range]

    if data_range.size == 0:
        print('\t\tWavelength range to correct to restframe is not present. '\
              'Cannot create the _res.fits file.')

    else:
        if np.mean(data_range) != 0.:
            data_range2 = continuum_det(x_range, data_range)
            rv = ccf(starname, x_range, data_range2)
            del data_range2
        else:
            rv = 0.

        if do_restframe is False:
            rv = 0.

        c = 299792.458
        print('\t\trv = %f km/s' % rv)

        if header0 is not None:
            if 'CUNIT1' in header0:
                if header0['CUNIT1'] == 'NM':
                    wave = wave*10.

        if (old_div(rv, c)) != (-1.):
            x_corr = old_div(wave, (1. + old_div(rv, c)))
            x_corr_1d = old_div(x, (1. + old_div(rv, c)))

        # Saves the plot of the spectrum corrected for restframe
        # in the three lines used by the CCF, plus the h-alpha line.

        plot_fit = True
        if plot_fit:
            plot_lines(x_corr_1d, data, lineas, name_lineas, starname)

        #################################################
        
        if hasattr(snr, "__len__"):
            if len(snr) == 0:
                snr = compute_snr(x_corr_1d, data)
                if np.isnan(snr):
                    snr = None

        if snr is None:
            snr = compute_snr(x_corr_1d, data)
            if np.isnan(snr):
                snr = None

        if len(x_corr.shape) == 1:
            data_corr = np.zeros((2, len(x_corr)))
            data_corr[0, :] = x_corr
            data_corr[1, :] = flux

        else:
            data_corr = np.zeros((2, x_corr.shape[0], x_corr.shape[1]))
            data_corr[0, :, :] = x_corr
            data_corr[1, :, :] = flux


        if header0 is None:
            fits.writeto('%s_res.fits' % starname, data=data_corr, overwrite=True)
            header0 = getheader('%s_res.fits' % starname)
            
        if make_1d:
            if flux_1d is None:
                flux_1d = data
            fits.writeto('%s_res_1d.fits' % starname, data=flux_1d, overwrite=True)
            header0_1d = getheader('%s_res_1d.fits' % starname)
            header0_1d.set('CRPIX1', 1.0)
            header0_1d.set('CDELT1', np.mean(np.diff(x_corr_1d)))
            header0_1d.set('CRVAL1', x_corr_1d[0])
            fits.writeto('%s_res_1d.fits' % starname, data=flux_1d, header=header0_1d,
                         overwrite=True)
            del header0_1d

        try:
            isize = len(snr)
        except:
            isize = 0

        after = 'NAXIS'
        if 'NAXIS2' in header0:
            after = 'NAXIS2'
        if 'NAXIS3' in header0:
            after = 'NAXIS3'

        if isize > 0:
            for o in range(isize):
                header0.set('SNR%d' % o, snr[o], after=after)
        else:
            if snr is not None:
                header0.set('SNR', snr, after=after)

        header0.set('RV', rv, after=after)
        header0['NAXIS'] = 2

        if snr is not None:
            try:
                fits.writeto('%s_res.fits' % starname, data=data_corr, header=header0,
                             overwrite=True)
            except Exception as e:
                print(e)
                fits.writeto('%s_res.fits' % starname, data=data_corr, header=header0[:50],
                             overwrite=True)

        del x_corr, data, x_corr_1d, data_corr

    del x_range, data_range, i_range

    return snr
