'''
Computes the CCF between the spectra with a template
(binary mask) to correct to restframe.

Created: 24/08/2015

To do: -Add an option  to save the CCF, for later use in
        bisector analysis.
       -Add more lines.
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyfits, os, collections
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
#from matplotlib import rcParams
#rcParams['font.family'] = 'serif'

plt.style.use(['classic'])

########################################################
########################################################
########################################################

def values(h, j):# funcion que extrae las longitudes de onda del header
    N = h['NAXIS' + str(j)];
    val = np.zeros(N);
    for i in range(0, N):
        val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)])+float(h['CRVAL'+str(j)]);
    return val;

########################################################
########################################################
########################################################


def fit_func(t, p0, p1, p2, p3):
    return p3 + p0*np.exp(-(t-p1)**2./(2.*p2**2.))


########################################################
########################################################
########################################################

def cont(t, p_0, p_1, p_2):
    return p_0 + p_1*t + p_2*t*t

########################################################
########################################################
########################################################

def continuum_det(x, y): # Substracts the continuum and normalizes
    nx = len(x)
    rejt = 0.98
    ini_cont = [max(y), 0., 0.]
    coefs,_ = curve_fit(cont, x, y, p0 = ini_cont)
    #ynorm = [0. + cont(x[i], coefs[0], coefs[1], coefs[2]) for i in range(nx)]
    ynorm = cont(x, *coefs)

    for jk in range(5):
        vecx = []
        vecy = []
        for i in range(nx-1):
            if (y[i] > (ynorm[i]*rejt)) and abs(y[i]-y[i+1]) < 0.1*y[i]:
                vecx.append(x[i])
                vecy.append(y[i])

        vecx = np.array(vecx)
        vecy = np.array(vecy)
        coefs,_ = curve_fit(cont, vecx, vecy, p0 = ini_cont)
        #ynorm = [0. + cont(x[i], coefs[0], coefs[1], coefs[2]) for i in range(nx)]
        ynorm = cont(x, *coefs)

        del vecx, vecy

    #ynorm = [y[i]/cont(x[i], coefs[0], coefs[1], coefs[2]) for i in range(nx)]
    #ynorm = [(ynorm[i] - 1.) for i in range(nx)]

    ynorm = y/(cont(x, *coefs))
    ynorm = ynorm - 1.0

    #del vecx, vecy

    return np.array(ynorm)

    del ynorm

########################################################
########################################################
########################################################

def ccf(nombre, x_range, data_range2):
    make_plot = True

    c = 299792.458

    x1 = np.arange(x_range[0]-200., x_range[0], x_range[1]-x_range[0])
    x2 = np.arange(x_range[-1], x_range[-1]+200, x_range[-1]-x_range[-2])
    xtem = np.hstack([x1, x_range, x2])

    lines1, lines2, flux_l = np.genfromtxt('./binary_masks/G2.mas', dtype=None, unpack=True)
    ilines = np.where((lines1>x_range[0]) & (lines2<x_range[-1]))[0]
    lines1_new = lines1[ilines]
    lines2_new = lines2[ilines]
    flux_l_new = flux_l[ilines]

    nlines = len(flux_l_new)

    ftem = np.zeros(xtem.size)

    for i in range(len(flux_l_new)):
        indices = np.where((xtem >= lines1_new[i]) & (xtem <= lines2_new[i]))[0]
        if indices.size>0:
            ftem[indices] = flux_l_new[i]
        del indices

    rv_temp, cc = pyasl.crosscorrRV(x_range, data_range2, xtem, ftem, -200., 200., 0.2)

    del lines1, lines2, flux_l, ilines, lines1_new, lines2_new, flux_l_new

    index_min = np.argmin(cc)
    x_min = rv_temp[index_min]
    y_min = cc[index_min]
    popt,_ = curve_fit(fit_func, rv_temp, cc, p0 = (y_min, x_min, 1., 0.))
    final_rv = popt[1]
    y_fit = fit_func(rv_temp, *popt)

    if make_plot:
        fontname = 'Courier New'
        matplotlib.rcParams.update({'font.family': fontname, 'font.weight': 'medium'})
        ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=10, weight='medium', stretch='normal')

        if not os.path.exists('./Spectra/plots_ccf'):
            os.makedirs('./Spectra/plots_ccf')
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        lines = [6021.8, 6024.06, 6027.06]


        i_range_plot = np.where((x_range >= 6020) & (x_range <= 6030))[0]
        x_range_plot = x_range[i_range_plot]
        data_range_plot = data_range2[i_range_plot]

        ax1.plot(x_range_plot, data_range_plot, color = 'darkgrey')
        for i in range(len(lines)):
            ax1.axvline(x = lines[i], color = 'red')

        xnew = x_range_plot/(1. + final_rv/c)
        ax1.plot(xnew, data_range_plot, color = 'black')

        ax2.plot(rv_temp, cc, color = 'black')
        ax2.plot(rv_temp, y_fit, color = 'red')
        ax2.set_xlim(popt[1] - 5.*popt[2], popt[1] + 5.*popt[2])
        sy, ey = ax2.get_ylim()
        ax2.set_ylim(sy, 1.0)
        sx, ex = ax2.get_xlim()
        delta_x = (ex-sx)
        sy, ey = ax2.get_ylim()
        delta_y = (ey-sy)
        ax2.text(sx + 2.*delta_x/5., ey - delta_y/5., 'RV = ' + str(round(final_rv, 2)) + ' km/s', fontname = fontname)
        ax1.set_xlabel(r'$\lambda$ ($\AA$)', fontsize = 'medium', fontname = fontname)
        ax2.set_xlabel(r'Velocity (km/s)', fontsize = 'medium', fontname = fontname)

        for ax in (ax1, ax2):
            [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
            for label in ax.get_yticklabels():
                label.set_fontproperties(ticks_font)
            for label in ax.get_xticklabels():
                label.set_fontproperties(ticks_font)

        repeats = collections.Counter(nombre)['/']
        for p in range(repeats):
            indice = nombre.index('/')
            new_name = nombre[indice + 1:]
            nombre = new_name
        plt.tight_layout()
        #fig.savefig('./Spectra/plots_ccf/' + nombre + '.ps')
        fig.savefig('./Spectra/plots_ccf/' + nombre + '_ccf.pdf')
        plt.close('all')

        del fig, ax1, ax2, sy, ey, delta_y, repeats, \
            i_range_plot, x_range_plot, data_range_plot, sx, ex, delta_x

    del x1, x2, xtem, ftem, rv_temp, cc, y_fit, popt

    return final_rv

########################################################
########################################################
########################################################

def plot_lines(x, data, lines, name_lines, nombre, savefigure = True, prev_fig = None, figclear = False):
    space = 5.
    fontname = 'Courier New'

    matplotlib.rcParams.update({'font.family': fontname, 'font.weight': 'medium'})
    ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=10, weight='medium', stretch='normal')


    if prev_fig == None:
        #fig, ax = plt.subplots()
        fig = plt.figure()
    else:
        fig = prev_fig
    if figclear == True:
        fig.clear
    for l in range(len(lines)):
        i_lineas = np.where((x > (lines[l] - space)) & (x < (lines[l] + space)))[0]
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

        ax.plot(x_lineas, data_lineas, color = 'black')
        ax.axvline(lines[l], color = 'red')
        ax.set_title(name_lines[l], fontname = fontname)

        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)

        del i_lineas, x_lineas, data_lineas

    plt.draw()
    plt.tight_layout()

    if savefigure == True:
        if not os.path.exists('./Spectra/plots_fit'):
            os.makedirs('./Spectra/plots_fit')
        repeats = collections.Counter(nombre)['/']
        for p in range(repeats):
            indice = nombre.index('/')
            new_name = nombre[indice + 1:]
            nombre = new_name
        #fig.savefig('./Spectra/plots_fit/' + nombre + '.ps')
        fig.savefig('./Spectra/plots_fit/' + nombre + '_fit_lines.pdf')
    else:
        plt.show()
    return fig

########################################################
########################################################
########################################################

def restframe(archivo):
    nombre = archivo[:-5]

    hdulist = pyfits.open(archivo)
    header = hdulist[0].header
    x, data = pyasl.read1dFitsSpec(archivo)

    check_corr = True

    lineas = [6021.8, 6024.06, 6027.06, 6562.81]
    name_lineas = ['6021.8', '6024.06', '6027.06', 'H-alpha']

    i_range = np.where((x>=5500) & (x <= 6050))[0]
    x_range = x[i_range]
    data_range = data[i_range]

    if len(data_range) == 0:
        print '\t\tWavelength range to correct to restframe is not present. Cannot create the _res.fits file.'

    else:
        if np.mean(data_range) != 0.:

            data_range2 = continuum_det(x_range, data_range)

            rv = ccf(nombre, x_range, data_range2)

            del data_range2

        else: rv = 0.

        c = 299792.458

        delta = header['CDELT1']
        x0 = header['CRVAL1']
        if (rv/c) != (-1.):
            delta = delta/(1. + rv/c)
            x0 = x0/(1. + rv/c)

        #x2 = [x[i]/(1. + rv/c) for i in range(len(x))]
        #x2 = np.array(x2)
        x2 = x/(1. + rv/c)
        z = rv/c

        print '\t\trv = %f km/s' % rv

        ##################################################
        # Check the restframe correction by plotting the three lines
        # used in the CCF, plus the h-alpha line.

        check_corr = False

        if check_corr == True:
            fig = plot_lines(x2, data, lineas, name_lineas, nombre, savefigure = False)
            while True:
                check = raw_input('\t\tDo you want to shift the lines? (if yes write "yes", else press enter) ')
                if check != 'yes':
                    z2 = 0.
                    x = x2
                    break
                if check == 'yes':
                    print 'Which line you want to use?'
                    for l in range(len(lineas)):
                        print '(%d) % f' % (l, lineas[l])
                    luse = raw_input('')
                    luse = int(luse)
                    for l in range(len(lineas)):
                        if luse == l:
                            wave = raw_input('\t\tCurrent %s line wavelength? ' % (name_lineas[l]))
                            wave = float(wave)
                            z2 = (wave - lineas[l])/lineas[l]
                            break
                    x3 = [x2[i]/(1. + z2) for i in range(len(x2))]
                    x2 = np.array(x3)
                    fig = plot_lines(x2, data, lineas, name_lineas, nombre, savefigure = False)#,prev_fig=fig,figclear=True)

                    if z2 != (-1.):
                        delta = delta/(1. + z2)
                        x0 = x0/(1. + z2)
                z = z + z2 + z*z2
            plt.close('all')

        #################################################
        # Saves the plot of the spectrum corrected for restframe
        # in the three lines used by the CCF, plus the h-alpha line.

        plot_fit = True

        if plot_fit == True:
            fig = plot_lines(x2, data, lineas, name_lineas, nombre)
            plt.close('all')

        #################################################

        #print '\t\tFinal value to use is z = %f' % z

        header['CRPIX1'] = 1.
        header['CRVAL1'] = x2[0]
        header['CDELT1'] = float(delta)

        pyfits.writeto(nombre + '_res.fits', data = data, header = header[:20], clobber = True)

        hdulist.close()

        del x, data, i_range, x_range, data_range, x2, header, hdulist, delta, x0, rv

if __name__ == '__main__':
    restframe('./Spectra/HD142_harps.fits')
