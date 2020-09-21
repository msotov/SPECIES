"""
Module to compute the Equivalent Widths (EW) of a set of lines,
for a given spectrum.
Uses regression methods to fit Gaussian-like profiles to the
absorption lines, and outputs the distribution of EWs
for every line (with a valid fit), based on the fit parameters
and uncertainties.

**Important**
The spectra has to be shifted to restframe before attempting to
compute the EWs, otherwise the lines won't be found.

Usage:

EWComputation.EW_calc(starname, wave, flux,
                      linelist='linelist.dat',
                      snr=100., makeplot=False,
                      path_plots='./EW/plots_EW',
                      from_dic=False, save_line_data=False)

Input:
    starname:          Name of the star, after which the output
                       files will be named.

    wave:              Numpy array.
                       Wavelength array in Angstroms.
                       Can be 1D or more, but has to have
                       the same shape as the flux array.

    flux:              Numpy array.
                       Flux array.
                       Can be 1D or more, but has to have
                       the same shape as the wave array.

    [Optional]

    linelist:          Name of the linelist file containing the
                       wavelengths of the lines (in angstrom) to be fitted.
                       Lines must be under the 'WL' column.

    snr:               str or numpy array.
                       Signal-to-noise (S/N) of the data. It can either
                       be just a number, or an array. If array,
                       must have the same dimensions as number of
                       orders in the wave and flux arrays.

    makeplot:          Boleean.
                       Generate or not plots with the fits for inspection.

    path_plots:        Str.
                       Path to save the plots.

    from_dic:          Boleean.
                       Use the fit parameters from a file. This file must have
                       been created in a previous run of the code,
                       using save_line_data=True.

    save_line_data:    Boleean.
                       Save the line fits to a dictionary file.

Output:
    EW file:           The output of the code will be a file
                       (stored under the EW/ directory)
                       with five columns:
                       [Wavelength, EW, EW_mean, EW_plus_error, EW_minus_error].
                       EW refers to the value first computed by the most
                           probable fit parameters.
                       EW_mean, EW_plus_error, EW_minus_error:
                           [50%, 84%-50%, 50%-16%] of the EW distribution
                           (computed using the parameters estimated for the
                           Gaussian fit, plus their uncertainties).

    EW plot:           Plot with the Gaussian fits for each line, including
                       the continuum fit. Generated only if makeplot=True.


Example:

import numpy as np
from astropy.io import fits
from EWComputation import EW_calc

hdu = fits.open('sun01_harps_res.fits')
data = hdu[0].data
header = hdu[0].header

wave = data[0]
flux = data[1]
snr = header['SNR']

EW_calc('sun01_harps', wave, flux,
        linelist='linelist.dat', snr=SN, makeplot=True)

"""


from __future__ import print_function
from __future__ import division
from builtins import range
import os
import pickle
import logging
import warnings
import sys
from dataclasses import dataclass
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import ascii
from astropy.convolution import convolve, Box1DKernel
from scipy import interpolate
from scipy.signal import convolve as scipy_convolve
import scipy.odr as ODR

warnings.simplefilter("ignore")
plt.style.use(['seaborn-muted'])
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


__author__ = "Maritza Soto (marisotov@gmail.com)"
__date__ = "2020-09-18"
__version__ = "2.0.0"


@dataclass
class gaussian:
    a: float
    m: float
    s: float

    def __call__(self, x):
        return self.a*np.exp(-(x-self.m)**2./(2.*self.s**2.))

    @property
    def params(self):
        return (self.a, self.m, self.s)

    @property
    def errors(self):
        return (self.err_a, self.err_m, self.err_s)

    def set(self, a, m, s):
        self.a = a
        self.m = m
        self.s = s

    def seterr(self, e_a, e_m, e_s):
        self.err_a = e_a
        self.err_m = e_m
        self.err_s = e_s

    def EW(self, x, cdelt):
        integral = self(x).sum()
        return -integral*cdelt*1000.0

    def EWdist(self, x, a, m, s, cdelt=None):
        n = a.size
        EW = np.ones(n)*np.nan
        if cdelt is None:
            cdelt = np.median(x[1:]-x[:-1])
        for i in range(n):
            self.set(a[i], m[i], s[i])
            EW[i] = self.EW(x, cdelt)
        inonan = np.where(~np.isnan(EW))[0]
        return EW[inonan]

    @staticmethod
    def distribution_fit(x, a, m, s):
        n = a.size
        y = np.zeros((n, x.size))
        y16 = np.zeros(x.size)
        y50 = np.zeros(x.size)
        y84 = np.zeros(x.size)

        for i in range(n):
            y[i] = a[i]*np.exp(-(x-m[i])**2./(2.*s[i]**2.))

        yT = y.T
        for i in range(x.size):
            y16[i], y50[i], y84[i] = np.percentile(yT[i], [16, 50, 84])
        del y, yT
        return y16, y50, y84


@dataclass
class Gaussian:
    ngauss: int

    def __call__(self, x, *components):
        y = np.zeros(x.size)
        for i in range(self.ngauss):
            y += gaussian(*components[3*i:3*i+3])(x)
        self.eval = y
        return y

    def add_component(self, x, *components):
        icomp = len(components)/3
        self.ngauss += 1
        y = self.eval
        for i in range(icomp):
            y += gaussian(*components[3*i:3*i+3])(x)
        self.eval = y
        return y

    @staticmethod
    def params(components):
        a = [components[i] for i in range(len(components)) if i%3 == 0]
        m = [components[i] for i in range(len(components)) if i%3 == 1]
        s = [components[i] for i in range(len(components)) if i%3 == 2]
        return np.array(a), np.array(m), np.array(s)

    @staticmethod
    def errors(components):
        a = [components[i] for i in range(len(components)) if i%3 == 0]
        m = [components[i] for i in range(len(components)) if i%3 == 1]
        s = [components[i] for i in range(len(components)) if i%3 == 2]
        return np.array(a), np.array(m), np.array(s)

    @staticmethod
    def line_closest_to(l, components):
        return np.argmin(np.abs(components-l))

    def set_components(self, components):
        self.components = [gaussian(*components[3*i:3*i+3]) for i in range(self.ngauss)]



class Spectrum:

    def __init__(self, line, wave, flux, snr):
        self.line = line
        self.wave = wave
        self.flux = flux
        if np.min(flux) < 0.0:
            self.flux = flux - np.min(flux)
        self.snr = snr
        self.rejt = 1.-1./snr
        self.is_visible = True
        self.flux_norm = np.ones(self.flux.size)*np.nan
        self.width = None
        self.dxarr = np.ones(wave.size)*np.nan

    def make_dxarr(self, coordinate_location='center'):
        dxarr = np.diff(self.wave)
        if self.wave.size <= 2:
            self.dxarr = np.ones(self.wave.size)*dxarr
        elif coordinate_location in ['left', 'center']:
            self.dxarr = np.concatenate([dxarr, dxarr[-1:]])
        elif coordinate_location in ['right']:
            self.dxarr = np.concatenate([dxarr[:1, dxarr]])

    def cdelt(self, tolerance=1e-8, approx=False):
        if np.all(np.isnan(self.dxarr)):
            self.make_dxarr()
        if approx or abs(self.dxarr.max()-self.dxarr.min())/abs(self.dxarr.min()) < tolerance:
            return self.dxarr.mean().flat[0]
        return np.mean(self.wave[1:]-self.wave[:-1])

    def check_data(self):
        wave_right = self.wave[np.where(self.wave <= self.line)[0]]
        wave_left = self.wave[np.where(self.wave >= self.line)[0]]
        if wave_right.size > 5 and wave_left.size > 5:
            del wave_right, wave_left
            return True
        del wave_right, wave_left
        return False

    @property
    def data(self):
        return [self.wave, self.flux]

    @property
    def data_norm(self):
        return [self.wave, self.flux_norm]

    @staticmethod
    def _consecutive(data, stepsize=1):
        return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

    def _find_emission_lines(self, wave, flux, flux_threshold=None):
        kernel = [1, 0, -1]
        dY = scipy_convolve(flux, kernel, 'valid')
        S = np.sign(dY)
        ddS = scipy_convolve(S, kernel, 'valid')
        candidates = np.where(dY > 0)[0] + (len(kernel)- 1)
        line_inds = sorted(set(candidates).intersection(np.where(ddS == -2)[0] + 1))
        if flux_threshold is not None:
            line_inds = np.array(line_inds)[flux[line_inds] > flux_threshold]
        line_inds_grouped = self._consecutive(line_inds, stepsize=1)
        if len(line_inds_grouped[0]) > 0:
            emission_inds = [inds[np.argmax(flux[inds])] for inds in line_inds_grouped]
        else:
            emission_inds = []
        return wave[emission_inds], emission_inds

    def _find_absorption_lines(self, flux_threshold=None, width=None):
        wave, flux = self.data_norm
        self.width = width
        kernel = [1, 0, -1]
        dY = scipy_convolve(flux, kernel, mode='valid')
        S = np.sign(dY)
        ddS = scipy_convolve(S, kernel, mode='valid')
        candidates = np.where(dY < 0)[0] + (len(kernel) - 1)
        line_inds = sorted(set(candidates).intersection(np.where(ddS == 2)[0] + 1))

        if flux_threshold is not None:
            line_inds = np.array(line_inds)[flux[line_inds] < -flux_threshold]

        # Now group them and find the max highest point.
        line_inds_grouped = self._consecutive(line_inds, stepsize=1)

        if len(line_inds_grouped[0]) > 0:
            absorption_inds = [inds[np.argmin(flux[inds])] for inds in
                               line_inds_grouped]
        else:
            absorption_inds = []
        if width is not None and absorption_inds:
            ii = np.where(np.abs(wave[absorption_inds]-self.line) <= width)[0]
            absorption_inds = np.array(absorption_inds)[ii]
            del ii
        return wave[absorption_inds], absorption_inds

    def normalize(self, order=2, detect_emission=False, passes=5):
        x, y = self.data
        p = np.poly1d(np.polyfit(x, y, order))
        ynorm = p(x)
        lines_em_final = []

        for ii in range(passes):
            dif = np.hstack((np.abs((y[:-1]-y[1:])/y[:-1]), [1.0]))
            i = np.where(((y-ynorm*self.rejt) > 0) & (dif < 0.1))[0]
            vecx = x[i]
            vecy = y[i]
            p = np.poly1d(np.polyfit(vecx, vecy, order))
            ynorm = p(x)
            if ii == 0 and detect_emission:
                lines_em, _ = self._find_emission_lines(x, y/ynorm-1.0,
                                                        flux_threshold=0.09/self.rejt**2.)
                if len(lines_em) > 0:
                    yy = vecy/p(vecx)-1.0
                    i_ex = []
                    lines_em_final = []
                    for ll in lines_em:
                        ii_ex = np.where(np.abs(vecx-ll) <= 0.25)[0].astype(int)
                        ii_ex2 = np.where((np.abs(vecx-ll) > 0.25))[0].astype(int)
                        if np.abs(np.mean(yy[ii_ex])-np.mean(yy[ii_ex2])) >= np.std(yy[ii_ex2]):
                            i_ex.append(np.where(np.abs(x-ll) <= 0.2)[0].astype(int))
                            lines_em_final.append(ll)
                    if i_ex:
                        i_ex = np.hstack(i_ex)
                        i_ex = np.unique(i_ex)
                        i_no_ab = np.array([ii_ for ii_ in range(x.size) if (ii_ in i_ex) is False])
                        x_no_ab = x[i_no_ab]
                        y_no_ab = y[i_no_ab]
                        p = np.poly1d(np.polyfit(x_no_ab, y_no_ab, order))
                        ynorm = p(x_no_ab)
                        for _ in range(passes):
                            dif = np.hstack((np.abs((y_no_ab[:-1] - y_no_ab[1:])/y_no_ab[:-1]),
                                             [1.0]))
                            i = np.where(((y_no_ab - ynorm * self.rejt) > 0) & (dif < 0.1))[0]
                            vecx = x_no_ab[i]
                            vecy = y_no_ab[i]
                            p = np.poly1d(np.polyfit(vecx, vecy, order))
                            ynorm = p(x_no_ab)
                        break

        yfit = p(x)
        del ynorm

        if np.any(yfit <= 0.0):
            i_nonzero = np.where(y > (min(y) + 0.05*(max(y)-min(y))))[0]
            xmed = min(x) + (max(x)-min(x))*0.5
            i_in_line = np.where((x >= (xmed-1.0)) & (x <= (xmed+1.0)))[0]

            if np.any(np.isin(i_in_line, i_nonzero)):
                logging.info('line not visible')
                yfit = np.zeros(len(y))
                self.is_visible = False

            del i_nonzero, i_in_line

        self.flux_norm = y/yfit-1.0

        return vecx, vecy, p


    def detect_lines(self):
        wave, flux = self.data_norm
        # Select only the regions where the spectra is larger than -1
        ivalid = np.where(flux > -1.0)[0]
        ysmooth = convolve(flux[ivalid], Box1DKernel(3))
        for n in range(3):
            gradient = np.gradient(ysmooth)
            ysmooth = convolve(gradient, Box1DKernel(3))
            if n == 1:
                tckn = interpolate.UnivariateSpline(wave[ivalid], ysmooth, k=5, s=0)

        tck = interpolate.UnivariateSpline(wave[ivalid], ysmooth, k=3, s=0)
        zeros = tck.roots()
        flux_zeros = interpolate.UnivariateSpline(wave[ivalid], flux[ivalid], k=5, s=0)(zeros)
        ifinal = np.where((tckn(zeros) > 0.0) & (flux_zeros < -0.02))[0]

        del ysmooth, gradient, tckn, tck, ivalid

        return zeros[ifinal]

    def get_exc_list(self, lines, width=0.15):
        x, y = self.data_norm
        if self.width is not None:
            ii = np.where(np.abs(x-self.line) <= self.width)[0]
            ycopy = np.copy(y)[ii]
            xcopy = np.copy(x)[ii]
            del ii
        else:
            ycopy = np.copy(y)
            xcopy = np.copy(x)
        for l in lines:
            iline = np.where(np.abs(xcopy-l) <= width)[0]
            ycopy[iline] = np.nan

        return np.where(np.isnan(ycopy))[0]


    def fit_gauss(self, guess=None, exc_list=None, lines=None):
        x, y = self.data_norm
        if self.width is not None:
            ii = np.where(np.abs(x-self.line) <= self.width)[0]
            xcopy = np.copy(x)[ii]
            ycopy = np.copy(y)[ii]
            del ii
        else:
            xcopy = np.copy(x)
            ycopy = np.copy(y)

        if lines is None:
            ncomp = 1
        else:
            ncomp = len(lines)

        if guess is None:
            print('Needs a guess!')
            return np.ones(ncomp*3, dtype=int), np.ones(ncomp*3)

        fgauss = Gaussian(ncomp)
        if exc_list is not None:
            ycopy[exc_list] = np.nan
            inonan = np.where(~np.isnan(ycopy))[0]
            xcopy = xcopy[inonan]
            ycopy = ycopy[inonan]

        ivalid = np.where(ycopy > -1.0)[0]

        def fgauss2(B, x):
            return fgauss(x, *B)

        weights = np.ones(xcopy[ivalid].size)
        ii = np.where(np.abs(xcopy[ivalid]-self.line) > 0.05)[0]
        weights[ii] = 0.9
        ii = np.where(np.abs(xcopy[ivalid]-self.line) > 0.10)[0]
        weights[ii] = 0.75
        ii = np.where(np.abs(xcopy[ivalid]-self.line) > 0.15)[0]
        weights[ii] = 0.6
        ii = np.where(np.abs(xcopy[ivalid]-self.line) > 0.20)[0]
        weights[ii] = 0.5
        ii = np.where(np.abs(xcopy[ivalid]-self.line) > 1.00)[0]
        weights[ii] = 0.2
        weights = weights/np.sum(weights)

        func = ODR.Model(fgauss2)
        mydata = ODR.Data(x=xcopy[ivalid], y=ycopy[ivalid], we=weights)
        myodr = ODR.ODR(mydata, func, beta0=guess)
        myoutput = myodr.run()
        s_rw = myoutput.beta
        del func, mydata, myodr, myoutput

        func = ODR.Model(fgauss2)
        mydata = ODR.Data(x=xcopy[ivalid], y=ycopy[ivalid])
        myodr = ODR.ODR(mydata, func, beta0=guess)
        myoutput = myodr.run()
        e_s_rw = myoutput.sd_beta

        del weights, ii, func, mydata, myodr, myoutput, xcopy, ycopy, ivalid

        return s_rw, e_s_rw


    def guess(self, lines=None):
        x, y = self.data_norm
        tck = interpolate.UnivariateSpline(x, y, s=0, k=5)

        if lines is None:
            lines = [self.line]

        g = []
        for l in lines:
            g.append(float(tck(l)))
            g.append(l)
            g.append(0.05)

        return g


def EW_for_line(l, wave_l, flux_l, rejt, snr_o):
    if snr_o < 20.:
        snr_o = 20.0
        rejt = 1.-1./snr_o
    try:
        sp = Spectrum(l, wave_l, flux_l, snr_o)
        if sp.check_data() is False:
            del sp
            print('\t\t\tLine %.2f is not visible in data. Skipping.' % l)
            logging.error('Line %.2f is not visible in data. Skipping.', l)
            return l, 0.0, 0.0, 0.0, 0.0, {}

        vecx, vecy, p = sp.normalize()
        x, y = sp.data_norm
        if np.any(np.isinf(y)):
            inoinf = np.where(~np.isinf(y))[0]
            x = x[inoinf]
            y = y[inoinf]
            del inoinf
        if np.all(np.isnan(y)) or np.all(np.isinf(y)):
            del sp, x, y, vecx, vecy
            print('\t\t\tLine %.2f has some invalid datapoints. Skipping.' % l)
            logging.error('Line %.2f has some invalid datapoints. Skipping.', l)
            return l, 0.0, 0.0, 0.0, 0.0, {}

        dicl = {}
        dicl['wave'] = wave_l
        dicl['flux'] = y
        dicl['rejt'] = rejt
        dicl['vecx'] = vecx
        dicl['vecy'] = vecy
        dicl['polinomial'] = p
        dicl['flux_no_norm'] = sp.flux

        lines, _ = sp._find_absorption_lines(flux_threshold=0.04, width=2.0)
        if len(np.where(np.abs(lines-l) <= 0.10)[0]) > 1:
            ii = np.where(np.abs(lines-l) > 0.10)[0]
            lines = lines[ii]
        if ~np.any(np.abs(lines-l) < 0.1):
            lines = np.append(lines, l)
            isort = np.argsort(lines)
            lines = lines[isort]
        guess = sp.guess(lines=lines)
        try:
            popt, pcov = sp.fit_gauss(guess=guess, lines=lines)
        except RuntimeError:
            lines_far_away = lines[np.where(np.abs(lines-l) > 1.0)[0]]
            lines_close = lines[np.where(np.abs(lines-l) <= 1.0)[0]]
            if len(lines_close) == 0:
                lines = [l]
                guess = sp.guess(lines=lines)
                exc_list = sp.get_exc_list(lines_far_away)
                popt, pcov = sp.fit_gauss(guess=guess, lines=lines, exc_list=exc_list)
            else:
                guess = sp.guess(lines=lines_close)
                exc_list = sp.get_exc_list(lines_far_away)
                popt, pcov = sp.fit_gauss(guess=guess, lines=lines_close, exc_list=exc_list)
                lines = lines_close

        g = Gaussian(len(lines))
        a, m, s = g.params(popt)
        err_a, err_m, err_s = g.errors(pcov)
        iline = g.line_closest_to(l, m)
        logging.info(r'Obtained a= %.2f +/- %.2f, mu=%.2f +/- %.2f, s=%.2f +/- %.2f',
                     a[iline], err_a[iline], m[iline], err_m[iline], s[iline], err_s[iline])
        if (np.abs(m[iline]-l) > 0.075) or (a[iline] > 0) or (a[iline] < -1) or\
            (s[iline] > 0.10) or\
            np.any(np.array([err_a[iline], err_m[iline], err_s[iline]]) > 0.12):
            if len(lines) == 1 and lines[0] == l:
                print('\t\t\tLine %.2f: Incorrect fit' % l)
                logging.error('Line %.2f: Incorrect fit', l)
                return l, 0.0, 0.0, 0.0, 0.0, {}
            logging.info('Doing the fit again')
            lines_far_away = lines[np.where(np.abs(lines-l) > 0.15)[0]]
            lines = [l]
            guess = sp.guess(lines=lines)
            exc_list = sp.get_exc_list(lines_far_away)
            popt, pcov = sp.fit_gauss(guess=guess, lines=lines, exc_list=exc_list)
            g = Gaussian(len(lines))
            a, m, s = g.params(popt)
            err_a, err_m, err_s = g.errors(pcov)
            iline = g.line_closest_to(l, m)
            logging.info(r'Obtained a= %.2f +/- %.2f, mu=%.2f +/- %.2f, s=%.2f +/- %.2f',
                         a[iline], err_a[iline], m[iline], err_m[iline], s[iline], err_s[iline])
            if (np.abs(m[iline]-l) > 0.075) or (a[iline] > 0) or (a[iline] < -1) or\
                (s[iline] > 0.10) or\
                np.any(np.array([err_a[iline], err_m[iline], err_s[iline]]) > 0.15):
                print('\t\t\tLine %.2f: Incorrect fit' % l)
                return l, 0.0, 0.0, 0.0, 0.0, {}

        g.set_components(popt)
        cdelt = sp.cdelt(approx=True)
        final_eqws = g.components[iline].EW(x, cdelt)

        print(
            '\t\t\tLine {:6.2f}: a={: 3.2f}, mu={:5.1f}, s={:4.3f}, EW={:6.2f}, rejt={:4.3f}, '
            'snr={:5.2f}'.format(
                l,
                a[iline], m[iline], s[iline],
                final_eqws,
                rejt,
                snr_o))
        logging.info(
            '{:6.2f}: a={: 3.2f}, mu={:5.1f}, s={:4.3f}, EW={:6.2f}, rejt={:4.3f}, '
            'snr={:5.2f}'.format(
                l,
                a[iline], m[iline], s[iline],
                final_eqws,
                rejt,
                snr_o))

        dist_a = np.random.normal(a[iline], err_a[iline], 1000)
        dist_m = np.random.normal(m[iline], err_m[iline], 1000)
        dist_s = np.random.normal(s[iline], err_s[iline], 1000)

        EWdist = g.components[iline].EWdist(x, dist_a, dist_m, dist_s, cdelt=cdelt)

        dicl['EW'] = final_eqws
        dicl['snr'] = snr_o
        dicl['EW_dist'] = EWdist
        dicl['parvalues'] = [a[iline], m[iline], s[iline]]
        dicl['parerrors'] = [err_a[iline], err_m[iline], err_s[iline]]
        dicl['popt'] = popt
        dicl['pcov'] = pcov
        dicl['lines'] = lines
        dicl['iline'] = iline
        dicl['dist_params'] = [dist_a, dist_m, dist_s]
        dicl['full_model'] = g(x, *popt)

        p = np.percentile(EWdist, [16, 50, 84])
        final_eqws_mean = p[1]
        final_eqws_err1 = p[2]-p[1]
        final_eqws_err2 = p[1]-p[0]

        del sp, g, lines, guess, popt, pcov, dist_a, dist_m, dist_s, EWdist, p
        del a, m, s, err_a, err_m, err_s, x, y

        return l, final_eqws, final_eqws_mean, final_eqws_err1, final_eqws_err2, dicl


    except Exception as e:
        logging.error(e)
        logging.error('Problem with line %.2f. Ignoring.', l)
        exc_type, _, exc_tb = sys.exc_info()
        logging.error('line %d: %s', exc_tb.tb_lineno, e)
        logging.error(exc_type)
        return l, 0.0, 0.0, 0.0, 0.0, {}


def EW_calc(starname, wave, flux, linelist='linelist.dat', snr=100., makeplot=False,
            path_plots='./EW/plots_EW',
            from_dic=False, save_line_data=False):

    print('\t\tCreating EW file for %s' % starname)
    if (from_dic is False) or os.path.isfile('EW/%s_line_data.pkl' % starname) is False:

        norders = wave.ndim

        # Select the lines
        lines = ascii.read('Spectra/%s' % linelist, comment='-')['WL']

        # Arrays to store the results for each line
        final_eqws = np.zeros(len(lines))
        final_eqws_mean = np.zeros(len(lines))
        final_eqws_err1 = np.zeros(len(lines))
        final_eqws_err2 = np.zeros(len(lines))

        dic = {}
        space = 3.0
        if norders == 1:
            logging.info('Only one order detected')

        for i, l in enumerate(lines):
            # Select a window of 3 AA around the line
            if norders is 1:
                iline = np.where((wave > l-space) & (wave < l+space))[0]
                wave_l = wave[iline]
                flux_l = flux[iline]
                rejt = 1.-(1./snr)
                snr_o = snr
                if len(wave_l) == 0:
                    logging.warning('Line not found in spectrum.')
                    continue

            else:
                iline_coords = np.where((wave > l-space) & (wave < l+space))
                # Check is the line is present in more than one order
                if np.unique(iline_coords[0]).size > 1:
                    logging.info('Line present in more than one order')
                    # Check that the order's limits are within the line
                    ix = np.unique(iline_coords[0])
                    order_limits = np.zeros(ix.size)
                    for x_i, xorder in enumerate(ix):
                        wave_nonzero = wave[xorder][np.nonzero(wave[xorder])]
                        if wave_nonzero[0] < (l-space):
                            order_limits[x_i] += 1
                        if wave_nonzero[-1] > (l+space):
                            order_limits[x_i] += 1


                    if len(np.where(order_limits == max(order_limits))[0]) > 1:
                        # Choose the redder order
                        order_mean = np.zeros(ix.size)
                        for x_i, xorder in enumerate(ix):
                            wave_nonzero = wave[xorder][np.nonzero(wave[xorder])]
                            order_mean[x_i] = np.mean(wave_nonzero)

                        iorder = ix[np.where(order_mean == max(order_mean))[0]]

                    else:
                        iorder = ix[np.where(order_limits == max(order_limits))[0]]


                    x_order = iline_coords[0][np.where(iline_coords[0] == iorder)[0]]
                    y_order = iline_coords[1][np.where(iline_coords[0] == iorder)[0]]

                    if hasattr(snr, "__len__"):
                        rejt = 1.-(1./snr[iorder])
                        snr_o = snr[iorder]
                    else:
                        rejt = 1.-(1./snr)
                        snr_o = snr

                    wave_l = wave[x_order, y_order]
                    flux_l = flux[x_order, y_order]

                else:
                    logging.info('Line in just one order')
                    wave_l = wave[iline_coords]
                    flux_l = flux[iline_coords]

                    if wave_l.size == 0:
                        logging.warning('Line not found in spectrum')
                        continue

                    if hasattr(snr, "__len__"):
                        iorder = np.unique(iline_coords[0])
                        rejt = 1.-(1./snr[iorder])
                        snr_o = snr[iorder]
                    else:
                        rejt = 1.-(1./snr)
                        snr_o = snr

                del iline_coords

            if hasattr(snr_o, "__len__"):
                snr_o = snr_o[0]
            if hasattr(rejt, "__len__"):
                rejt = rejt[0]

            if snr_o > 500.:
                snr_o = 500.0
                rejt = 1.-(1./snr_o)


            l, final_eqw, final_eqw_mean, final_eqw_err1, final_eqw_err2, dicl = EW_for_line(l,
                                                                                             wave_l,
                                                                                             flux_l,
                                                                                             rejt,
                                                                                             snr_o)
            if final_eqw != 0.0:
                final_eqws[i] = final_eqw
                final_eqws_mean[i] = final_eqw_mean
                final_eqws_err1[i] = final_eqw_err1
                final_eqws_err2[i] = final_eqw_err2
                dic['%f' % l] = dicl

            del dicl, wave_l, flux_l


        ascii.write([lines, final_eqws, final_eqws_mean, final_eqws_err1, final_eqws_err2],\
                    'EW/%s.txt' % starname, format='fixed_width_no_header', delimiter=' ',\
                    overwrite=True)

        if save_line_data:
            fl = open('EW/%s_line_data.pkl' % starname, 'wb')
            pickle.dump(dic, fl)
            fl.close()

    else:
        print('\t\tReading EW information from EW/%s_line_data_test.pkl' % starname)
        dic = pickle.load(open('EW/%s_line_data_test.pkl' % starname, 'rb'))

    if makeplot:
        plot_lines(starname, dic, path_plots)

    del dic



def plot_lines(starname, dic, path_plots):

    class PlotGrid:

        def __init__(self, ncols, nrows, path_plots='.', islines=True):
            self.ncols = ncols
            self.nrows = nrows
            if islines:
                self.fig = plt.figure(figsize=(6*6, 2*nrows))
            else:
                self.fig = plt.figure(figsize=(3*6, 2*nrows))
            self.grid = gridspec.GridSpec(nrows, ncols, figure=self.fig)
            self.coordinates = np.arange(ncols*nrows).reshape((nrows, ncols))
            self.islines = islines
            self.path_plots = path_plots

        def __call__(self, i):
            irow, icol = np.where(self.coordinates == i)
            if self.islines:
                gs0 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=self.grid[irow, icol],
                                                       hspace=0)
                axi1 = plt.Subplot(self.fig, gs0[0, :])
                self.fig.add_subplot(axi1)
                axi = plt.Subplot(self.fig, gs0[1, :])
                self.fig.add_subplot(axi)

                return axi1, axi, irow, icol

            ax_err = self.fig.add_subplot(self.nrows, self.ncols, i+1)
            return ax_err, irow, icol

        def finish_plot(self, starname):
            if self.islines:
                self.fig.subplots_adjust(bottom=1./(self.nrows*4), top=1-1./(self.nrows*6),
                                         left=0.015, right=0.99, wspace=0.1)
                self.fig.savefig(os.path.join(os.path.relpath(self.path_plots, '.'),
                                              '%s.pdf' % starname), format='pdf')
            else:
                self.fig.subplots_adjust(bottom=0.02, top=0.99, left=0.05, right=0.98, wspace=0.2)
                self.fig.savefig(os.path.join(os.path.relpath(self.path_plots, '.'),
                                              '%s_error.pdf' % starname), format='pdf')

    try:
        ncols = 6
        lines = list(dic.keys())
        nrows = int(np.ceil(len(lines)/6.))

        LinesPlot = PlotGrid(ncols, nrows, path_plots=path_plots, islines=True)
        ErrorPlot = PlotGrid(ncols, nrows, path_plots=path_plots, islines=False)

        for i, l in enumerate(lines):
            d = dic[l]

            ####### Plot the fits to the lines ########

            wave_l = d['wave']
            flux_l = d['flux'] + 1.0
            values = d['parvalues']
            errors = d['parerrors']

            ax_cont, ax, irow, icol = LinesPlot(i)

            # Continuum
            vecx = d['vecx']
            vecy = d['vecy']
            flux_no_norm = d['flux_no_norm']
            mean_flux = np.mean(flux_no_norm)
            ax_cont.plot(wave_l, flux_no_norm - mean_flux + 1.0, lw=0.5, color='black')
            p = d['polinomial']
            xcont = np.linspace(wave_l[0], wave_l[-1], 100)
            ycont = p(xcont)
            ax_cont.plot(xcont, ycont - mean_flux + 1.0, lw=0.5, color='steelblue')
            ax_cont.plot(vecx, vecy - mean_flux + 1.0, ls='None', marker='x', color='orangered', markersize=1.2)
            del vecx, vecy, p, xcont, ycont

            # Line

            ax.plot(wave_l, flux_l, lw=0.4, color='black')
            ax.axvline(float(l), color='red', ls=':', lw=0.5)
            ax.plot(wave_l, d['full_model']+1.0, color='green', lw=0.5)
            xfit = np.hstack((np.linspace(min(wave_l), float(l)-1.0, 50),
                              np.linspace(float(l)-1.0, float(l)+1.0, 200),
                              np.linspace(float(l)+1.0, max(wave_l), 50)))
            yfit = values[0]*np.exp(-(xfit-values[1])**2/ (2.0*values[2]**2))
            ax.plot(xfit, yfit+1.0, color='steelblue', lw=0.5)
            yfit16, yfit50, yfit84 = plot_fit_distribution(xfit, values, errors)
            ax.fill_between(xfit, yfit16+1.0, y2=yfit50+1.0,
                            alpha=0.4, edgecolor='white', lw=0.1, color='orangered')
            ax.fill_between(xfit, yfit50+1.0, y2=yfit84+1.0,
                            alpha=0.4, edgecolor='white', lw=0.1, color='orangered')
            ax.plot(xfit, yfit50+1.0, color='red', lw=0.3)

            del xfit, yfit, yfit16, yfit50, yfit84

            ax.axhline(1.0, ls=':', color='gray')

            _ = [
                ax.axvline(
                    l_,
                    color='gray',
                    ls='-',
                    alpha=0.5,
                    lw=0.5) for l_ in d['lines']]


            # Set limits, label sizes, number of ticks, and text to be added

            sx, ex = min(wave_l), max(wave_l)
            sy, ey = ax_cont.get_ylim()
            ax_cont.set_xlim(sx, ex)
            ax_cont.xaxis.set_ticks(np.arange(np.ceil(sx), ex, 1))
            if ex-sx > 6.0:
                ax_cont.xaxis.set_ticks(np.arange(np.ceil(sx), ex, 2))
            ax_cont.yaxis.set_ticks(np.linspace(min(flux_no_norm-mean_flux+1.0),
                                                max(flux_no_norm-mean_flux+1.0), 5)[1:])
            ax_cont.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
            ax_cont.tick_params(labelsize='x-small', labelbottom=False)
            if icol == 0:
                ax_cont.set_ylabel('Flux', fontsize='small')
            del flux_no_norm

            sy, ey = ax.get_ylim()
            sy = max(sy, 0.0)
            ew = d['EW']/1000.
            ew16, ew50, ew84 = np.percentile(d['EW_dist'], [16, 50, 84])/1000.
            midpt = float(l)
            ax.fill_between([midpt-ew50/2.0, midpt+ew50/2.0], [0, 0],
                            [1.0, 1.0], color='g', alpha=0.3, edgecolor='white', lw=0.1)
            ax.text(sx+(ex-sx)/20., sy+(ey-sy)/7.,
                    r'EW$\,=\,%.1f^{+%.1f}_{-%.1f}\,$mA' % (ew50*1000., (ew84-ew50)*1000.,
                                                           (ew50-ew16)*1000.),
                    fontsize='x-small', bbox=dict(facecolor='white', alpha=0.8,
                                                  edgecolor='None', linestyle='None'))

            ax.text(ex-(ex-sx)/4, sy+(ey-sy)/6.,
                    r'$A\,=\,%.3f\,\pm\,%.3f$''\n'r'$\mu\,=\,%.1f\,\pm\,%.3f$''\n'
                    r'$\sigma\,=\,%.3f\,\pm\,%.3f$' %
                    (values[0], errors[0], values[1], errors[1], values[2], errors[2]),
                    fontsize='x-small', bbox=dict(facecolor='white', alpha=0.8,
                                                  edgecolor='None', linestyle='None'))

            ax.set_ylim(max(sy, 0.0), min(1.30, ey))
            ax.set_xlim(sx, ex)

            ax.xaxis.set_ticks(np.arange(np.ceil(sx), ex, 1))
            ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))

            sy, ey = ax.get_ylim()
            ax.yaxis.set_ticks(np.linspace(sy, ey, 6))
            ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
            ax.set_ylim(sy, ey)

            ax.tick_params(labelsize='x-small')
            if icol == 0:
                ax.set_ylabel('Norm. Flux', fontsize='small')
            else:
                ax.set_ylabel(' ')
            if irow == nrows-1:
                ax.set_xlabel('Wavelength', fontsize='small')
            elif (irow == nrows-2) and (icol >= (ncols-((ncols*nrows)%len(lines)))):
                ax.set_xlabel('Wavelength', fontsize='small')
            else:
                ax.set_xlabel(' ')

            ####### Plot the error distributions ########

            ax_err, irow, icol = ErrorPlot(i)

            inonan = np.where(~np.isnan(d['EW_dist']))[0]
            ax_err.hist(d['EW_dist'][inonan], bins=30, histtype='stepfilled', color='steelblue')
            ax_err.axvline(ew16*1000., color='orange')
            ax_err.axvline(ew50*1000., color='orange')
            ax_err.axvline(ew84*1000., color='orange')
            ax_err.tick_params(labelsize='small')

            start_x, end_x = ax_err.get_xlim()
            ax_err.xaxis.set_ticks(np.linspace(start_x, end_x, 5))
            ax_err.axvline(ew*1000., color='orangered')

            del wave_l, flux_l, d

        LinesPlot.finish_plot(starname)
        ErrorPlot.finish_plot(starname)

        plt.close('all')
        del LinesPlot, ErrorPlot

    except Exception as e:
        logging.error('Error while plotting:')
        _, _, exc_tb = sys.exc_info()
        logging.error('line %d: %s', exc_tb.tb_lineno, e)
        logging.error('Skipping...')


def plot_fit_distribution(wave, parvalues, parerrors):
    n = 1000
    y = np.zeros((n, wave.size))
    y16 = np.zeros(wave.size)
    y50 = np.zeros(wave.size)
    y84 = np.zeros(wave.size)
    a = np.random.normal(parvalues[0], parerrors[0], n)
    m = np.random.normal(parvalues[1], parerrors[1], n)
    s = np.random.normal(parvalues[2], parerrors[2], n)

    for i in range(n):
        y[i] = a[i]*np.exp(-(wave-m[i])**2./(2.*s[i]**2.))

    yT = y.T
    for i in range(wave.size):
        y16[i], y50[i], y84[i] = np.percentile(yT[i], [16, 50, 84])
    del y, a, m, s, yT
    return y16, y50, y84
