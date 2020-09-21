from __future__ import division
from builtins import map
from builtins import str
from builtins import range
import os
import re
import pickle
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from past.utils import old_div
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from PyAstronomy import pyasl
from astropy.io import fits, ascii
#from .interpol_function import interpol
from AtmosInterpol import interpol
from uncertainties import unumpy, ufloat, umath

plt.style.use(['seaborn-muted'])
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


#******************************************************************************
#******************************************************************************

class Vsini:
    """
    spec_window is the spectral analysis window, a 1x2 numpy array
    gauss is the instrumental broadening parameter
    v_macro is the macroturbulence velocity
    line_file is the name of the file containing the chosen lines
    line is which of the lines on the previous file to work on
    SN is the signal-to-noise ratio of the spectrum

    # x_vel is the velocity shift to be applied to the x-axis
    # x_wl is the wavelengths shift to be applied to the x-axis
    # y_add is the additive shift to be applied to the spectrum
    # y_mult is the multiplicative shift to be applied to the spectrum
    # perf_radius is the number of points around the line center where
        # to evaluate the performance of the synthetic spectrum
    # bwing_w is the weight to be applied to the blue side of the line
        # when evaluating the performance
    # rwing_w is the weight to be applied to the red side of the line
        # when evaluating the performance
    # center_w is the weight to be applied to the line center when
        # evaluating the performance
    # Maximum number of points around the performance radius that are
        # allowed to be a bad fit (1 S/N sigma lower than observed signal)
        # If this limit is exceeded, the variable badfit_status will return
        # True after running find()
        # For high precision spectrum, set this to a very low number
    """
    def __init__(self, spec_window, gauss, v_macro, line_file, line, SN,\
                 **kwargs):

        self.name = kwargs.get('star_name', 'Unnamed star')
        self.vshift = kwargs.get('x_vel', 0.0)
        self.xshift = kwargs.get('x_wl', 0.0)
        self.yadd = kwargs.get('y_add', 0.0)
        self.ymult = kwargs.get('y_mult', 1.0)
        self.radius = kwargs.get('perf_radius', 10)
        self.bwing_w = kwargs.get('bwing_w', 3.0)
        self.rwing_w = kwargs.get('rwing_w', 5.0)
        self.center_w = kwargs.get('center_w', 25.0)
        self.badfit_tol = kwargs.get('badfit_tol', 10)

        self.c = 2.998E18
        self.am = arr_manage()
        self.spec = spec_window
        self.gauss = gauss
        self.v_m = v_macro
        self.lines = np.loadtxt(line_file, skiprows=1, usecols=(0, 1))
        try:
            self.Z = self.lines[line, 1]
            self.line_center = self.lines[line, 0]
        except IndexError:
            self.Z = self.lines[1]
            self.line_center = self.lines[0]
        self.spec_sigma = 1./SN

        self.data = np.loadtxt('./Spectra/%s_%d.dat' % (self.name, line))
        self.data_new = self.data
        self.line_number = line

        # Other attributes that will be properly assigned in other functions
        self.data_target = []
        self.center_index = 0
        self.ci0 = 0
        self.ci1 = 0
        self.MOOG = None
        self.check = None
        self.pts = 15
        self.pace = np.array([2.0, 2.0])
        self.a_guess = np.array([-0.100, 0.100])
        self.v_guess = np.array([0.5, 25.0])
        self.min_i = 3
        self.max_i = 21
        self.limits = np.array([0.01, 0.001])
        self.plot = True
        self.v_low_limit = 0.5
        self.save = None
        self.silent = True
        self.best_a = np.nan
        self.best_v = np.nan
        self.it = 0
        self.finish = False
        self.badfit_status = False
        self.it2 = 0
        self.best_v_antes = np.nan
        self.best_a_antes = np.nan
        self.v_grid = []
        self.S = []
        self.S_v = []
        self.yfit_v = []
        self.intern_u = 0
        self.best_v_ind = 0
        self.a_grid = []
        self.S_a = []
        self.yfit_a = []
        self.best_a_ind = 0
        self.go_v = 0
        self.go_a = 0
        self.v_change = None
        self.a_change = None
        self.v_width = np.nan
        self.a_width = np.nan



    def perf_new(self, v, a, mode='vsini'):
        """
        The performance function: first it creates the params.txt file, then runs
        moog in silent mode, interpolates the generated model to the points of
        the observed spectrum, and then simply calculates the sum of squared
        differences, weighted by the inverse of the observed spectrum to the power
        of alpha.
        """

        data_old = np.copy(self.data)
        data_n = np.copy(self.data)

        self.data_new[:, 0] = data_n[:, 0] + self.xshift - data_n[:, 0] * \
            (old_div(self.c, (self.vshift*1E13 + self.c)) - 1.0)
        self.data_new[:, 1] = data_n[:, 1] * self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0], self.spec[1],
                                                self.data_new)
        self.center_index = self.am.find_index(self.line_center, self.data_target[:, 0])
        self.ci0 = self.center_index - self.radius
        self.ci1 = self.center_index + self.radius+1

        if 2.*self.radius > len(self.data_target[:, 1]):
            self.radius = int(np.floor(old_div(len(self.data_target[:, 0]), 2)) - 1)
            self.ci0 = self.center_index - self.radius
            self.ci1 = self.center_index + self.radius+1

        if self.ci1 > len(self.data_target[:, 0]):
            resto = int(np.ceil((self.ci1 - len(self.data_target[:, 0]))))
            self.radius -= resto
            self.ci0 = self.center_index - self.radius
            self.ci1 = self.center_index + self.radius+1

        if self.ci0 < 0:
            self.radius -= (self.radius - self.center_index)
            self.ci0 = self.center_index - self.radius
            self.ci1 = self.center_index + self.radius+1


        if mode == 'vsini':
            S = np.inf * np.ones(v.size)

            self.MOOG.abunds = a
            self.MOOG = self.MOOG.change_vsini(v)

            for k, vsini in enumerate(v):
                model_v = self.MOOG.model_vsini[str(vsini)]
                if ~all(np.isnan(model_v.T[1])):
                    model_interp = np.interp(self.data_target[self.ci0:self.ci1, 0],\
                                             model_v.T[0], model_v.T[1])
                    w = np.zeros(2 * self.radius + 1, float)
                    if self.ci1 > len(self.data_target[:, 0]):
                        w = np.zeros(2 * self.radius, float)
                    w[:self.radius-3] = self.bwing_w
                    w[self.radius+4:] = self.rwing_w
                    w[self.radius-3:self.radius+4] = self.center_w
                    S[k] = np.sum(w * (self.data_target[self.ci0:self.ci1, 1] - \
                            model_interp)**2.) / np.sum(w)

                    del model_interp, w
                del model_v

        else:
            S = np.inf * np.ones(a.size)

            self.MOOG.vsini = v
            self.MOOG.model_ab = {}

            for k, val in enumerate(a):
                self.MOOG = self.MOOG.change_ab(val)

                model_a = self.MOOG.model_ab[str(val)]
                if ~all(np.isnan(model_a[:, 1])):
                    model_interp = np.interp(self.data_target[self.ci0:self.ci1, 0],\
                                             model_a[:, 0], model_a[:, 1])


                    w = np.zeros(2 * self.radius + 1, float)
                    if self.ci1 > len(self.data_target[:, 0]):
                        w = np.zeros(2 * self.radius, float)
                    w[:self.radius-2] = self.bwing_w
                    w[self.radius+3:] = self.rwing_w
                    w[self.radius-2:self.radius+3] = self.center_w

                    S[k] = np.sum(w * (self.data_target[self.ci0:self.ci1, 1] - \
                            model_interp)**2) / np.sum(w)

                    del model_interp, w

                del model_a

        self.data = data_old
        del data_old, data_n, self.data_target
        return S


    def perf(self, p):
        data_old = np.copy(self.data)
        data_n = np.copy(self.data)

        self.data_new[:, 0] = data_n[:, 0] + self.xshift - data_n[:, 0] * \
            (old_div(self.c, (self.vshift*1E13 + self.c)) - 1.0)
        self.data_new[:, 1] = data_n[:, 1] * self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0], self.spec[1],
                                                self.data_new)

        # Running MOOGSILENT
        self.MOOG.vsini = p[0]
        self.MOOG.abunds = p[1]
        self.MOOG = self.MOOG.run()

        # Evaluating the performance in a radius around the center of the line
        self.center_index = self.am.find_index(self.line_center,
                                               self.data_target[:, 0])
        self.ci0 = self.center_index - self.radius
        self.ci1 = self.center_index + self.radius+1

        if 2.*self.radius > len(self.data_target[:, 0]):
            self.radius = int(np.floor(old_div(len(self.data_target[:, 0]), 2)) - 1)
            self.ci0 = self.center_index - self.radius
            self.ci1 = self.center_index + self.radius+1

        model_interp = np.interp(self.data_target[self.ci0:self.ci1, 0],
                                 self.MOOG.model[:, 0],
                                 self.MOOG.model[:, 1])


        # Checking the fit on line wings
        self.check = self.data_target[self.ci0:self.ci1, 1] - model_interp
        self.check = len(np.where(self.check > 1.*self.spec_sigma)[0])

        # Creating the weights vector
        w = np.zeros(2 * self.radius + 1, float)
        if self.ci1 > len(self.data_target[:, 0]):
            w = np.zeros(2 * self.radius, float)
        w[:self.radius] = self.bwing_w
        w[self.radius+1:] = self.rwing_w
        w[self.radius] = self.center_w

        S = old_div(np.sum(w * (self.data_target[self.ci0:self.ci1, 1] - \
                   model_interp)**2), np.sum(w))

        self.data = data_old

        del data_old, data_n, model_interp, w

        return S


    def find(self, **kwargs):
        """
        -N:  Number of points to try for each iteration
        -pace:  Narrowing factor when going to the next iteration
             pace[0] = narrowing factor for vsini
             pace[1] = narrowing factor for abundance
        -a_guess: Initial guess range for abundance. It has to be a numpy array of
             length = 2
        -v_guess: Initial guess range for vsini. It has to be a numpy array of
             length = 2
        -min_i: Minimum number of iterations
        -max_i: Maximum number of iterations
        -limits: Convergence limits: a numpy array with length 2, corresponding to the
             limits of vsini and abundance, respectively
        -plot: Plot the spectral line fit at the end?
        -v_low_limit: Lower limit of estimation of vsini
        -save: Set 'save' to a filename with an extension (e.g. png, eps)
             Overrides 'plot' to False
        """
        self.pts = kwargs.get('N', 15)
        self.pace = kwargs.get('pace', np.array([2.0, 2.0]))
        self.a_guess = kwargs.get('a_guess', np.array([-0.100, 0.100]))
        self.v_guess = kwargs.get('v_guess', np.array([0.5, 25.0]))
        self.min_i = kwargs.get('min_i', 3)
        self.max_i = kwargs.get('max_i', 21)
        self.limits = kwargs.get('limits', np.array([0.01, 0.001]))
        self.plot = kwargs.get('plot', True)
        self.v_low_limit = kwargs.get('v_low_limit', 0.5)
        self.save = kwargs.get('save', None)

        if 'save' in kwargs:
            self.plot = False
        self.silent = False

        self.best_a = np.mean(self.a_guess)
        self.best_v = np.mean(self.v_guess)
        self.it = 1
        self.finish = False
        self.badfit_status = False

        MOOG = Driver(synth_interval=self.spec,\
                      abunds=np.array([[self.Z, self.best_a],]),\
                      obs_wl=self.data[:, 0], obs_flux=self.data[:, 1],\
                      gauss=self.gauss, macro_v=self.v_m,\
                      star_name=self.name, plot=self.plot,\
                      savefig=self.save,\
                      y_shift_add=self.yadd,\
                      y_shift_mult=self.ymult,\
                      wl_shift=self.xshift,\
                      line_number=self.line_number)

        self.MOOG = MOOG
        self.it2 = [0, 0]

        while ~self.finish and self.it < self.max_i and self.v_guess[1] < 100.:

            self.MOOG.it = self.it

            self.best_v_antes = self.best_v
            self.best_a_antes = self.best_a

            # Evaluating vsini
            self.v_grid = np.linspace(self.v_guess[0], self.v_guess[1], self.pts)
            self.S = []
            self.S = self.perf_new(self.v_grid, self.best_a, mode='vsini')
            self.S_v = self.S

            tck = UnivariateSpline(self.v_grid, self.S, k=4, s=0.0)#, s = 0.05)
            yfit = tck.__call__(self.v_grid)
            self.yfit_v = yfit
            self.intern_u = [False, False]

            try:

                z = tck.derivative().roots()
                tck2 = tck.derivative(n=2)
                z2 = tck2.__call__(z)

                i_s = np.where(z2 > 0.)[0]
                if i_s.size == 1:
                    self.best_v = z[i_s[0]]
                    self.intern_u[0] = True
                    self.it2[0] += 1
                    best_l = np.searchsorted(self.v_grid, self.best_v)
                    best_u = best_l + 1
                    try:
                        dif_l = self.best_v - self.v_grid[best_l]
                        dif_u = self.v_grid[best_u] - self.best_v
                        if dif_l <= dif_u:
                            self.best_v_ind = best_l
                        else:
                            self.best_v_ind = best_u
                    except IndexError:
                        self.best_v_ind = best_l
                else:
                    self.best_v_ind = np.where(self.S == min(self.S))[0][0]
                    self.best_v = self.v_grid[self.best_v_ind]

                del z, tck2, z2

            except ValueError:

                self.best_v_ind = np.where(self.S == min(self.S))[0][0]
                self.best_v = self.v_grid[self.best_v_ind]
            del tck, yfit

            # Evaluating abundance
            self.a_grid = np.linspace(self.a_guess[0], self.a_guess[1], self.pts)
            self.a_grid = self.a_grid[np.argsort(self.a_grid)]
            self.S = []
            self.S = self.perf_new(self.best_v, self.a_grid, mode='abundance')
            self.S_a = self.S
            tck = UnivariateSpline(self.a_grid, self.S, k=4, s=0.1)
            yfit = tck.__call__(self.a_grid)
            z = tck.derivative()

            self.yfit_a = yfit

            try:
                z = tck.derivative().roots()
                tck2 = tck.derivative(n=2)
                z2 = tck2.__call__(z)

                i_s = np.where(z2 > 0.)[0]
                if i_s.size == 1:
                    self.best_a = z[i_s[0]]
                    self.intern_u[1] = True
                    self.it2[1] += 1
                    best_l = np.searchsorted(self.a_grid, self.best_a)
                    best_u = best_l + 1
                    try:
                        dif_l = self.best_a - self.a_grid[best_l]
                        dif_u = self.a_grid[best_u] - self.best_a
                        if dif_l <= dif_u:
                            self.best_a_ind = best_l
                        else:
                            self.best_a_ind = best_u
                    except IndexError:
                        self.best_a_ind = best_l
                else:
                    self.best_a_ind = np.where(self.S == min(self.S))[0][0]
                    self.best_a = self.a_grid[self.best_a_ind]

                del z, tck2, z2

            except ValueError:
                self.best_a_ind = np.where(self.S == min(self.S))[0][0]
                self.best_a = self.a_grid[self.best_a_ind]

            del tck, yfit

            self.go_v = True
            self.go_a = True

            # Checking if the best values are too near the edges of the guess
            if self.best_v_ind == 0 or self.best_v_ind == (self.pts-1) or \
                    self.best_v_ind == 1 or self.best_v_ind == (self.pts-2):
                self.go_v = False
            elif self.best_a_ind == 0 or self.best_a_ind == (self.pts-1) or \
                    self.best_a_ind == 1 or self.best_a_ind == (self.pts-2):
                self.go_a = False

            # Calculating changes
            self.v_change = np.abs(self.best_v-np.mean(self.v_guess))
            self.a_change = np.abs(self.best_a-np.mean(self.a_guess))
            if ~self.silent:
                if (self.it > self.min_i) and self.go_v and self.go_a\
                    and (min(self.S) <= self.spec_sigma)\
                    and(np.abs(self.best_a - self.best_a_antes) <= 0.01)\
                    and (np.abs(self.best_v - self.best_v_antes) <= 0.01):
                    self.finish = True
                    break
                elif (self.it > self.min_i) and self.go_v and self.go_a and min(self.S) < 1e-4:
                    self.finish = True
                    break
                elif self.it > self.min_i and self.intern_u[0]\
                        and self.intern_u[1] and self.go_v and \
                        self.go_a and self.it2[0] >= 2 and self.it2[1] >= 2:
                    self.finish = True
                    break
            else:
                if self.it > self.min_i and self.intern_u[0]\
                                        and self.intern_u[1]:
                    self.finish = True
                    break

            # Setting the new guess. If one of the best values are too near the
            # edges of the previous guess, it will not narrow its new guess range.
            self.v_width = self.v_guess[1]-self.v_guess[0]
            self.a_width = self.a_guess[1]-self.a_guess[0]

            if self.go_v:
                self.v_guess = np.array([self.best_v-\
                    old_div(self.v_width, self.pace[0]), self.best_v+\
                        old_div(self.v_width, self.pace[0])])
            else:
                self.v_guess = np.array([self.best_v-old_div(self.v_width, 2),\
                    self.best_v+old_div(self.v_width, 2)])

            if self.go_a:
                self.a_guess = np.array([self.best_a-\
                    old_div(self.a_width, self.pace[1]), self.best_a+\
                        old_div(self.a_width, self.pace[1])])
            if np.abs(self.a_guess[1] - self.a_guess[0]) < 0.05:
                self.a_guess = np.array([self.best_a - 0.025,\
                    self.best_a + 0.025])

            # Checking if the v_guess contains vsini lower than v_low_limit.
            # If True, it will add a value to the array so that the lower limit
            # is equal to the v_low_limit
            if self.v_guess[0] < self.v_low_limit and ~self.silent:
                self.v_guess[0] += self.v_low_limit-self.v_guess[0]

            if self.a_guess[0] < -3.0 and ~self.silent:
                self.a_guess[0] = -3.0#+ (-3.0 - self.a_guess[0])

            self.it += 1

        # Finalizing the routine
        self.S = self.perf(np.array([self.best_v, self.best_a]))

        # Trigger bad fit warning
        if self.check > self.badfit_tol:
            self.badfit_status = True

        del MOOG

        return self

#******************************************************************************
#******************************************************************************

class Driver:
    """
    The MOOG driver object.

    Parameters
    ----------
    synth_interval : sequence
        The synthesis wavelength interval lower and upper limits in angstrons.
        Example: (6000, 6100).
    abunds : ``numpy.array``
        The atomic number (first column) and the abundance (second column) of
        the elements to be synthetisized.
        Example: numpy.array([[26, -0.05], [32, 0.01]])
    step: float, optional
        The wavelength step size of the synthesis. Default is 0.1.
    opac: float, optional
        Wavelength point to consider opacity contributions from neighboring
        transitions. Default is 2.0.
    wl_shift: float, optional
        Wavelength shift to be applied to the observed spectrum. Default is 0.
    v_shift: float, optional
        Doppler velocity shift to be applied to the observed spectrum. Default
        is 0.
    y_shift_add: float, optional
        Additive shift to be applied to the observed spectrum. Default is 0.
    y_shift_mult: float, optional
        Multiplicative factor to be applied to the observed spectrum. Default
        is 1.0 (no modification).
    gauss: float, optional
        Value of the 1 sigma dispersion of the Gaussian smoothing to be applied
        to the synthetic spectrum. Default is 0.
    lorentz: float, optional
        Default is 0.
    eps: float, optional
        Limb darkening coefficient. Default is 0.6.
    macro_v: float, optional
        Macroturbulence velocity of the star. Default is 0.
    vsini: float, optional
        The projected rotational velocity of the star. Default is 0.
    linelist_in: str, optional
        Name of the line list input file. Default is 'lines.dat'.
    observed_in: str, optional
        Name of the input file containing the observed spectrum. Default is
        'spectrum.dat'.
    atmosphere: int, optional
        Default is 1.
    molecules: int, optional
        Default is 1.
    trudamp: int, optional
        Default is 1.
    lines: int, optional
        Default is 1.
    flux: int, optional
        Default is 0.
    damping: int, optional
        Default is 0.
    star_name: str, optional
        Self-explanatory. Default is 'Unnamed star'.
    """

    def __init__(self, synth_interval, abunds, obs_wl, obs_flux, step=0.01, opac=2.0,
                 wl_shift=0.0, v_shift=0.0, y_shift_add=0.0, y_shift_mult=1.0,
                 gauss=0.0, lorentz=0.0, eps=0.6, macro_v=0.0, vsini=0.0,
                 observed_in='spectrum.dat',
                 atmosphere=1, molecules=1, trudamp=1, lines=1, flux=0,
                 damping=0, star_name='Unnamed star', plot=True, savefig=False,
                 line_number=0):

        self.name = star_name
        self.plot_switch = plot
        self.savefig = savefig
        # Output files
        self.standard_out = './output/%s_l.out' % self.name
        self.summary_out = './output/%s_li.out' % self.name
        self.smoothed_out = './output/%s_s.out' % self.name
        self.smoothed_out_new = './output/%s_sn.out' % self.name

        # Input files
        self.model_in = './atm_models/%s_v.atm' % self.name
        self.lines_in = './MOOG_linelist/lines.%s_v.txt' % self.name
        self.observed_in = './Spectra/%s_%d.dat' % (self.name, line_number)

        # Output files
        self.standard_out_moog = './output/%s_l.out' % self.name
        self.summary_out_moog = './output/%s_li.out' % self.name
        self.smoothed_out_moog = './output/%s_s.out' % self.name
        self.smoothed_out_new_moog = './output/%s_sn.out' % self.name

        # Input files
        self.model_in_moog = './atm_models/%s_v.atm' % self.name
        self.lines_in_moog = './MOOG_linelist/lines.%s_v.txt' % self.name
        self.observed_in_moog = './Spectra/%s_%d.dat' % (self.name, line_number)

        self.lines_ab = np.loadtxt(self.lines_in, usecols=(0,), skiprows=1)

        # Synthesis parameters
        self.syn_start = synth_interval[0]
        self.syn_end = synth_interval[1]
        self.wl_start = synth_interval[0]
        self.wl_end = synth_interval[1]
        self.step = step
        self.opac = opac
        self.wl_shift = wl_shift
        if int(v_shift) != 0:
            raise NotImplementedError('Doppler shift in the observed spectrum'
                                      'is not implemented yet.')
        self.v_shift = v_shift
        self.y_shift_add = y_shift_add
        self.y_shift_mult = y_shift_mult
        self.gauss = gauss
        self.lorentz = lorentz
        self.dark = eps
        self.macro_v = macro_v
        self.vsini = vsini
        self.N, self.n_cols = np.shape(abunds)
        assert(self.n_cols == 2), 'Number of columns in `abunds` must be 2.'
        if self.N == 1:
            self.Z = abunds[0][0]
            self.abunds = abunds[0][1]
        elif self.N > 1:
            self.Z = abunds[:, 0]
            self.abunds = abunds[:, 1]

        # MOOG synth options
        self.atm = atmosphere
        self.mol = molecules
        self.tru = trudamp
        self.lin = lines
        self.flu = flux
        self.dam = damping

        # Reading the observed spectrum
        if isinstance(observed_in, str):
            self.obs_wl = obs_wl + wl_shift
            self.obs_flux = obs_flux * y_shift_mult + y_shift_add
        elif observed_in is None:
            self.observed_in = observed_in
        else:
            raise TypeError('observed_in must be ``str`` or ``None``.')

        self.data = np.array([self.obs_wl, self.obs_flux]).T
        self.it = 0
        self.c = 2.998E5  # km/s
        self.model_ab = {}
        self.model_vsini = {}
        self.model = []
        self.index = 0
        self.start_index = 0
        self.end_index = 0

    def create_batch(self):
        """
        Writes the MOOG driver file batch.par
        """
        with open('./MOOGFEB2017/%s_synth.par' % self.name, 'w') as f:
            f.truncate()
            f.write("synth\n")
            f.write("standard_out  '%s'\n" % self.standard_out_moog)
            f.write("summary_out  '%s'\n" % self.summary_out_moog)
            f.write("smoothed_out  '%s'\n" % self.smoothed_out_moog)
            f.write("model_in  '%s'\n" % self.model_in_moog)
            f.write("lines_in  '%s'\n" % self.lines_in_moog)
            f.write("observed_in  '%s'\n" % self.observed_in_moog)
            f.write("atmosphere  %i\n" % self.atm)
            f.write("molecules  %i\n" % self.mol)
            f.write("trudamp  %i\n" % self.tru)
            f.write("lines    %i\n" % self.lin)
            f.write("flux/int  %i\n" % self.flu)
            f.write("damping    %i\n" % self.dam)
            f.write("freeform  0\n")
            f.write("plot    3\n")
            f.write("abundances  %i 1\n" % self.N)
            if self.N > 1:
                for k in range(self.N):
                    f.write("   %i %f\n" % (self.Z[k], self.abunds[k]))
            else:
                f.write("   %i %f\n" % (self.Z, self.abunds))
            f.write("isotopes   0  1\n")
            f.write("synlimits\n")
            f.write(" %.2f %.2f %.2f %.1f\n" % (self.syn_start, self.syn_end,
                                                self.step, self.opac))
            f.write("obspectrum  5\n")
            f.write("plotpars  1\n")
            f.write(" %.2f %.2f 0.05 1.05\n" % (self.wl_start, self.wl_end))
            f.write(" %.4f  %.4f  %.3f  %.3f\n" % (self.v_shift, self.wl_shift,
                                                   self.y_shift_add,
                                                   self.y_shift_mult))
            f.write(" gm  %.3f  0.0  %.1f  %.2f  %.1f" % (self.gauss,
                                                          self.dark,
                                                          self.macro_v,
                                                          self.lorentz))
        del f

    def change_vsini(self, grid_v):

        self.create_batch()
        os.system('MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/%s_synth.par\n\nEOF' % self.name)
        try:
            model = np.loadtxt(self.smoothed_out, skiprows=2)
            synth_wl = model[:, 0]
            synth_flux = model[:, 1]
        except:
            model = []
            synth_wl = np.arange(self.syn_start, self.syn_end, self.step)
            synth_flux = np.zeros(synth_wl.size)

        inonan = np.where(np.isfinite(synth_flux))[0]
        if len(np.unique(synth_flux[inonan])) > 1:

            synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

        self.model_vsini = {}

        for vsini in grid_v:
            self.vsini = vsini
            if self.vsini > 0.0:
                conv_flux = pyasl.rotBroad(synth_wl, synth_flux, self.dark, self.vsini)
                self.model_vsini[str(vsini)] = np.array([synth_wl, \
                                               conv_flux/max(conv_flux)]).T
                del conv_flux

            else:
                self.model_vsini[str(vsini)] = np.array([synth_wl,
                                                         np.nan*np.ones(synth_wl.size)]).T

        del model, synth_wl, synth_flux
        return self

    def change_ab(self, a):
        self.abunds = a

        self.create_batch()
        os.system('MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/%s_synth.par\n\nEOF' % self.name)
        try:
            model = np.loadtxt(self.smoothed_out, skiprows=2)
            synth_wl = model[:, 0]
            synth_flux = model[:, 1]
        except:
            model = []
            synth_wl = np.arange(self.syn_start, self.syn_end, self.step)
            synth_flux = np.zeros(synth_wl.size)

        if self.vsini > 0.0:
            inonan = np.where(np.isfinite(synth_flux))[0]
            if len(np.unique(synth_flux[inonan])) > 1:
                synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                    self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

            conv_flux = pyasl.rotBroad(synth_wl, synth_flux, self.dark, self.vsini)

            self.model_ab[str(a)] = np.array([synth_wl, old_div(conv_flux, max(conv_flux))]).T
            del conv_flux

        else:
            self.model_ab[str(a)] = np.array([synth_wl, np.nan * np.ones(synth_wl.size)]).T

        del synth_wl, synth_flux, model
        return self

    def run(self):
        """
        Used to run MOOG silent.
        """
        self.create_batch()

        os.system('MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/%s_synth.par\n\nEOF' % self.name)

        try:
            self.model = np.loadtxt(self.smoothed_out, skiprows=2)
            synth_wl = self.model[:, 0]
            synth_flux = self.model[:, 1]
        except:
            self.model = []
            synth_wl = np.arange(self.syn_start, self.syn_end, self.step)
            synth_flux = np.zeros(synth_wl.size)

        if self.vsini > 0.0:
            inonan = np.where(np.isfinite(synth_flux))[0]
            if len(np.unique(synth_flux[inonan])) > 1:
                synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                    self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

            conv_flux = pyasl.rotBroad(synth_wl, synth_flux, self.dark, self.vsini)

            self.model = np.array([synth_wl, old_div(conv_flux, max(conv_flux))]).T

            del conv_flux
        del synth_wl, synth_flux

        return self

    def rot_prof(self, vz):
        """
        This function creates a rotational profile based on Gray (2005).

        Parameters
        ----------
        vz : ``numpy.array``
            The Doppler velocities from the spectral line center.

        Returns
        -------
        profile : ``numpy.array``
            The rotational profile.
        """

        n = len(vz)
        profile = np.zeros(n, float)
        m = np.abs(vz) < self.vsini
        profile[m] = old_div((2.*(1.-self.dark)*(1.-(old_div(vz[m], self.vsini)) ** 2.)
                              ** 0.5 + 0.5 * np.pi * self.dark *
                              (1. - (old_div(vz[m], self.vsini)) ** 2.)), \
                     (np.pi * self.vsini * (1. - self.dark / 3.)))

        del m
        return profile

    @staticmethod
    def smart_cut(wl, flux, obs_wl, obs_flux):
        """
        smart_cut() is used to prepare the synthetic spectrum for a convolution
        with the rotational profile.
        """
        ind0 = np.where(flux == min(flux))[0][0]
        n = len(wl)
        if ind0 < old_div((n - 1), 2):
            if (ind0 + 1) % 2 == 0:
                wl = wl[1:2 * ind0]
                flux = flux[1:2 * ind0]
                obs_flux = obs_flux[1:2 * ind0]
                obs_wl = obs_wl[1:2 * ind0]
            else:
                wl = wl[0:2 * ind0 + 1]
                flux = flux[0:2 * ind0 + 1]
                obs_flux = obs_flux[0:2 * ind0 + 1]
                obs_wl = obs_wl[0:2 * ind0 + 1]
        elif ind0 > old_div((n - 1), 2):
            if (ind0 + 1) % 2 == 0:
                wl = wl[2*(ind0 - old_div((n - 1), 2)) + 1:-1]
                flux = flux[2 * (ind0 - old_div((n - 1), 2)) + 1:-1]
                obs_flux = obs_flux[2 * (ind0 - old_div((n - 1), 2)) + 1:-1]
                obs_wl = obs_wl[2 * (ind0 - old_div((n - 1), 2)) + 1:-1]
            else:
                wl = wl[2 * (ind0 - old_div((n - 1), 2)):]
                flux = flux[2 * (ind0 - old_div((n - 1), 2)):]
                obs_flux = obs_flux[2 * (ind0 - old_div((n - 1), 2)):]
                obs_wl = obs_wl[2 * (ind0 - old_div((n - 1), 2)):]
        return wl, flux, obs_wl, obs_flux



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


class arr_manage:
    """
    Used to perform a series of specific management routines to
    numpy-arrays and data arrays and files.
    """
    def __init__(self):
        self.index = 0
        self.start_index = 0
        self.end_index = 0
    def find_index(self, target, array):
        """
        This routine finds the index of a value closest to the target in
        a numpy-array.
        """
        self.index = np.searchsorted(array, target, side='left')
        return self.index
    def x_set_limits(self, start, end, data2d):
        """
        This routine returns a section of an array given the start and end
        values. These values do not need to be the exact ones found in the
        array.
        """
        self.start_index = np.searchsorted(data2d[:, 0], start, side='left')
        self.end_index = np.searchsorted(data2d[:, 0], end, side='left')
        return data2d[self.start_index:self.end_index]



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_broadening(starname, Teff, met, logg, micro, ab_ni,\
                    err_T, err_logg, err_met, err_ni, snr, alias):
    vmac, err_vmac = calc_vmac((Teff, err_T), (logg, err_logg))
    #vmac, err_vmac = np.ones(len(vmac))*0.1, np.ones(len(vmac))*0.1
    vsini, err_vsini, weight_vsini = calc_vsini(starname, Teff, met, logg, micro, vmac,\
                                                ab_ni, err_met, err_ni, snr, alias)

    if (vsini.size == 0) or np.all(vsini == 0.0) or np.all(err_vsini == 0.0):
        i = np.where((vmac != 0.0) & (err_vmac != 0.0))[0]
        vmac_final = np.median(vmac[i])
        err_vmac_final = old_div(np.median(err_vmac[i]), np.sqrt(float(len(i))))
        vsini_final = 0.0
        err_vsini_final = 0.0

    else:

        i = np.where((vmac != 0.0) & (vsini != 0.0) & (err_vmac != 0.0) & (err_vsini != 0.0))[0]

        #vmac_final = np.median(vmac[i])
        #err_vmac_final = old_div(np.median(err_vmac[i]), np.sqrt(float(len(i))))
        vm = np.median(unumpy.uarray(vmac[i], err_vmac[i]))
        vmac_final = vm.n
        err_vmac_final = vm.s
        #vsini_final = np.average(vsini[i], weights=weight_vsini[i])
        #err_vsini_final = np.sqrt(1./(1./(np.sum(err_vsini[i]**2.) + err_vmac_final**2.)))
        v = np.median(unumpy.uarray(vsini[i], err_vsini[i])) + ufloat(0.0, err_vmac_final)
        vsini_final = v.n
        err_vsini_final = v.s
        

    del i, vmac, err_vmac, vsini, err_vsini

    if np.isnan(err_vsini_final):
        err_vsini_final = 0.0

    return vsini_final, err_vsini_final, vmac_final, err_vmac_final


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_vmac(Teff, logg):
    """
    Computes vmac
    """
    vmac_sun = np.array([3.0, 3.2, 3.1, 3.6, 2.9])
    vmac = vmac_sun - 0.00707*Teff[0] + 9.2422*10**(-7.)*Teff[0]**2. + \
                        10.0 - 1.81*(logg[0] - 4.44) - 0.05
    err_vmac = np.ones(5)*np.sqrt(0.1**2. + (9.2422*10**(-7.)*2*Teff[0]-0.00707)**2.*Teff[1]**2. +\
                        1.81**2.*logg[1]**2. + (logg[0] - 4.44)**2.*0.26**2. + 0.03**2.)
    
    if np.mean(vmac) < 0.5:
        #From Brewer et al. 2016:

        teff = ufloat(Teff[0], Teff[1])
        log_g = ufloat(logg[0], logg[1])

        if log_g.n >= 4.0:
            vmac = 2.202*umath.exp(0.0019*(teff - 5777.)) + 1.30
        elif (4.0 > log_g.n >= 3.0):
            vmac = 1.166*umath.exp(0.0028*(teff - 5777.)) + 3.30
        else:
            vmac = 4.0 + ufloat(0.0, 0.25)
            
        logging.info(vmac.n + np.zeros(5))
        logging.info(vmac.s + np.zeros(5))

        return vmac.n + np.zeros(5), vmac.s + np.zeros(5)#, err_vmac

    # E.1 from Melendez et al. 2012:
    #vmac=np.zeros(5)+(13.499-0.00707*Teff[0]+9.2422*10**(-7)*Teff[0]**2.)
    #err_vmac=np.zeros(5)+np.sqrt((2.*9.2422*10**(-7)*Teff[0]-0.00707)**2.*Teff[1]**2.)

    # E.2 from Melendez et al. 2012
    #vmac=np.zeros(5)+(3.50+(Teff[0]-5777.)/650.)
    #err_vmac=np.zeros(5)+np.sqrt((1./650.)**2.*Teff[1]**2.)

    #From Brewer et al. 2016:

    #teff = ufloat(Teff[0], Teff[1])
    #log_g = ufloat(logg[0], logg[1])

    #if log_g.n >= 4.0:
    #    vmac = 2.202*umath.exp(0.0019*(teff - 5777.)) + 1.30
    #elif (4.0 > log_g.n >= 3.0):
    #    vmac = 1.166*umath.exp(0.0028*(teff - 5777.)) + 3.30
    #else:
    #    vmac = 4.0 + ufloat(0.0, 0.25)

    #return vmac.n + np.zeros(5), vmac.s + np.zeros(5)#, err_vmac

    return vmac, err_vmac


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_vsini(starname, Teff, met, logg, micro, v_macro, ab_ni,
               err_met, err_ni, snr=100., alias='test'):

    def f_polinomial(x, a, b):
        return a*x + b

    def f_gauss(x, a, b, c):
        return a*np.exp(old_div(-(x - b)**2., (2.*c**2.)))

    def f_total(x, a, b, c, d, e):
        return f_gauss(x, a, b, c) + f_polinomial(x, d, e)

    def continuum_det2(x, y, snr): # Substracts the continuum and normalizes
        rejt = 1.-1./snr
        p = np.poly1d(np.polyfit(x, y, 2))
        ynorm = p(x)

        for _ in range(3):
            dif = np.hstack((np.abs((y[:-1]-y[1:])/y[:-1]), [1.0]))
            i = np.where(((y-ynorm*rejt) > 0) & (dif < 0.1))[0]
            vecx = x[i]
            vecy = y[i]

            p = np.poly1d(np.polyfit(vecx, vecy, 2))
            ynorm = p(x)

        yfit = p(x)
        del p, vecx, vecy

        if np.any(yfit <= 0.0):
            i_nonzero = np.where(y > (min(y) + 0.05*(max(y)-min(y))))[0]
            xmed = min(x) + (max(x)-min(x))*0.5
            i_in_line = np.where((x >= (xmed-1.0)) & (x <= (xmed+1.0)))[0]

            if np.any(np.isin(i_in_line, i_nonzero)):
                logging.info('line not visible')
                yfit = np.zeros(len(y))

            del i_nonzero, i_in_line

        return np.array(y/yfit)

    
    # Estimate error in vsini by seeing how the vsini estimation changes
    # when changing the flux data by its noise
        
    def best_S(vrot, S):
        tck = UnivariateSpline(vrot.v_grid, S, k=4, s=0.0)#, s = 0.05)
        new_grid = np.linspace(vrot.v_grid[0], vrot.v_grid[-1], 100)
        intern_u = [False, False]
        Sfit = tck(new_grid)

        try:
            z = tck.derivative().roots()
            tck2 = tck.derivative(n=2)
            z2 = tck2.__call__(z)

            i_s = np.where(z2 > 0.)[0]
            if i_s.size == 1:
                best_v = z[i_s[0]]
                intern_u[0] = True
                #self.it2[0] += 1
                best_l = np.searchsorted(vrot.v_grid, best_v)
                best_u = best_l + 1
                try:
                    dif_l = best_v - vrot.v_grid[best_l]
                    dif_u = vrot.v_grid[best_u] - best_v
                    if dif_l <= dif_u:
                        best_v_ind = best_l
                    else:
                        best_v_ind = best_u
                except IndexError:
                    best_v_ind = best_l
                Smin = S[best_v_ind]
            else:
                #best_v_ind = np.where(S == min(S))[0][0]
                #best_v = vrot.v_grid[best_v_ind]
                Sfit = tck(new_grid)
                best_v_ind = np.argmin(Sfit)
                best_v = new_grid[best_v_ind]
                Smin = Sfit[best_v_ind]

            del z, tck2, z2

        except ValueError:
            #best_v_ind = np.where(S == min(S))[0][0]
            #best_v = vrot.v_grid[best_v_ind]
            best_v_ind = np.argmin(Sfit)
            best_v = new_grid[best_v_ind]
            Smin = Sfit[best_v_ind]

        del tck, Sfit, new_grid
    
        return best_v, Smin
    

    def get_S(vrot, ynoise):
        if 2.*vrot.radius > len(ynoise[:, 1]):
            vrot.radius = int(np.floor(old_div(len(ynoise[:, 0]), 2)) - 1)
            vrot.ci0 = vrot.center_index - vrot.radius
            vrot.ci1 = vrot.center_index + vrot.radius+1

        if vrot.ci1 > len(ynoise[:, 0]):
            resto = int(np.ceil((vrot.ci1 - len(ynoise[:, 0]))))
            vrot.radius -= resto
            vrot.ci0 = vrot.center_index - vrot.radius
            vrot.ci1 = vrot.center_index + vrot.radius+1

        if vrot.ci0 < 0:
            vrot.radius -= (vrot.radius - vrot.center_index)
            vrot.ci0 = vrot.center_index - vrot.radius
            vrot.ci1 = vrot.center_index + vrot.radius+1
    
        v = vrot.v_grid
        a = vrot.best_a
        S = np.inf * np.ones(v.size)

        vrot.MOOG.abunds = a
        #self.MOOG = self.MOOG.change_vsini(v)

        for k, vsini in enumerate(v):
            model_v = vrot.MOOG.model_vsini[str(vsini)]
            if ~all(np.isnan(model_v.T[1])):
                model_interp = np.interp(ynoise[vrot.ci0:vrot.ci1, 0],\
                                         model_v.T[0], model_v.T[1])
                w = np.zeros(2 * vrot.radius + 1, float)
                if vrot.ci1 > len(ynoise[:, 0]):
                    w = np.zeros(2 * vrot.radius, float)
                w[:vrot.radius-3] = vrot.bwing_w
                w[vrot.radius+4:] = vrot.rwing_w
                w[vrot.radius-3:vrot.radius+4] = vrot.center_w
                S[k] = np.sum(w * (ynoise[vrot.ci0:vrot.ci1, 1] - \
                                   model_interp)**2.) / np.sum(w)

                del model_interp, w
            del model_v
        
        return S

    def noise_estimation(vrot, snr, N=1000):
        vsini = np.nan*np.ones(N)
        ydata = vrot.data_target
        minSnoise = np.nan*np.ones(N)
    
        noise_array = np.zeros((N,len(ydata)))
        for i in range(len(ydata)):
            noise_array.T[i] = np.random.normal(0.0, np.abs(ydata[:,1][i]/snr), N)
        
        for i in range(N):
            ynoise = np.copy(ydata)
            ynoise[:,1] = ydata[:,1]+noise_array[i]
            Snoise = get_S(vrot, ynoise)
            vsini[i], minSnoise[i] = best_S(vrot, Snoise)
        return vsini, minSnoise
    
    
    # Set lines
    # Creates dic with characteristics of lines
    lines1 = {'name' : 'FeI', 'wave' : 6027.05, 'Z' : 26, 'EP' : 4.076, 'loggf' : -1.09}
    lines2 = {'name' : 'FeI', 'wave' : 6151.62, 'Z' : 26, 'EP' : 2.176, 'loggf' : -3.30}
    lines3 = {'name' : 'FeI', 'wave' : 6165.36, 'Z' : 26, 'EP' : 4.143, 'loggf' : -1.46}
    lines4 = {'name' : 'FeI', 'wave' : 6705.10, 'Z' : 26, 'EP' : 4.607, 'loggf' : -0.98}
    lines5 = {'name' : 'NiI', 'wave' : 6767.77, 'Z' : 28, 'EP' : 1.826, 'loggf' : -2.17}
    lines_o = (lines1, lines2, lines3, lines4, lines5)

    # Create model atmosphere
    interpol(starname, Teff, logg, met, micro, alias+'_v')
    inst = starname[starname.index('_')+1:]
    lines_ab, dev_ab = calc_ab(starname, err_met, err_ni, alias)
    line_file = './MOOG_linelist/lines.%s_v.txt' % alias

    try:
        x, data = pyasl.read1dFitsSpec('./Spectra/%s_res.fits' % starname)
    except:
        hdu = fits.open('./Spectra/%s_res.fits' % starname)
        d = hdu[0].data
        x = d[0]
        data = d[1]
        hdu.close()
        del hdu, d


    new_data = np.array([x, data]).T
    resolution = None

    try:
        resolution = fits.getval('./Spectra/%s_res.fits' % starname, 'R', 0)
    except (IndexError, KeyError):
        pass

    # For each line create vsini object
    vsini_lines = np.zeros(5)
    err_vsini_lines = np.zeros(5)
    weight_lines = np.zeros(5)
    ab_keys = list(map(float, list(lines_ab.keys())))
    data_lines = {}

    available_lines = np.loadtxt(line_file, skiprows=1, usecols=(0))

    lines = []
    for l in lines_o:
        if l['wave'] in available_lines:
            lines.append(l)

    for l, li in enumerate(lines):
        data_o = new_data[:]

        w = li['wave']

        if w in ab_keys:
            info_line = {}
            info_line['vmac'] = v_macro[l]
            kwargs = {'star_name' : alias}
            j = (x > (w - 3.0)) & (x < (w + 3.0))
            #j = (x > (w - 2.0)) & (x < (w + 2.0))

            x_l = x[j]
            y_l = data[j]
            
            if x_l.size == 0 or y_l.size == 0:
                continue

            k = (x > (w - 0.5)) & (x < (w+0.5))
            #k = (x > (w - 1.0)) & (x < (w+1.0))

            i_r = int(np.floor(old_div(len(x[k]), 2)))
            kwargs['perf_radius'] = i_r

            try:
                new_y_l = continuum_det2(x_l, y_l, snr)
                popt, _ = curve_fit(f_total, x_l, new_y_l,
                                    p0=(old_div(-max(y_l), min(y_l)), w, 0.2, 0.0, max(y_l)))
                if popt[1] > (w + 0.3) or popt[1] < (w - 0.3):
                    x_shift = 0.0
                else:
                    x_shift = w - popt[1]

                new_x_l = x_l + x_shift
                del popt

            except (RuntimeError, TypeError, ValueError):
                try:
                    popt, _ = curve_fit(f_total, x[j], data[j],
                                        p0=(old_div(-max(data[k]), min(data[k])),
                                            w, 0.2, 0.0, max(data[j])))
                    y_polinomial = f_polinomial(x_l, *popt[3:])

                    if popt[1] > (w + 0.3) or popt[1] < (w - 0.3):
                        x_shift = 0.0
                    else:
                        x_shift = w - popt[1]

                    new_x_l = x_l + x_shift
                    new_y_l1 = old_div(y_l, f_total(x_l, *popt))

                    popt2, _ = curve_fit(f_polinomial, new_x_l, new_y_l1, p0=(0.0, 1.0))
                    y_polinomial2 = f_polinomial(new_x_l, *popt2)

                    new_y_l = old_div(y_l, y_polinomial/y_polinomial2)

                    del popt2, y_polinomial, y_polinomial2, new_y_l1, popt

                except RuntimeError:
                    new_x_l = x_l[:]
                    new_y_l = y_l[:]

            kwargs['badfit_tol'] = 50

            ascii.write([new_x_l, new_y_l], './Spectra/%s_%d.dat' % (alias, l),\
                        format='fixed_width_no_header', overwrite=True, delimiter='\t')
            SN = snr

            del j, x_l, y_l, new_x_l, new_y_l

            if resolution is None:
                if inst == 'harps':
                    gauss = w/115000.
                elif inst in ['feros', 'feros_o']:
                    gauss = w/48000.
                elif inst == 'uves':
                    gauss = w/110000.
                elif inst in ['hires', 'HIRES']:
                    gauss = w/67000.
                elif inst == 'coralie':
                    gauss = w/60000.
                elif inst == 'psf':
                    gauss = w/38000.
                else:
                    x_ = ascii.read('./Spectra/%s_%d.dat' % (alias, l))['col1']
                    R = min(np.mean(x_)/np.mean(x_[1:]-x_[:-1]), 150000.)      
                    gauss = w/R
                    #print(R)
                    del x_, R

            else:
                gauss = w/float(resolution)

            spec_window = np.array([w-1.0, w+1.0])

            if v_macro[l] > 0.0:
                vmacro = v_macro[l]
            else:
                vmacro = 0.1
                info_line['vmac'] = vmacro

            vrot = Vsini(spec_window, gauss, vmacro, line_file, l, SN, **kwargs)
            logging.info('Working on line %.3f, line abundance is %.3f, a_guess is [%.3f, %.3f]', \
                          w, lines_ab['%.3f' % w], \
                          (lines_ab['%.3f' % w]) - 2.5*max(0.1, dev_ab['%.3f' % w]),\
                          (lines_ab['%.3f' % w]) + 2.5*max(0.1, dev_ab['%.3f' % w]))

            kwargs2 = {'a_guess' : np.array([(lines_ab['%.3f' % w])\
                                              -2.5*max(0.1, dev_ab['%.3f' % w]),
                                             (lines_ab['%.3f' % w])\
                                              + 2.5*max(0.1, dev_ab['%.3f' % w])]),\
                       'v_guess' : np.array([0.1, 25.]),\
                       'save' : True,\
                       'N' : 30,\
                       'v_low_limit' : 0.1,\
                       'max_i' : 30}


            vrot = vrot.find(**kwargs2)
            info_line['data'] = vrot.MOOG.data
            info_line['model'] = vrot.MOOG.model
            info_line['vsini'] = vrot.best_v
            info_line['v_grid'] = vrot.v_grid
            info_line['S_v'] = vrot.S_v
            info_line['a_grid'] = vrot.a_grid
            info_line['S_a'] = vrot.S_a
            info_line['abundance'] = vrot.best_a
            info_line['yfit_v'] = vrot.yfit_v
            info_line['yfit_a'] = vrot.yfit_a

            #Compare fit with a straight line
            ilm = np.where(np.abs(vrot.MOOG.model[:, 0]-w) <= 1.0)[0]
            il = np.where(np.abs(vrot.MOOG.data[:, 0]-w) <= 1.0)[0]
            model_interp = UnivariateSpline(vrot.MOOG.model[:, 0][ilm],
                                            vrot.MOOG.model[:, 1][ilm],
                                            s=0, k=5)(vrot.MOOG.data[:, 0][il])
            S_model = np.sum((vrot.MOOG.data[:, 1][il]-model_interp)**2.)
            S_line = np.sum((vrot.MOOG.data[:, 1][il]-1.0)**2.)
            if S_line < S_model:
                vrot.badfit_status = True

            del model_interp, S_model, S_line

            info_line['badfit'] = vrot.badfit_status
            data_lines[str(w)] = info_line

            del info_line
            if ~vrot.badfit_status:
                #vsini_lines[l] = vrot.best_v
                #err_vsini_lines[l] = np.sqrt(vrot.S)
                weight_lines[l] = 1./vrot.S
                vsini_dist, minSnoise = noise_estimation(vrot, snr, N=1000)
                p = np.percentile(vsini_dist, [16, 50, 84])
                vsini_lines[l] = p[1]
                err_vsini_lines[l] = max(p[1]-p[0], p[2]-p[1])
                data_lines[str(w)]['vsini_dist'] = vsini_dist
                del vsini_dist

            del vrot
            del kwargs, spec_window, kwargs2
            os.system('rm -f ./Spectra/%s_%d.dat' % (alias, l))


        new_data = data_o
        del data_o

    f = open('./plots_broadening/%s_data_lines.pkl' % starname, 'wb')
    pickle.dump(data_lines, f)
    f.close()

    plot_paper(starname, data_lines)
    plot_dist(starname, data_lines)

    del new_data, data_lines, lines1, lines2, lines3, lines4, lines5, lines, f

    os.system('rm -f ./atm_models/%s_v.atm' % alias)
    os.system('rm -f ./output/%s_l.out' % alias)
    os.system('rm -f ./output/%s_li.out' % alias)
    os.system('rm -f ./output/%s_s.out' % alias)
    os.system('rm -f ./output/%s_sn.out' % alias)
    os.system('rm -f ./MOOG_linelist/lines.%s_v.txt' % alias)
    os.system('rm -f ./MOOGFEB2017/abfind_%s_v_2.par' % alias)
    os.system('rm -f ./MOOGFEB2017/abfind_%s_v.par' % alias)
    os.system('rm -f ./MOOGFEB2017/%s_synth.par' % alias)
    os.system('rm -f ./output/%s.dat' % alias)
    os.system('rm -f ./output/%s_o.dat' % alias)
    
    # Correct for the solar values, so that vsini_sun = 1.9 km/s
    #vsini_lines = vsini_lines - 1.28

    return vsini_lines, err_vsini_lines, weight_lines


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_ab(starname, err_met, err_ni, alias='test'):

    moog_linelist(starname, alias)

    cmd = 'cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/abfind_%s_v.par' % (alias)
    os.system(cmd)

    with open('./MOOGFEB2017/abfind_%s_v.par' % (alias), 'r') as par:
        with open('./MOOGFEB2017/abfind_%s_v_2.par' % (alias), 'w') as par_out:
            for linea in par:
                columnas = linea.strip()
                m = re.search(r'standard_out\s*(\S*).*', columnas)
                if m:
                    linea = "standard_out  './output/%s.dat'\n" % (alias)
                m = re.search(r'summary_out\s*(\S*).*', columnas)
                if m:
                    linea = "summary_out   './output/%s_o.dat'\n" % (alias)
                m = re.search(r'model_in\s*(\S*).*', columnas)
                if m:
                    linea = "model_in      './atm_models/%s_v.atm'\n" % (alias)
                m = re.search(r'lines_in\s*(\S*).*', columnas)
                if m:
                    linea = "lines_in      './MOOG_linelist/lines.%s_v.txt'\n" % (alias)
                par_out.writelines(linea)

    cmd = 'cp ./MOOGFEB2017/abfind_%s_v_2.par ./MOOGFEB2017/abfind_%s_v.par' % (alias, alias)
    os.system(cmd)

    os.system('MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/abfind_%s_v.par\n\nEOF' % alias)

    ab = {}
    dev = {}
    with open('./output/%s_o.dat' % alias) as output:
        for linea in output:
            linea = linea.strip()
            m = re.search(r'[a-z]', linea)
            if m is None:
                m = re.search(r'[\d]', linea)
                if m:
                    linea = linea.split()
                    ID = int(float(linea[1]))
                    if ID == 26:
                        ab_id = 7.50# + met
                        dev_id = err_met
                    else:
                        ab_id = 6.22# + ab_ni
                        dev_id = err_ni
                    ab[linea[0]] = float(linea[6]) - ab_id# - met
                    dev[linea[0]] = dev_id
            del m

    del cmd, par, par_out, output
    return ab, dev


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def moog_linelist(starname, alias='test'):
    linelist = np.genfromtxt('./Spectra/linelist_vsini.dat', dtype=None, skip_header=2,\
                             names=('line', 'excit', 'loggf', 'num', 'ion'))
    line = linelist['line']
    excit = linelist['excit']
    loggf = linelist['loggf']
    num = linelist['num']
    ion = linelist['ion']

    file_ew = np.genfromtxt('./EW/%s_vsini.txt' % starname, dtype=None,\
                            names=('line', 'ew', 'ew_e', 'ew_err1', 'ew_err2'))

    line_ew_b = file_ew['line']
    ew_b = file_ew['ew']
    ew_err_b = np.maximum(file_ew['ew_err1'], file_ew['ew_err2'])

    #Take only the lines that 10. <= EW <=150 and the error in the EW is lower than the EW
    ilines = np.where((ew_b >= 10.) & (ew_b <= 150.) & (old_div(ew_err_b, ew_b) <= 1.0))[0]
    line_ew = line_ew_b[ilines]
    ew = ew_b[ilines]
    ew_err = ew_err_b[ilines]

    with open('MOOG_linelist/lines.%s_v.txt' % alias, 'w') as output:
        output.write(' %s_vsini.txt\n' % starname)
        for i, l in enumerate(line_ew):
            index = np.where(line == l)[0]
            if len(index) == 1:
                index = int(index[0])
                output.write('%s%7.2f%s%4.1f%s%5.2f%s%6.3f%s%7.4f\n' %\
                                  (' '*2, line[index], ' '*4, ion[index], ' '*7, excit[index], \
                                   ' '*5, loggf[index], ' '*23, ew[i]))
    del output

    del linelist, line, excit, loggf, num, ion, file_ew, line_ew_b, ew_b, ew_err_b,\
        ilines, line_ew, ew, ew_err

#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def compute_snr(x, data, w):
    """
    Computes the S/N for spectra using a set of
    ranges where there shouldn't be any lines,
    and only continuum.
    """

    k1 = np.where(x < (w - 0.3))[0]
    k2 = np.where(x > (w + 0.3))[0]

    sn = []
    if k1.size > 0:
        sn.append(old_div(np.mean(data[k1]), np.std(data[k1])))
    if k2.size > 0:
        sn.append(old_div(np.mean(data[k2]), np.std(data[k2])))

    del k1, k2
    return np.mean(sn)


#******************************************************************************

def plot_paper(starname, data_lines):
    ticks_font = matplotlib.font_manager.FontProperties(style='normal', size=9,
                                                        weight='medium', stretch='normal')

    def _plot_grid(ax, value, w, grid, S, yfit, type_g='vsini'):
        ax.plot(grid, S, '.-', color='dimgrey')
        ax.plot(grid, yfit, color='tomato')
        ax.axvline(value, color='tomato')
        ax.locator_params(nbins=4)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
        if type_g == 'vsini':
            ax.set_xlabel(r'$v\sin i$ (km/s)', fontsize=9)
        else:
            ax.set_xlabel('abundance', fontsize=9)
        ax.set_ylabel('$S$', fontsize=9)

        _ = [i.set_linewidth(0.5) for i in ax.spines.values()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(ex - 5.*(ex-sx)/11., ey - (ey-sy)/10., \
                r'$\lambda$ = %s $\AA$' % (w),\
                style='italic', fontsize=8, backgroundcolor='white')
        del sx, sy, ex, ey
        return ax


    def _add_plot(ax, data, model, w, vsini, vmac, badfit):
        x_limits = [data[:, 0][0], data[:, 0][-1]]
        i_m1 = np.where(model[:, 0] < data[:, 0][0])[0]
        i_m2 = np.where(model[:, 0] > data[:, 0][-1])[0]
        model[:, 0][i_m1] = 1.0
        model[:, 0][i_m2] = 1.0
        del i_m1, i_m2
        if badfit:
            r = ax.patch
            r.set_facecolor('red')
            r.set(alpha=0.2)
            del r
        ax.plot(data[:, 0], data[:, 1], marker='.', ls='None', color='dimgrey')
        ax.plot(model[:, 0], model[:, 1], color='tomato')

        ax.tick_params(axis='both', which='major', labelsize=8)

        _ = [i.set_linewidth(0.5) for i in ax.spines.values()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        ax.set_xlabel(r'$\lambda$ ($\AA$)', fontsize=9)
        ax.set_ylabel(r'$F_\lambda$ $d\lambda$', fontsize=9)

        ax.set_xlim(x_limits)
        ax.locator_params(axis='x', nbins=6)
        ax.locator_params(axis='y', nbins=5)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(sx + (ex-sx)/13., sy + (ey-sy)/10., \
                r'$\lambda$ = %s $\AA$'
                '\n'
                r'$v \sin{i}$ = %.2f'
                '\n'
                r'$v_{macro}$ = %.2f' % (w, vsini, vmac),\
                style='italic', fontsize=8)
        del x_limits, sx, sy, ex, ey
        return ax

    try:
        lines = sorted(data_lines.keys())
        n = len(lines)
        cols, rows = 3, n

        fig, ax = plt.subplots(rows, cols, figsize=(cols*4, rows*2.5))

        if rows == 1:
            d = data_lines[lines[0]]
            ax[0] = _plot_grid(ax[0], d['abundance'], lines[0], d['a_grid'],
                               d['S_a'], d['yfit_a'], type_g='abundance')
            ax[1] = _plot_grid(ax[1], d['vsini'], lines[0], d['v_grid'],
                               d['S_v'], d['yfit_v'], type_g='vsini')
            ax[2] = _add_plot(ax[2], d['data'], d['model'], lines[0], d['vsini'],
                              d['vmac'], d['badfit'])
            del d
        else:
            for r in range(rows):
                d = data_lines[lines[r]]
                ax[r][0] = _plot_grid(ax[r][0], d['abundance'], lines[r], d['a_grid'],
                                      d['S_a'], d['yfit_a'], type_g='abundance')
                ax[r][1] = _plot_grid(ax[r][1], d['vsini'], lines[r], d['v_grid'],
                                      d['S_v'], d['yfit_v'], type_g='vsini')
                ax[r][2] = _add_plot(ax[r][2], d['data'], d['model'], lines[r],
                                     d['vsini'], d['vmac'], d['badfit'])
                del d

        fig.subplots_adjust(hspace=0.3, wspace=0.25, bottom=0.05, left=0.07, right=0.98, top=0.98)
        fig.savefig('./plots_broadening/%s_vsini_paper.pdf' % starname)
        plt.close('all')

        del lines, n, cols, rows, fig, ax

    except:
        pass
        
def plot_dist(starname, data_lines):
    try:
        lines = sorted(data_lines.keys())
        n = len(lines)
        fig, ax = plt.subplots(1, n, figsize=(n*3, 3))
        if n == 1:
            d = data_lines[lines[0]]
            if 'vsini_dist' in d:
                ax.hist(d['vsini_dist'], bins=40)
                p = np.percentile(d['vsini_dist'], [16, 50, 84])
                ax.axvline(p[1], label=r'%.2f $\pm$ %.2f km/s' % (p[1], max(p[1]-p[0], p[2]-p[1])), color='orange')
                ax.legend()
                ax.set_xlabel('vsini (km/s)')
            del d
        else:
            for i in range(n):
                d = data_lines[lines[i]]
                if 'vsini_dist' in d:
                    ax[i].hist(d['vsini_dist'], bins=40)
                    p = np.percentile(d['vsini_dist'], [16, 50, 84])
                    ax[i].axvline(p[1], label=r'%.2f $\pm$ %.2f km/s' % (p[1], max(p[1]-p[0], p[2]-p[1])), color='orange')
                    ax[i].legend()      
                    ax[i].set_xlabel('vsini (km/s)')
                del d
        fig.subplots_adjust(hspace=0.3, wspace=0.25, bottom=0.2, left=0.03, right=0.98, top=0.95)
        fig.savefig('./plots_broadening/%s_vsini_dist.pdf' % starname)
        plt.close('all')
    except Exception as e:
        print(e)
        pass  
