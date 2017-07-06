import os, re, time
from astropy.io import fits, ascii
from interpol_function import interpol
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['font.family'] = 'serif'
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import convolve1d
from PyAstronomy import pyasl

plt.style.use(['classic'])
fontname = 'Courier New'
matplotlib.rcParams.update({'font.family': fontname, 'font.weight': 'medium'})
ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=5, weight='medium', stretch='normal')


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

class Vsini(object):


    # spec_window is the spectral analysis window, a 1x2 numpy array
    # gauss is the instrumental broadening parameter
    # v_macro is the macroturbulence velocity
    # line_file is the name of the file containing the chosen lines
    # line is which of the lines on the previous file to work on
    # SN is the signal-to-noise ratio of the spectrum
    def __init__(self, spec_window, gauss, v_macro, line_file, line, SN,\
                 **kwargs):

        # Writes star's name at the plot
        if ('star_name' in kwargs):
            self.name = kwargs['star_name']
        else:
            self.name = 'Unnamed star'

        # x_vel is the velocity shift to be applied to the x-axis
        if ('x_vel' in kwargs):
            self.vshift = kwargs['x_vel']
        else:
            self.vshift = 0.0

        # x_wl is the wavelengths shift to be applied to the x-axis
        if ('x_wl' in kwargs):
            self.xshift = kwargs['x_wl']
        else:
            self.xshift = 0.0

        # y_add is the additive shift to be applied to the spectrum
        if ('y_add' in kwargs):
            self.yadd = kwargs['y_add']
        else:
            self.yadd = 0.0

        # y_mult is the multiplicative shift to be applied to the spectrum
        if ('y_mult' in kwargs):
            self.ymult = kwargs['y_mult']
        else:
            self.ymult = 1.0

        # perf_radius is the number of points around the line center where
        # to evaluate the performance of the synthetic spectrum
        if ('perf_radius' in kwargs):
            self.radius = kwargs['perf_radius']
        else:
            self.radius = 5

        # bwing_w is the weight to be applied to the blue side of the line
        # when evaluating the performance
        if ('bwing_w' in kwargs):
            self.bwing_w = kwargs['bwing_w']
        else:
            self.bwing_w = 1.0

        # bwing_w is the weight to be applied to the red side of the line
        # when evaluating the performance
        if ('rwing_w' in kwargs):
            self.rwing_w = kwargs['rwing_w']
        else:
            self.rwing_w = 2.5

        # center_w is the weight to be applied to the line center when
        # evaluating the performance
        if ('center_w' in kwargs):
            self.center_w = kwargs['center_w']
        else:
            self.center_w = 10.0

        # Maximum number of points around the performance radius that are
        # allowed to be a bad fit (1 S/N sigma lower than observed signal)
        # If this limit is exceeded, the variable badfit_status will return
        # True after running find()
        # For high precision spectrum, set this to a very low number
        if ('badfit_tol' in kwargs):
            self.badfit_tol = kwargs['badfit_tol']
        else:
            self.badfit_tol = 10

        self.c = 2.998E18
        self.am = arr_manage()
        self.spec = spec_window
        self.gauss = gauss
        self.v_m = v_macro
        self.lines = np.loadtxt(line_file, skiprows=1, usecols=(0,1))
        self.Z = self.lines[line,1]
        self.line_center = self.lines[line,0]
        self.spec_sigma = 1./SN

        self.data = np.loadtxt('./Spectra/%s_%d.dat' % (self.name, line))
        self.data_new = self.data
        self.line_number = line


    '''
    The performance function: first it creates the params.txt file, then runs
    moog in silent mode, interpolates the generated model to the points of
    the observed spectrum, and then simply calculates the sum of squared
    differences, weighted by the inverse of the observed spectrum to the power
    of alpha.
    '''
    def perf_new(self, v, a, mode = 'vsini'):

        data_old = np.copy(self.data)
        data_n = np.copy(self.data)

        self.data_new[:,0] = data_n[:,0] + self.xshift - data_n[:,0] * \
            (self.c / (self.vshift*1E13 + self.c) - 1.0)
        self.data_new[:,1] = data_n[:,1] * self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0], self.spec[1],
                                                self.data_new)
        self.center_index = self.am.find_index(self.line_center, self.data_target[:,0])
        self.ci0 = self.center_index - self.radius
        self.ci1 = self.center_index + self.radius+1

        if 2.*self.radius > len(self.data_target[:,0]):
            self.radius = int(np.floor(len(self.data_target[:,0])/2) - 1)
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


            for k,vsini in enumerate(v):
                model_v = self.MOOG.model_vsini[str(vsini)]
                if all(np.isnan(model_v[:,1])) == False:
                    model_interp = np.interp(self.data_target[self.ci0:self.ci1,0],\
                                             model_v[:,0], model_v[:,1])


                    w = np.zeros(2 * self.radius + 1, float)
                    w[:self.radius] = self.bwing_w
                    w[self.radius+1:] = self.rwing_w
                    w[self.radius] = self.center_w
                    S[k] = np.sum(w * (self.data_target[self.ci0:self.ci1,1] - \
                            model_interp)**2.) / np.sum(w)


                    del model_interp, w
                del model_v

        else:
            S = np.inf * np.ones(a.size)

            self.MOOG.vsini = v
            self.MOOG.model_ab = {}

            for k,a in enumerate(a):
                self.MOOG = self.MOOG.change_ab(a)

                model_a = self.MOOG.model_ab[str(a)]
                if all(np.isnan(model_a[:,1])) == False:
                    model_interp = np.interp(self.data_target[self.ci0:self.ci1,0],\
                                             model_a[:,0], model_a[:,1])


                    w = np.zeros(2 * self.radius + 1, float)
                    w[:self.radius] = self.bwing_w
                    w[self.radius+1:] = self.rwing_w
                    w[self.radius] = self.center_w

                    S[k] = np.sum(w * (self.data_target[self.ci0:self.ci1,1] - \
                            model_interp)**2) / np.sum(w)


                    del model_interp, w

                del model_a


        self.data = data_old

        del data_old, data_n, self.data_target

        return S


    def perf(self,p):

        data_old = np.copy(self.data)
        data_n = np.copy(self.data)

        self.data_new[:,0] = data_n[:,0] + self.xshift - data_n[:,0] * \
            (self.c / (self.vshift*1E13 + self.c) - 1.0)
        self.data_new[:,1] = data_n[:,1] * self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0], self.spec[1],
                                                self.data_new)

        # Running MOOGSILENT
        self.MOOG.vsini = p[0]
        self.MOOG.abunds = p[1]
        self.MOOG = self.MOOG.run()

        # Evaluating the performance in a radius around the center of the line
        self.center_index = self.am.find_index(self.line_center,
                                               self.data_target[:,0])
        self.ci0 = self.center_index - self.radius
        self.ci1 = self.center_index + self.radius+1

        if 2.*self.radius > len(self.data_target[:,0]):
            self.radius = int(np.floor(len(self.data_target[:,0])/2) - 1)
            self.ci0 = self.center_index - self.radius
            self.ci1 = self.center_index + self.radius+1

        model_interp = np.interp(self.data_target[self.ci0:self.ci1,0],
                                      self.MOOG.model[:,0],
                                      self.MOOG.model[:,1])


        # Checking the fit on line wings
        self.check = self.data_target[self.ci0:self.ci1,1] - model_interp
        self.check = len(np.where(self.check > 1.*self.spec_sigma)[0])

        # Creating the weights vector
        w = np.zeros(2 * self.radius + 1, float)
        w[:self.radius] = self.bwing_w
        w[self.radius+1:] = self.rwing_w
        w[self.radius] = self.center_w

        S = np.sum(w * (self.data_target[self.ci0:self.ci1,1] - \
                   model_interp)**2) / np.sum(w)

        self.data = data_old

        del data_old, data_n, model_interp, w

        return S


    def find(self, **kwargs):

        # Number of points to try for each iteration
        if ('N' in kwargs):
            self.pts = kwargs['N']
        else:
            self.pts = 15

        # Narrowing factor when going to the next iteration
        # pace[0] = narrowing factor for vsini
        # pace[1] = narrowing factor for abundance
        if ('pace' in kwargs):
            self.pace = kwargs['pace']
        else:
            self.pace = np.array([1.5,2.0])

        # Initial guess range for abundance. It has to be a numpy array of
        # length = 2
        if ('a_guess' in kwargs):
            self.a_guess = kwargs['a_guess']
        else:
            self.a_guess = np.array([-0.100,0.100])

        # Initial guess range for vsini. It has to be a numpy array of
        # length = 2
        if ('v_guess' in kwargs):
            self.v_guess = kwargs['v_guess']
        else:
            self.v_guess = np.array([0.5,10.0])

        # Minimum number of iterations
        if ('min_i' in kwargs):
            self.min_i = kwargs['min_i']
        else:
            self.min_i = 3

        # Maximum number of iterations
        if ('max_i' in kwargs):
            self.max_i = kwargs['max_i']
        else:
            self.max_i = 21

        # Convergence limits: a numpy array with length 2, corresponding to the
        # limits of vsini and abundance, respectively
        if ('limits' in kwargs):
            self.limits = kwargs['limits']
        else:
            self.limits = np.array([0.01,0.001])

        # Plot the spectral line fit at the end?
        if ('plot' in kwargs):
            self.plot = kwargs['plot']
        else:
            self.plot = True

        # Lower limit of estimation of vsini
        if ('v_low_limit' in kwargs):
            self.v_low_limit = kwargs['v_low_limit']
        else:
            self.v_low_limit = 0.5

        # Set 'save' to a filename with an extension (e.g. png, eps)
        # Overrides 'plot' to False
        if ('save' in kwargs):
            self.save = kwargs['save']
            self.plot = False
        else:
            self.save = None

        self.silent = False

        self.best_a = np.mean(self.a_guess)
        self.best_v = np.mean(self.v_guess)
        self.it = 1
        self.finish = False
        self.badfit_status = False

        MOOG = Driver(synth_interval = self.spec,\
                      abunds = np.array([[self.Z, self.best_a],]),\
                      obs_wl = self.data[:,0], obs_flux = self.data[:,1],\
                      data = self.data,
                      gauss = self.gauss, macro_v = self.v_m,\
                      star_name = self.name, plot = self.plot,\
                      savefig = self.save,\
                      y_shift_add = self.yadd,\
                      y_shift_mult = self.ymult,\
                      wl_shift = self.xshift,\
                      line_number = self.line_number)

        self.MOOG = MOOG
        self.it2 = [0,0]

        while self.finish == False and self.it < self.max_i:

            self.MOOG.it = self.it

            # Evaluating vsini
            self.v_grid = np.linspace(self.v_guess[0],self.v_guess[1],self.pts)
            self.S = []
            self.S = self.perf_new(self.v_grid, self.best_a, mode = 'vsini')
            self.S_v = self.S

            tck = UnivariateSpline(self.v_grid, self.S, k = 4, s = 0.05)
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
            self.a_grid = np.linspace(self.a_guess[0],self.a_guess[1],self.pts)
            self.a_grid = self.a_grid[np.argsort(self.a_grid)]
            self.S = []
            self.S = self.perf_new(self.best_v, self.a_grid, mode = 'abundance')

            self.S_a = self.S


            tck = UnivariateSpline(self.a_grid, self.S, k = 4, s = 0.1)
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
            if self.silent == False:
                if self.it > self.min_i and self.intern_u[0] == True \
                        and self.intern_u[1] == True and self.go_v == True and \
                        self.go_a == True and self.it2[0]>=2 and self.it2[1]>=2:
                    self.finish = True
                    break
            else:
                if self.it > self.min_i and self.intern_u[0] == True \
                                        and self.intern_u[1] == True:
                    self.finish = True
                    break


            # Setting the new guess. If one of the best values are too near the
            # edges of the previous guess, it will not narrow its new guess range.
            self.v_width = self.v_guess[1]-self.v_guess[0]
            self.a_width = self.a_guess[1]-self.a_guess[0]


            if self.go_v == True:
                self.v_guess = np.array([self.best_v-\
                    self.v_width/2/self.pace[0], self.best_v+\
                        self.v_width/2/self.pace[0]])
            else:
                self.v_guess = np.array([self.best_v-self.v_width/2,\
                    self.best_v+self.v_width/2])

            if self.go_a == True:
                self.a_guess = np.array([self.best_a-\
                    self.a_width/2/self.pace[1], self.best_a+\
                        self.a_width/2/self.pace[1]])
                if np.abs(self.a_guess[1] - self.a_guess[0]) < 0.05:
                    self.a_guess = np.array([self.best_a - 0.025,\
                            self.best_a + 0.025])

            # Checking if the v_guess contains vsini lower than v_low_limit.
            # If True, it will add a value to the array so that the lower limit
            # is equal to the v_low_limit
            if self.v_guess[0] < self.v_low_limit and self.silent == False:
                self.v_guess += self.v_low_limit-self.v_guess[0]

            else:
                self.best_a += np.random.normal(scale=0.001)
                self.a_guess = np.array([self.best_a-self.a_width/2,\
                    self.best_a+self.a_width/2])

            self.it += 1

        # Finalizing the routine

        self.S = self.perf(np.array([self.best_v,self.best_a]))

        # Trigger bad fit warning
        if self.check > self.badfit_tol:
            self.badfit_status = True

        del MOOG

        return self


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


class Driver(object):

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

    model_in: str, optional
        Name of the atmosphere model input file. Default is 'star.mod'.

    linelist_in: str, optional
        Name of the line list input file. Default is 'lines.dat'.

    observed_in: str, optional
        Name of the input file containing the observed spectrum. Default is
        'spectrum.dat'.

    std_out: str, optional
        Name of the standard output file. Default is 'vm_long.out'.

    summary_out: str, optional
        Default is 'vm_li.out'.

    smoothed_out: str, optional
        Name of the output file containing the smoothed synthetic spectrum.
        Default is 'vm_smooth.out'.

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

    def __init__(self, synth_interval, abunds, obs_wl, obs_flux, data, step=0.01, opac=2.0,
                 wl_shift=0.0, v_shift=0.0, y_shift_add=0.0, y_shift_mult=1.0,
                 gauss=0.0, lorentz=0.0, eps=0.6, macro_v=0.0, vsini=0.0,
                 model_in='star.mod', linelist_in='lines.dat',
                 observed_in='spectrum.dat', std_out='vm_long.out',
                 summary_out='vm_li.out', smoothed_out='vm_smooth.out',
                 atmosphere=1, molecules=1, trudamp=1, lines=1, flux=0,
                 damping=0, star_name='Unnamed star', plot=True, savefig=False,
                 line_number = 0):

        self.name = star_name
        self.plot_switch = plot
        self.savefig = savefig
        # Output files
        self.standard_out = './output/%s_long.out' % self.name
        self.summary_out = './output/%s_li.out' % self.name
        self.smoothed_out = './output/%s_smooth.out' % self.name
        self.smoothed_out_new = './output/%s_smooth_new.out' % self.name

        # Input files
        self.model_in = './atm_models/%s_v.atm' % self.name
        self.lines_in = './MOOG_linelist/lines.%s_v.ares' % self.name
        self.observed_in = './Spectra/%s_%d.dat' % (self.name, line_number)

        # Output files
        self.standard_out_moog = '../output/%s_long.out' % self.name
        self.summary_out_moog = '../output/%s_li.out' % self.name
        self.smoothed_out_moog = '../output/%s_smooth.out' % self.name
        self.smoothed_out_new_moog = '../output/%s_smooth_new.out' % self.name

        # Input files
        self.model_in_moog = '../atm_models/%s_v.atm' % self.name
        self.lines_in_moog = '../MOOG_linelist/lines.%s_v.ares' % self.name
        self.observed_in_moog = '../Spectra/%s_%d.dat' % (self.name, line_number)

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

    def create_batch(self):
        """
        Writes the MOOG driver file batch.par
        """
        f = open('./MOOGFEB2017/%s_synth.par' % self.name, 'w')
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
        f.write(" %.1f %.1f %.2f %.1f\n" % (self.syn_start, self.syn_end,
                                            self.step, self.opac))
        f.write("obspectrum  5\n")
        f.write("plotpars  1\n")
        f.write(" %.2f %.2f 0.05 1.05\n" % (self.wl_start, self.wl_end))
        f.write(" %.4f  %.4f  %.3f  %.3f\n" % (self.v_shift, self.wl_shift,
                                               self.y_shift_add,
                                               self.y_shift_mult))
        f.write(" r  %.3f  0.0  %.1f  %.2f  %.1f" % (self.gauss,
                                                     self.dark,
                                                     self.macro_v,
                                                     self.lorentz))
        f.close()
        del f

    def change_vsini(self, grid_v):

        self.create_batch()
        os.system('bash run_moog_synth.bash %s_synth.par' % self.name)
        model = np.loadtxt(self.smoothed_out, skiprows = 2)
        synth_wl = model[:,0]
        synth_flux = model[:,1]

        synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
            self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)


        self.model_vsini = {}

        for vsini in grid_v:
            self.vsini = vsini
            if self.vsini > 0.0:
                conv_flux = pyasl.fastRotBroad(synth_wl, synth_flux, self.dark, self.vsini)

                self.model_vsini[str(vsini)] = np.array([synth_wl, \
                                               conv_flux/max(conv_flux)]).T

                del conv_flux

            else:

                self.model_vsini[str(vsini)] = np.array([syth_wl, np.nan * np.ones(synth_wl.size)]).T

        del model, synth_wl, synth_flux

        return self

    def change_ab(self, a):
        self.abunds = a

        self.create_batch()
        os.system('bash run_moog_synth.bash %s_synth.par' % self.name)
        model = np.loadtxt(self.smoothed_out, skiprows = 2)
        synth_wl = model[:,0]
        synth_flux = model[:,1]

        if self.vsini > 0.0:
            synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

            conv_flux = pyasl.fastRotBroad(synth_wl, synth_flux, self.dark, self.vsini)

            self.model_ab[str(a)] = np.array([synth_wl, conv_flux/max(conv_flux)]).T
            del conv_flux

        else:

            self.model_ab[str(a)] = np.array([syth_wl, np.nan * np.ones(synth_wl.size)]).T

        del synth_wl, synth_flux, model

        return self

    def run(self):

        """
        Used to run MOOG silent.
        """
        self.create_batch()

        os.system('bash run_moog_synth.bash %s_synth.par' % self.name)

        self.model = np.loadtxt(self.smoothed_out, skiprows=2)
        synth_wl = self.model[:,0]
        synth_flux = self.model[:,1]

        if self.vsini > 0.0:
            synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

            conv_flux = pyasl.fastRotBroad(synth_wl, synth_flux, self.dark, self.vsini)

            self.model = np.array([synth_wl, conv_flux/max(conv_flux)]).T

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
        profile[m] = (2. * (1. - self.dark) * (1. - (vz[m] /
                                                     self.vsini) ** 2.)
                      ** 0.5 + 0.5 * np.pi * self.dark *
                      (1. - (vz[m] / self.vsini) ** 2.)) / \
                     (np.pi * self.vsini * (1. - self.dark / 3.))

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
        if ind0 < (n - 1) / 2:
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
        elif ind0 > (n - 1) / 2:
            if (ind0 + 1) % 2 == 0:
                wl = wl[2*(ind0 - (n - 1) / 2) + 1:-1]
                flux = flux[2 * (ind0 - (n - 1) / 2) + 1:-1]
                obs_flux = obs_flux[2 * (ind0 - (n - 1) / 2) + 1:-1]
                obs_wl = obs_wl[2 * (ind0 - (n - 1) / 2) + 1:-1]
            else:
                wl = wl[2 * (ind0 - (n - 1) / 2):]
                flux = flux[2 * (ind0 - (n - 1) / 2):]
                obs_flux = obs_flux[2 * (ind0 - (n - 1) / 2):]
                obs_wl = obs_wl[2 * (ind0 - (n - 1) / 2):]
        return wl, flux, obs_wl, obs_flux



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


class arr_manage(object):

    """
    Used to perform a series of specific management routines to
    numpy-arrays and data arrays and files.
    """

    def find_index(self,target,array):

        """
        This routine finds the index of a value closest to the target in
        a numpy-array.
        """

        self.index = np.searchsorted(array, target, side = 'left')
        return self.index

    def x_set_limits(self,start,end,data2d):

        """
        This routine returns a section of an array given the start and end
        values. These values do not need to be the exact ones found in the
        array.
        """

        self.start_index = np.searchsorted(data2d[:,0], start, side = 'left')
        self.end_index = np.searchsorted(data2d[:,0], end, side = 'left')
        return data2d[self.start_index:self.end_index]



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_broadening(starname, Teff, met, logg, micro, ab_ni, err_T, err_logg, err_met, err_ni):

    vmac, err_vmac = calc_vmac((Teff, err_T), (logg, err_logg))
    vmac += 0.382 + 0.027
    vsini, err_vsini = calc_vsini(starname, Teff, met, logg, micro, vmac, ab_ni, err_met, err_ni)

    #print vmac, err_vmac, vsini, err_vsini

    if len(vsini) == 0 or np.all(vsini == 0.0) or np.all(err_vsini == 0.0):
        i = np.where((vmac != 0.0) & (err_vmac != 0.0))[0]
        vmac_final = np.mean(vmac[i])
        err_vmac_final = np.mean(err_vmac[i])/np.sqrt(float(len(i)))
        vsini_final = 0.0
        err_vsini_final = 0.0

    else:

        i = np.where((vmac != 0.0) & (vsini != 0.0) & (err_vmac != 0.0) & (err_vsini != 0.0))[0]

        vmac_final = np.mean(vmac[i])
        err_vmac_final = np.mean(err_vmac[i])/np.sqrt(float(len(i)))

        vsini_final = np.mean(vsini[i]) + 0.1 + 0.02
        err_vsini_final = np.sqrt(1./(1./(np.sum(err_vsini[i]**2.) + err_vmac_final**2.)))

    del i, vmac, err_vmac, vsini, err_vsini

    return vsini_final, err_vsini_final, vmac_final, err_vmac_final


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_vmac(Teff, logg):
    '''
    Computes vmac
    '''
    vmac_sun = np.array([3.0, 3.2, 3.1, 3.6, 2.9])
    vmac = np.zeros(5)
    err_vmac = np.zeros(5)
    for i in range(5):
        v = vmac_sun[i] - 0.00707*Teff[0] + 9.2422*10**(-7.)*Teff[0]**2. + \
                        10.0 - 1.81*(logg[0] - 4.44) - 0.05
        err_v = np.sqrt(0.1**2. + (9.2422*10**(-7.)*2*Teff[0]-0.00707)**2.*Teff[1]**2. + \
                        1.81**2.*logg[1]**2. + (logg[0] - 4.44)**2.*0.26**2. + 0.03**2.)
        vmac[i] = v
        err_vmac[i] = err_v

        del v, err_v

    return vmac, err_vmac


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_vsini(starname, Teff, met, logg, micro, v_macro, ab_ni, err_met, err_ni):

    def f_polinomial(x, a, b):
        return a*x + b

    def f_gauss(x, a, b, c):
        return a*np.exp(-(x - b)**2./(2.*c**2.))

    def f_total(x, a, b, c, d, e):
        return f_gauss(x, a, b, c) + f_polinomial(x, d, e)

    def cont(t, p_0, p_1, p_2):
        return p_0 + p_1*t + p_2*t*t

    def continuum_det(x, y):
        nx = len(x)
        rejt = 0.98
        ini_cont = [max(y), 0., 0.]
        coefs,_ = curve_fit(cont, x, y, p0 = ini_cont)
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
            ynorm = cont(x, *coefs)

            del vecx, vecy

        ynorm = y/(cont(x, *coefs))

        del coefs

        return ynorm

    # Set lines
    # Creates dic with characteristics of lines
    lines1 = {'name' : 'FeI', 'wave' : 6027.05, 'Z' : 26, 'EP' : 4.076, 'loggf' : -1.09}
    lines2 = {'name' : 'FeI', 'wave' : 6151.62, 'Z' : 26, 'EP' : 2.176, 'loggf' : -3.30}
    lines3 = {'name' : 'FeI', 'wave' : 6165.36, 'Z' : 26, 'EP' : 4.143, 'loggf' : -1.46}
    lines4 = {'name' : 'FeI', 'wave' : 6705.10, 'Z' : 26, 'EP' : 4.607, 'loggf' : -0.98}
    lines5 = {'name' : 'NiI', 'wave' : 6767.77, 'Z' : 28, 'EP' : 1.826, 'loggf' : -2.17}
    lines = (lines1, lines2, lines3, lines4, lines5)

    # Create model atmosphere
    interpol(starname + '_v', Teff, logg, met, micro)

    inst = starname[starname.index('_')+1:]

    lines_ab, dev_ab = calc_ab(starname, Teff, met, logg, micro, ab_ni, err_met, err_ni)

    line_file = './MOOG_linelist/lines.%s_v.ares' % starname

    x, data = pyasl.read1dFitsSpec('./Spectra/%s_res.fits' % starname)

    new_data = np.array([x, data]).T

    # For each line create vsini object
    vsini_lines = np.zeros(5)
    err_vsini_lines = np.zeros(5)
    ab_keys = map(float,lines_ab.keys())
    data_lines = {}
    for l in range(len(lines)):

        data_o = new_data

        w = lines[l]['wave']

        if w in ab_keys:
            info_line = {}
            info_line['vmac'] = v_macro[l]

            kwargs = {'star_name' : starname}

            j = (x > (w - 0.8)) & (x < (w + 0.8))
            #kwargs['perf_radius'] = 35

            x_l = x[j]
            y_l = data[j]

            k = (x > (w - 0.5)) & (x < (w+0.5))
            #k2 = ~k

            i_r = int(np.floor(len(x[k])/2))
            kwargs['perf_radius'] = i_r

            try:

                new_y_l = continuum_det(x_l, y_l)

                popt,_ = curve_fit(f_total, x_l, y_l, p0 = (-max(y_l)/min(y_l), w, 0.2, 0.0, max(y_l)))
                #print w, popt[1]
                if popt[1] > (w + 0.3) or popt[1] < (w - 0.3):
                    x_shift = 0.0
                else:
                    x_shift = w - popt[1]

                new_x_l = x_l + x_shift

            except RuntimeError:

                popt,_ = curve_fit(f_total, x[j], data[j], p0 = (-max(data[k])/min(data[k]), w, 0.2, 0.0, max(data[j])))
                y_polinomial = f_polinomial(x_l, *popt[3:])
                #print w, popt[1]

                if popt[1] > (w + 0.3) or popt[1] < (w - 0.3):
                    x_shift = 0.0
                else:
                    x_shift = w - popt[1]

                new_x_l = x_l + x_shift
                new_y_l1 = y_l/f_total(x_l, *popt)

                popt2,_ = curve_fit(f_polinomial, new_x_l, new_y_l1, p0 = (0.0, 1.0))
                y_polinomial2 = f_polinomial(new_x_l, *popt2)

                new_y_l = y_l/y_polinomial/y_polinomial2

                del popt2, y_polinomial, y_polinomial2, new_y_l1


            kwargs['badfit_tol'] = 30

            '''
            fig,(ax1,ax2) = plt.subplots(2,1)
            ax1.plot(x_l, y_l)
            #ax1.plot(x_l, y_polinomial)
            #ax1.plot(x_l, f_total(x_l, *popt))
            ax1.axvline(w-0.3)
            ax1.axvline(w+0.3)
            ax1.axvline(w)
            ax2.plot(new_x_l, new_y_l)
            #ax2.plot(new_x_l, y_polinomial2)
            fig.savefig('./plots_broadening/%s_%f.pdf' % (starname, w))
            plt.close('all')
            '''

            #press = raw_input()

            ascii.write([new_x_l, new_y_l], './Spectra/%s_%d.dat' % (starname, l),\
                        format = 'fixed_width_no_header', overwrite = True, delimiter = '\t')
            SN = compute_snr(new_x_l, new_y_l, w)

            #del j, x_l, y_l, popt, y_polinomial, x_shift, new_x_l, new_y_l, new_y_l1, popt2, y_polinomial2
            del j, x_l, y_l, popt, new_x_l, new_y_l

            if inst == 'harps':
                gauss = w/115000.# 2.*np.sqrt(2.*np.log(2.))*w/120000.
            elif inst == 'feros' or inst == 'feros_o':
                gauss = w/48000.
            elif inst == 'uves':
                gauss = w/110000.
            elif inst == 'hires':
                gauss = w/67000.
            elif inst == 'coralie':
                gauss = w/60000.
            else:
                gauss = w/60000.

            spec_window = np.array([w-1.0, w+1.0])
            vrot = Vsini(spec_window, gauss, v_macro[l], line_file, l, SN, **kwargs)


            kwargs2 = {'a_guess' : np.array([(lines_ab['%.3f' % w]) - 1.5*max(0.1,dev_ab['%.3f' % w]), (lines_ab['%.3f' % w]) + 1.5*max(0.1,dev_ab['%.3f' % w])]),\
                       'v_guess' : np.array([0.5, 10.]),\
                       'save' : True,\
                       'N' : 30,\
                       'v_low_limit' : 0.1,\
                       'max_i' : 30}


            vrot = vrot.find(**kwargs2)
            if vrot.best_v == kwargs2['v_low_limit']:
                vrot.badfit_status = True

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
            info_line['badfit'] = vrot.badfit_status
            data_lines[str(w)] = info_line
            del info_line
            if vrot.badfit_status == False:
                vsini_lines[l] = vrot.best_v
                err_vsini_lines[l] = np.sqrt(vrot.S)

            del vrot
            del kwargs, spec_window, kwargs2
            os.system('rm -f ./Spectra/%s_%d.dat' % (starname, l))


        new_data = data_o

    plot_vsini(starname, data_lines)
    plot_grids(starname, data_lines)
    plot_paper(starname, data_lines)

    del new_data, data_lines, data_o, lines1, lines2, lines3, lines4, lines5, lines

    os.system('rm -f ./atm_models/' + starname + '_v.atm')
    os.system('rm -f ./output/%s_long.out' % starname)
    os.system('rm -f ./output/%s_li.out' % starname)
    os.system('rm -f ./output/%s_smooth.out' % starname)
    os.system('rm -f ./output/%s_smooth_new.out' % starname)
    os.system('rm -f ./MOOG_linelist/lines.%s_v.ares' % starname)
    os.system('rm -f ./MOOGFEB2017/abfind_%s_v_2.par' % starname)
    os.system('rm -f ./MOOGFEB2017/abfind_%s_v.par' % starname)
    os.system('rm -f ./MOOGFEB2017/%s_synth.par' % starname)
    os.system('rm -f ./output/%s_v.test' % starname)
    os.system('rm -f ./output/%s_v_out.test' % starname)

    return vsini_lines, err_vsini_lines


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def calc_ab(starname, Teff, met, logg, micro, ab_ni, err_met, err_ni):

    moog_linelist(starname)

    cmd = 'cp ./MOOGFEB2017/abfind.par ./MOOGFEB2017/abfind_%s_v.par' % (starname)
    os.system(cmd)

    par = open('./MOOGFEB2017/abfind_%s_v.par' % (starname), 'r')
    par_out = open('./MOOGFEB2017/abfind_%s_v_2.par' % (starname), 'w')
    for linea in par:
        columnas = linea.strip()
        m = re.search(r'standard_out\s*(\S*).*', columnas)
        if m:
            linea = "standard_out  '../output/%s_v.test'\n" % (starname)
        m = re.search(r'summary_out\s*(\S*).*', columnas)
        if m:
            linea = "summary_out   '../output/%s_v_out.test'\n" % (starname)
        m = re.search(r'model_in\s*(\S*).*', columnas)
        if m:
            linea = "model_in      '../atm_models/%s_v.atm'\n" % (starname)
        m = re.search(r'lines_in\s*(\S*).*', columnas)
        if m:
            linea = "lines_in      '../MOOG_linelist/lines.%s_v.ares'\n" % (starname)
        par_out.writelines(linea)
    par.close()
    par_out.close()

    cmd = 'cp ./MOOGFEB2017/abfind_%s_v_2.par ./MOOGFEB2017/abfind_%s_v.par' % (starname, starname)
    os.system(cmd)

    os.system('bash run_moog_abfind.bash abfind_%s_v.par' % starname)

    output = open('./output/' + starname + '_v_out.test')

    ab = {}
    dev = {}

    for linea in output:
        linea = linea.strip()
        m = re.search(r'[a-z]', linea)
        if m == None:
            m = re.search(r'[\d]', linea)
            if m:
                linea = linea.split()
                ID = int(float(linea[1]))
                if ID == 26:
                    ab_id = 7.50 + met
                    dev_id = err_met
                else:
                    ab_id = 6.22 + ab_ni
                    dev_id = err_ni
                ab[linea[0]] = float(linea[6]) - ab_id - met
                dev[linea[0]] = dev_id# float(linea[7])
        del m


    output.close()

    del cmd, par, par_out, output

    return ab, dev


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def moog_linelist(starname):
    linelist = np.genfromtxt('./Spectra/linelist_vsini.dat', dtype = None, skip_header = 2,\
                             names = ('line', 'excit', 'loggf', 'num', 'ion'))
    line = linelist['line']
    excit = linelist['excit']
    loggf = linelist['loggf']
    num = linelist['num']
    ion = linelist['ion']

    file_ew = np.genfromtxt('./EW/%s_vsini.ares' % starname, dtype = None, usecols = (0,4,5),\
                            names = ('line', 'ew', 'ew_err'))

    line_ew = file_ew['line']
    ew = file_ew['ew']

    output = open('./MOOG_linelist/lines.%s_v.ares' % starname, 'w')
    output.writelines(' %s_vsini.ares\n' % starname)
    for i in range(len(line_ew)):
        index = np.where(line == line_ew[i])[0]
        if len(index) > 0:
            index = int(index)
            output.writelines('  %7.2f    %4.1f       %5.2f     %6.3f                       %7.4f\n' %\
                              (line[index], ion[index], excit[index], loggf[index], ew[i]))
        del index

    output.close()
    del linelist, line, excit, loggf, num, ion, file_ew, line_ew, ew, output


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
    if len(k1) > 0:
        sn.append(np.mean(data[k1])/np.std(data[k1]))
    if len(k2) > 0:
        sn.append(np.mean(data[k2])/np.std(data[k2]))

    del k1, k2
    return np.mean(sn)


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def plot_vsini(starname, data_lines):

    def add_plot(ax, data, model, w, vsini, vmac, badfit):
        x_limits = [data[:,0][0], data[:,0][-1]]
        i_m1 = np.where(model[:,0] < data[:,0][0])[0]
        i_m2 = np.where(model[:,0] > data[:,0][-1])[0]
        model[:,0][i_m1] = 1.0
        model[:,0][i_m2] = 1.0
        del i_m1, i_m2
        if badfit:
            r = ax.patch
            r.set_facecolor('red')
            r.set(alpha = 0.2)
            del r
        ax.plot(data[:,0], data[:,1], marker = '.', ls = 'None', color = 'dimgrey')
        ax.plot(model[:,0], model[:,1], color = 'tomato')

        ax.tick_params(axis='both', which='major', labelsize=7)

        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        ax.set_xlabel(r'$\lambda$ ($\AA$)', fontsize = 'medium', fontname = fontname)
        ax.set_ylabel(r'$F_\lambda$ $d\lambda$', fontsize = 'medium', fontname = fontname)

        ax.set_xlim(x_limits)
        ax.locator_params(axis = 'x', nbins=6)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(sx + (ex-sx)/10., sy + (ey-sy)/10., \
                '$\lambda$ = %s $\AA$,\n$v \sin{i}$ = %.2f,\n$v_{macro}$ = %.2f' % (w, vsini, vmac),\
                style='italic',fontsize=7)
        del x_limits, sx, sy, ex, ey
        return ax

    try:

        lines = sorted(data_lines.keys())
        n = len(lines)

        if n == 1:
            cols, rows = 1,1
        elif n == 2:
            cols, rows = 2,1
        elif n == 3:
            cols, rows = 3,1
        elif n == 4:
            cols, rows = 2,2
        elif n == 5:
            cols, rows = 2,3

        fig, ax = plt.subplots(rows, cols, figsize = (cols*4, rows*3))

        if cols == 1 and rows == 1:
            d = data_lines[lines[0]]
            ax = add_plot(ax, d['data'], d['model'], lines[c], d['vsini'], d['vmac'], d['badfit'])
            del d

        elif rows == 1 and cols != 1:
            for c in range(cols):
                ax_l = ax[c]
                d = data_lines[lines[c]]
                ax_l = add_plot(ax_l, d['data'], d['model'], lines[c], d['vsini'], d['vmac'], d['badfit'])
                ax[c] = ax_l
                del d
        else:
            i = 0
            for r in range(rows):
                for c in range(cols):
                    try:
                        ax_l = ax[r][c]
                        d = data_lines[lines[i]]
                        ax_l = add_plot(ax_l, d['data'], d['model'], lines[i], d['vsini'], d['vmac'], d['badfit'])
                        ax[r][c] = ax_l
                        del d
                        i += 1
                    except IndexError:
                        fig.delaxes(ax[r][c])
                        plt.draw()

        plt.tight_layout()
        fig.savefig('./plots_broadening/%s_vsini_new.pdf' % starname)
        plt.close('all')

        del fig, ax, lines, n

    except:
        pass


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def plot_grids(starname, data_lines):

    def plot_grid(ax, value, w, grid, S, yfit, type_g = 'vsini', badfit = False):
        ax.plot(grid, S, '.-', color = 'steelblue')
        ax.plot(grid, yfit, color = 'orange')
        ax.axvline(value, color = 'steelblue')
        ax.locator_params(nbins=4)
        ax.tick_params(axis='both', which = 'major', labelsize = 5)
        if type_g == 'vsini':
            ax.set_xlabel('$v\sin i$ (km/s)', fontsize = 6, fontname = fontname)
            ax.get_yaxis().set_visible(False)
        else:
            ax.set_xlabel('abundance', fontsize = 6, fontname = fontname)
            ax.set_ylabel('$S$', fontsize = 6)

        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(ex - 4.*(ex-sx)/10., ey - (ey-sy)/10., \
                '$\lambda$ = %s $\AA$' % (w),\
                style='italic',fontsize=5, backgroundcolor = 'white')
        del sx, sy, ex, ey
        return ax

    try:

        lines = sorted(data_lines.keys())
        n = len(lines)
        cols, rows = 2, n

        fig, ax = plt.subplots(rows, cols, figsize = (cols*2, rows*2))

        if rows == 1:
            d = data_lines[lines[0]]
            ax[0] = plot_grid(ax[0], d['abundance'], lines[0], d['a_grid'], d['S_a'], d['yfit_a'], type_g = 'abundance', badfit = d['badfit'])
            ax[1] = plot_grid(ax[1], d['vsini'], lines[0], d['v_grid'], d['S_v'], d['yfit_v'], type_g = 'vsini', badfit = d['badfit'])
            del d
        else:
            for r in range(rows):
                d = data_lines[lines[r]]
                ax[r][0] = plot_grid(ax[r][0], d['abundance'], lines[r], d['a_grid'], d['S_a'], d['yfit_a'], type_g = 'abundance', badfit = d['badfit'])
                ax[r][1] = plot_grid(ax[r][1], d['vsini'], lines[r], d['v_grid'], d['S_v'], d['yfit_v'], type_g = 'vsini', badfit = d['badfit'])
                del d

        #plt.tight_layout()
        fig.subplots_adjust(hspace=0.3, wspace=0.1, bottom=0.05, left=0.15, right=0.95, top=0.98)
        fig.savefig('./plots_broadening/%s_a_vsini.pdf' % starname)
        plt.close('all')

        del lines, n, cols, rows, fig, ax

    except:
        pass


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

def plot_paper(starname, data_lines):

    ticks_font = matplotlib.font_manager.FontProperties(family=fontname, style='normal', size=9, weight='medium', stretch='normal')

    def plot_grid(ax, value, w, grid, S, yfit, type_g = 'vsini', badfit = False):
        ax.plot(grid, S, '.-', color = 'dimgrey')
        ax.plot(grid, yfit, color = 'tomato')
        ax.axvline(value, color = 'tomato')
        ax.locator_params(nbins=4)
        ax.tick_params(axis='both', which = 'major', labelsize = 10)
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
        if type_g == 'vsini':
            ax.set_xlabel('$v\sin i$ (km/s)', fontsize = 9, fontname = fontname)
            #ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%e'))
            #ax.get_yaxis().set_visible(False)
        else:
            ax.set_xlabel('abundance', fontsize = 9, fontname = fontname)
        ax.set_ylabel('$S$', fontsize = 9)

        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(ex - 5.*(ex-sx)/11., ey - (ey-sy)/10., \
                '$\lambda$ = %s $\AA$' % (w),\
                style='italic',fontsize=8, backgroundcolor = 'white')
        del sx, sy, ex, ey
        return ax


    def add_plot(ax, data, model, w, vsini, vmac, badfit):
        x_limits = [data[:,0][0], data[:,0][-1]]
        i_m1 = np.where(model[:,0] < data[:,0][0])[0]
        i_m2 = np.where(model[:,0] > data[:,0][-1])[0]
        model[:,0][i_m1] = 1.0
        model[:,0][i_m2] = 1.0
        del i_m1, i_m2
        if badfit:
            r = ax.patch
            r.set_facecolor('red')
            r.set(alpha = 0.2)
            del r
        ax.plot(data[:,0], data[:,1], marker = '.', ls = 'None', color = 'dimgrey')
        ax.plot(model[:,0], model[:,1], color = 'tomato')

        ax.tick_params(axis='both', which='major', labelsize=8)

        [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        for label in ax.get_yticklabels():
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels():
            label.set_fontproperties(ticks_font)
        ax.set_xlabel(r'$\lambda$ ($\AA$)', fontsize = 9, fontname = fontname)
        ax.set_ylabel(r'$F_\lambda$ $d\lambda$', fontsize = 9, fontname = fontname)

        ax.set_xlim(x_limits)
        ax.locator_params(axis = 'x', nbins=6)
        ax.locator_params(axis = 'y', nbins=5)
        sy, ey = ax.get_ylim()
        sx, ex = ax.get_xlim()
        ax.text(sx + (ex-sx)/13., sy + (ey-sy)/10., \
                '$\lambda$ = %s $\AA$\n$v \sin{i}$ = %.2f\n$v_{macro}$ = %.2f' % (w, vsini, vmac),\
                style='italic',fontsize=8)
        del x_limits, sx, sy, ex, ey
        return ax

    try:
        lines = sorted(data_lines.keys())
        n = len(lines)
        cols, rows = 3, n

        fig, ax = plt.subplots(rows, cols, figsize = (cols*4, rows*2.5))

        if rows == 1:
            d = data_lines[lines[0]]
            ax[0] = plot_grid(ax[0], d['abundance'], lines[0], d['a_grid'], d['S_a'], d['yfit_a'], type_g = 'abundance', badfit = d['badfit'])
            ax[1] = plot_grid(ax[1], d['vsini'], lines[0], d['v_grid'], d['S_v'], d['yfit_v'], type_g = 'vsini', badfit = d['badfit'])
            ax[2] = add_plot(ax[2], d['data'], d['model'], lines[c], d['vsini'], d['vmac'], d['badfit'])
            del d
        else:
            for r in range(rows):
                d = data_lines[lines[r]]
                ax[r][0] = plot_grid(ax[r][0], d['abundance'], lines[r], d['a_grid'], d['S_a'], d['yfit_a'], type_g = 'abundance', badfit = d['badfit'])
                ax[r][1] = plot_grid(ax[r][1], d['vsini'], lines[r], d['v_grid'], d['S_v'], d['yfit_v'], type_g = 'vsini', badfit = d['badfit'])
                ax[r][2] = add_plot(ax[r][2], d['data'], d['model'], lines[r], d['vsini'], d['vmac'], d['badfit'])
                del d

        fig.subplots_adjust(hspace=0.3, wspace=0.25, bottom=0.05, left=0.07, right=0.98, top=0.98)
        fig.savefig('./plots_broadening/%s_vsini_paper.pdf' % starname)
        plt.close('all')

        del lines, n, cols, rows, fig, ax

    except:
        pass


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


if __name__ == '__main__':
    import time
    start = time.time()
    vsini, err_vsini, vmac, err_vmac = calc_broadening('sun05_harps', 5754.160000, -0.015000, 4.416480, 0.654240, -0.016000, 39.748644, 0.269645, 0.080712, 0.034000)
    print vsini, err_vsini, vmac, err_vmac
    end = time.time()
    print end-start, (end-start)/60.
