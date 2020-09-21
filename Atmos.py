from __future__ import division
from builtins import range
import os
import math
import re
from dataclasses import dataclass
from past.utils import old_div
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import ascii
from scipy.stats import sigmaclip, linregress
import scipy.odr as ODR
#from interpol_function import interpol
from uncertainties import unumpy
import matplotlib.pyplot as plt
from AtmosInterpol import interpol

@dataclass
class Tolerance:
    ab: float
    ep: float
    dif: float
    rw: float

@dataclass
class AtmosQuantity:
    name: str
    value: float
    hold: bool
    ranges: list
    bounds: tuple
    width: float
    tol: float
    change: float


class atmos:

    def __init__(self, star, hold, init_vals, debug, file_debug, in_errors,
                 set_boundaries, vals_boundaries, tolerance=[0.001, 0.001, 0.001, 0.001],
                 alias='test', one_round=False, read_mode='linearregression'):

        self.starname = star
        self.alias = alias
        self.ini = init_vals
        if len(init_vals) < 4:
            self.ini = [0.0, 5500., 4.36, 1.23]
        self.debug = debug
        self.file_debug = file_debug
        self.change = 'metallicity'
        self.parameters = ['metallicity', 'temperature', 'gravity', 'velocity']
        #self.changepar = [200.0, 0.2, 0.2]
        self.tol = Tolerance(*tolerance)
        self.one_round = one_round
        self.read_mode = read_mode

        self.moog = [0.0, 0.0, 0.0, 0.0]
        self.hold = [i in hold for i in self.parameters]

        boundaries = [(-3.0, 1.0), (3500., 9000.), (0.5, 4.9), (0.0, 5.0)]
        if set_boundaries and vals_boundaries is not None:
            for i_p, p in enumerate(self.parameters):
                if p in vals_boundaries:
                    imin, imax = boundaries[i_p]
                    boundaries[i_p] = (max(vals_boundaries[p][0], imin),
                                       min(vals_boundaries[p][1], imax))

        self.metallicity = AtmosQuantity('metallicity', self.ini[0], self.hold[0],
                                         [-999.0, -999.0, -999.0], boundaries[0],
                                         0.25, self.tol.ab, 0.0)
        self.temperature = AtmosQuantity('temperature', self.ini[1], self.hold[1],
                                         [-999.0, -999.0], boundaries[1], 50., self.tol.ep, 200.0)
        self.gravity = AtmosQuantity('gravity', self.ini[2], self.hold[2],
                                     [-999.0, -999.0], boundaries[2], 0.25, self.tol.dif, 0.2)
        self.velocity = AtmosQuantity('velocity', self.ini[3], self.hold[3],
                                      [-999.0, -999.0], boundaries[3], 0.25, self.tol.rw, 0.2)

        self.nfailed = 0
        self.exception = 1
        self.change_antes = 'metallicity'

        if in_errors:
            self.n_repeat = 100
        else:
            self.n_repeat = 200

        self.params = []
        self.nit = 0
        self.nbreak = 0
        self.nit_total = 0
        self.nout = 0


    @property
    def values(self):
        met = self.metallicity.value
        T = self.temperature.value
        logg = self.gravity.value
        micro = self.velocity.value
        return (met, T, logg, micro)

    @property
    def boundaries(self):
        met = self.metallicity.bounds
        T = self.temperature.bounds
        logg = self.gravity.bounds
        micro = self.velocity.bounds
        return (met, T, logg, micro)

    def write_debug_moog(self):
        if self.debug:
            f = self.file_debug
            f.info('Ran %s with: feh=%.2f, T=%.0f, logg=%.2f, micro=%.2f\n'\
                    '\t\t Obtained: ab=%.3f, ep=%.3f, dif=%.3f, rw=%.3f, nfailed=%d',\
                    self.change, *list(self.values), *self.moog, self.nfailed)
            self.file_debug = f
            del f

    def write_log(self, message):
        if self.debug:
            f = self.file_debug
            f.info(message)
            self.file_debug = f
            del f

    def add_param(self):
        self.params.append(list(self.values))

    def check_correct_vals(self):
        output = self.moog
        if abs(output[0] - list(self.values)[0]) <= self.tol.ab and \
            abs(output[1]) <= self.tol.ep and \
            abs(output[2]) <= self.tol.dif and \
            abs(output[3]) <= self.tol.rw:
            self.write_log('Found right parameters')
            del output
            return -1
        del output
        return 0

    def check_nout(self):
        nout = self.nout
        vals = list(self.values)
        boundaries = list(self.boundaries)

        for v, b in zip(vals[:-1], boundaries[:-1]):
            if v < b[0] or v > b[1]:
                nout += 1
        if nout >= 3:
            self.write_log('[Fe/H], T and log g are out of the possible ranges. '\
                                     'Cannot find final parameters.')
            self.exception = 2
            del nout, vals, boundaries
            return -1
        del nout, vals, boundaries
        return 0

    @property
    def call_params(self):
        return self.metallicity, self.temperature, self.gravity, self.velocity

    def check_nfailed(self):
        if self.nfailed > 0:
            for p in list(self.call_params):
                if p.hold is False:
                    new_p = np.random.normal(p.value, p.width)
                    if new_p > p.bounds[1]:
                        new_p = p.bounds[1]
                    if new_p < p.bounds[0]:
                        new_p = p.bounds[1]
                    p.value = new_p
                    del new_p
                    p.ranges = [-999.0, -999.0]
                    if p.name == 'metallicity':
                        p.ranges = [-999.0, -999.0, p.value]


    def check_nbreak(self):
        if self.nbreak > 5:
            self.exception = 2
            self.write_log('Failed more than 5 times in the models.')
            return -1
        return 0


    def check_params_rep(self):
        params = self.params
        vals = list(self.values)
        if vals in params:
            self.write_log('Parameters have already been used in another iteration.')
            n = self.next_change(self.change_antes)
            self.change = n
            self.nit += 1
        del params, vals

    def check_nrepeat(self):
        if self.nit >= self.n_repeat:
            self.write_log('Parameters were repeated more than %d times.' % (self.n_repeat))
            self.exception = 2
            return -1
        return 0


    def check_nit_total(self):
        if self.nit_total >= 500000:
            self.write_log('More than 500000 iterations for the same star. Stopping.')
            self.exception = 2
            return -1
        return 0


    def check_hold(self, xmetal):
        for i_p, p in enumerate(self.call_params):
            if p.hold:
                self.moog[i_p] = 0.0
                if i_p == 0:
                    self.moog[i_p] = xmetal

    @property
    def show_hold(self):
        return [p.hold for p in self.call_params]


    def new_iteration(self, xmetal):
        self.nit_total += 1
        self.check_hold(xmetal)
        self.nout = 0
        self.add_param()
        self.change_antes = self.change


    def moog_output(self, output, nfail):
        self.moog = output
        self.nfailed = nfail


    def new_values(self, new_vals):
        for i_p, p in enumerate(self.call_params):
            p.value = new_vals[i_p]

    @staticmethod
    def next_change(change_ini):
        c = ['metallicity', 'temperature', 'pressure', 'velocity']
        i = c.index(change_ini)
        if i == 3:
            i = -1
        return c[i+1]

    def new_change(self, change_ini=None):
        if change_ini is None:
            change_ini = self.change
        self.change = self.next_change(change_ini)

    def check_nfailed_it(self, change_ini):
        if self.nfailed > 0:
            self.new_change(change_ini)
            self.write_log('Failed in metallicity. Change=%s' % (self.change))
            self.nbreak += 1

    def change_metallicity(self, new_val):
        self.metallicity.value = new_val

    def check_met(self):
        met = self.metallicity.ranges
        return (met[0] == met[1]) and (met[1] == met[2]) and (met[0] != self.metallicity.value)

    def select_param(self, name):
        for p in (self.metallicity, self.temperature, self.gravity, self.velocity):
            if p.name == name:
                break
        return p

    def change_parameter(self, name_par, moog_output, range_m, decimals):
        m_val = self.moog[moog_output]
        p = self.select_param(name_par)
        ext = p.ranges
        val = p.value

        if m_val > p.tol:
            p.ranges[0] = val
            if val < ext[1]:
                p.value = round(np.mean(ext), decimals)
            else:
                p.value = round(mult(val, range_m, 'upper'), decimals)

        else:
            p.ranges[1] = val
            if ext[0] != -999. and val > ext[0]:
                p.value = round(np.mean(ext), decimals)
            else:
                p.value = round(mult(val, range_m, 'floor'), decimals)

    def runMOOG(self, atmos_values, fesun=7.50):
        m, T, g, vt = atmos_values
        interpol(self.starname, T, g, m, vt, self.alias, fesun=fesun)
        cmd = str("MOOGSILENT > temp.log 2>&1 <<EOF\nMOOGFEB2017/ab_%s.par\n\nEOF" % self.alias)
        os.system(cmd)
        ab, ep, dif, rw, nfailed = compute_average_abundance(self.starname,\
                                                             w=False, alias=self.alias,\
                                                             mode=self.read_mode)
        ab = ab - fesun
        return ab, ep, dif, rw, nfailed

    def check_boundaries(self, atmos_values):
        boundaries = list(self.boundaries)

        for v, b in zip(atmos_values, boundaries):
            if v > b[1] or v < b[0]:
                return False
        return True

    def check_if_hold(self, atmos_values):
        #params = np.copy(np.array(atmos_values))
        for i, p in enumerate(self.call_params):
            if p.hold:
                atmos_values[i] = self.ini[i]
        return atmos_values[0], atmos_values[1:]

    def update_moog(self, moog, nfailed):
        for i_p, p in enumerate(self.call_params):
            if p.hold:
                self.moog[i_p] = 0.0
                if i_p == 0:
                    self.moog[i_p] = self.ini[0]
            else:
                self.moog[i_p] = moog[i_p]
        self.nfailed = nfailed
        return self.moog


    def objective_function_vec(self, X, met):
        boundaries = self.check_boundaries([met, X[0], X[1], X[2]])
        if boundaries:
            ab, ep, dif, rw, nfailed = self.runMOOG([met, X[0], X[1], X[2]])
            ab, ep, dif, rw = self.update_moog([ab, ep, dif, rw], nfailed)
            return ep, rw, ab, ab-dif
        return 10.**20., 10.**20., 10.**20., 10.**20.

    def objective_function(self, X, met):
        boundaries = self.check_boundaries([met, X[0], X[1], X[2]])
        if boundaries:
            ab, ep, dif, rw, nfailed = self.runMOOG([met, X[0], X[1], X[2]])
            ab, ep, dif, rw = self.update_moog([ab, ep, dif, rw], nfailed)
            return 5*((3.5* ep)**2.+(1.3*rw)**2.)+2*(dif)**2.
        return 10.**20.


    def simplex(self, S, met):
        Xm = np.array([0, 0, 0], dtype=float)
        Xr = np.array([0, 0, 0], dtype=float)
        Xe = np.array([0, 0, 0], dtype=float)
        Xc = np.array([0, 0, 0], dtype=float)

        Xm = np.mean(S[:3, 1], axis=0)
        Xr = 2*Xm - S[3][1]
        met, Xr = self.check_if_hold(np.array([met, Xr[0], Xr[1], Xr[2]]))
        fr = self.objective_function(Xr, met)

        if S[0][0] <= fr < S[2][0]:
            S[3][1] = Xr
            S[3][0] = fr

        elif fr < S[0][0]:
            Xe = 3*Xm - 2*S[3][1]
            met, Xe = self.check_if_hold(np.array([met, Xe[0], Xe[1], Xe[2]]))
            fe = self.objective_function(Xe, met)

            if fe < fr:
                S[3][1] = Xe
                S[3][0] = fe
            else:
                S[3][1] = Xr
                S[3][0] = fr

        else:
            Xc = 0.5*(Xm + S[3][1])
            met, Xc = self.check_if_hold(np.array([met, Xc[0], Xc[1], Xc[2]]))
            fc = self.objective_function(Xc, met)

            if fc <= S[3][0]:
                S[3][1] = Xc
                S[3][0] = fc
            else:
                for i in range(1, 4):
                    S[i][1] = 0.5*(S[0][1]+S[i][1])
                    met, S[i][1] = self.check_if_hold(np.array([met, S[i][1][0],\
                                                                S[i][1][1], S[i][1][2]]))
                    S[i][0] = self.objective_function(S[i][1], met)

        S = S[np.argsort(S.T[0])]
        del Xm, Xr, Xe, Xc
        return S

    def nelder_optimizer(self, it_simp, it_res_simp):
        counter = 0
        met, T, logg, vmic = self.ini
        for i in range(it_res_simp):
            log_string = '{:>2s} {:>8s} {:>8s} {:>5s} {:>5s} '\
                         '{:>8s} {:>8s} {:>8s} {:>8s}'.format('It', 'S', 'T', 'logg', 'vt',
                                                              'slp1', 'slp2',
                                                              'af1-af2', 'af1-met')
            self.write_log(log_string)
            self.write_log('{:-^68s}'.format('-'))
            xin = np.array([T, logg, vmic])
            metin = met
            slp1, slp2, af1, af2 = self.objective_function_vec(xin, metin)
            l1 = self.temperature.change if slp1 > 0 else -self.temperature.change
            l2 = self.gravity.change if (af1-af2) > 0 else -self.gravity.change
            l3 = self.velocity.change if slp2 > 0 else -self.velocity.change

            met, X0 = self.check_if_hold(np.array([met, T, logg, vmic]))
            met, X1 = self.check_if_hold(np.array([met, T+l1, logg, vmic]))
            met, X2 = self.check_if_hold(np.array([met, T, logg+l2, vmic]))
            met, X3 = self.check_if_hold(np.array([met, T, logg, vmic+l3]))

            f0 = self.objective_function(X0, met)
            f1 = self.objective_function(X1, met)
            f2 = self.objective_function(X2, met)
            f3 = self.objective_function(X3, met)

            S = np.array([[f0, X0], [f1, X1], [f2, X2], [f3, X3]])

            S = S[np.argsort(S.T[0])]
            if np.any(np.isnan(S.T[0].astype(float))):
                self.write_log('One of the values of S is nan. Stopping the computation.')
                return S[0][1], met, 2
            count_simp = 0
            slp1, slp2, af1, af2 = self.objective_function_vec(S[0][1], met)

            while (np.abs(slp1) > self.tol.ep or\
                  np.abs(slp2) > self.tol.rw or\
                  np.abs(af1 - af2) > self.tol.dif) and\
                  count_simp < it_simp:
                Santes = np.copy(S)
                S = self.simplex(S, met)
                if Santes == S:
                    self.write_log('No change in S from previous value. Stopping cycle')
                    break
                slp1, slp2, af1, af2 = self.objective_function_vec(S[0][1], met)
                self.new_values([af1, S[0][1][0], S[0][1][1], S[0][1][2]])
                count_simp += 1

                for j in range(4):
                    if j == 0:
                        log_string = '{:2d} {: 7.5f} {:7.3f} {:4.3f} {:4.3f} {: 6.5f} {: 6.5f} '\
                                     '{: 6.5f} {: 6.5f}'.format(count_simp, S[0][0], S[0][1][0],
                                                                S[0][1][1], S[0][1][2], slp1,
                                                                slp2, af1-af2, af1-met)
                    else:
                        log_string = '{:2d} {: 7.5f} {:7.3f} {:4.3f} {:4.3f}'.format(count_simp,
                                                                                     S[j][0],
                                                                                     S[j][1][0],
                                                                                     S[j][1][1],
                                                                                     S[j][1][2])
                    self.write_log(log_string)

            if self.check_correct_vals() == -1:
                counter += 1

            T = S[0][1][0]
            logg = S[0][1][1]
            vmic = S[0][1][2]
            met = af1
            self.new_values([met, T, logg, vmic])

            if counter > 1:
                self.write_log('###############################################')
                self.write_log('Coverged at:')
                log_string = '{:>5s} {:>8s} {:>5s} {:>5s} {:>5s} {:>8s} {:>8s} '\
                             '{:>8s} {:>8s}'.format('Cycle', 'T', 'logg', 'vt', 'met', 'slp1',
                                                    'slp2', 'af1-af2', 'af1-met')
                self.write_log(log_string)
                log_string = '{:5d} {:7.3f} {:4.3f} {:4.3f} {:4.3f} {: 6.5f} {: 6.5f} '\
                             '{: 6.5f} {: 6.5f}'.format(i+1, T, logg, vmic, met,
                                                        slp1, slp2, af1-af2, af1-met)
                self.write_log(log_string)

                return S[0][1], met, 1

            if (counter == 1) and self.one_round:
                self.write_log('###############################################')
                self.write_log('Coverged only one time at:')
                log_string = '{:>5s} {:>8s} {:>5s} {:>5s} {:>5s} {:>8s} '\
                             '{:>8s} {:>8s} {:>8s}'.format('Cycle', 'T', 'logg', 'vt', 'met',
                                                           'slp1', 'slp2', 'af1-af2', 'af1-met')
                self.write_log(log_string)
                log_string = '{:5d} {:7.3f} {:4.3f} {:4.3f} {:4.3f} {: 6.5f} {: 6.5f} '\
                             '{: 6.5f} {: 6.5f}'.format(i+1, T, logg, vmic, met, slp1, slp2,
                                                        af1-af2, af1-met)
                self.write_log(log_string)

                return S[0][1], met, 1

            self.write_log('###############################################')
            self.write_log('New Cycle')
            self.write_log('{:>8s} {:>5s} {:>5s} {:>6s}'.format('T', 'logg', 'vt', 'met'))
            self.write_log('{:7.3f} {:4.3f} {:4.3f} {: 4.3f}'.format(T, logg, vmic, af1))

        return S[0][1], met, 2


#******************************************************************************

def runODR(x, y, weights, deg=1):
    isort = np.argsort(x)
    func = ODR.polynomial(deg)
    mydata = ODR.Data(x=x[isort], y=y[isort], we=weights[isort])
    myodr = ODR.ODR(mydata, func)
    myoutput = myodr.run()
    beta = myoutput.beta[1]
    del isort, func, mydata, myodr, myoutput
    return beta

#@profile
def compute_average_abundance(starname, w=False, alias='test', mode='linearregression'):
    nfailed = 0
    flag = -1
    ep = dif = rw = final_Fe = -99

    abund = {'ab':{'FeI':-99, 'FeII': -999},\
             'lines':{'FeI':[], 'FeII':[]}}

    names = ['FeI', 'FeII']
    names_file = ['Fe I', 'Fe II']

    failed1 = re.compile(r'OH NO! ANOTHER FAILED ITERATION!')
    failed2 = re.compile(r'CANNOT DECIDE ON A LINE WAVELENGTH STEP FOR')
    line1 = re.compile(r'[a-z]')
    line2 = re.compile(r'[\d]+   [\d]+.+')
    epline = re.compile(r'E.P. correlation')
    rwline = re.compile(r'R.W. correlation')
    abline = re.compile(r'average abundance')

    with open('./output/%s_out.test' % alias, 'r') as filemoog:
        for line in filemoog:
            line = line.strip()
            if failed1.search(line) or failed2.search(line):
                nfailed += 1
            for p in range(2):
                m = re.search(r'Abundance Results for Species (%s\s)\.*' % names_file[p], line)
                if m:
                    flag = p

            m = line1.search(line)
            if m is None:
                m = line2.search(line)
                if m:
                    abund['lines'][names[flag]].append(line)

            m = epline.search(line)
            if m and flag == 0:
                ep = float(line.split()[4])

            m = rwline.search(line)
            if m and flag == 0:
                rw = float(line.split()[4])

            m = abline.search(line)
            if m:
                abund['ab'][names[flag]] = float(line.split()[3])

    del failed1, failed2, line1, line2, epline, rwline, abline, filemoog

    if mode == 'linearregression':

        for p in names:
            ab = np.array([fe.split()[6] for fe in abund['lines'][p]], dtype=float)
            if p == 'FeI':
                iclip = sigma_clip(ab, maxiters=1)
                a_list = np.array([list(map(fe.split().__getitem__, [2, 5]))\
                                   for fe in abund['lines'][p]], dtype=float).T
                ep_list = a_list[0]
                rw_list = a_list[1]
                isort_ep = np.argsort(ep_list[~iclip.mask])
                isort_rw = np.argsort(rw_list[~iclip.mask])

                ep, _, _, _, _ = linregress(ep_list[~iclip.mask][isort_ep],\
                                            ab[~iclip.mask][isort_ep])
                rw, _, _, _, _ = linregress(rw_list[~iclip.mask][isort_rw],\
                                            ab[~iclip.mask][isort_rw])
                abund['ab'][p] = np.mean(ab[~iclip.mask])
                del iclip, ep_list, rw_list, isort_ep, isort_rw, a_list

            else:
                abund['ab'][p] = np.median(ab)
            del ab

    elif mode == 'odr':

        filename = './EW/%s.txt' % starname
        filelines = ascii.read(filename, include_names=('col1', 'col2', 'col4', 'col5'))
        file_wave = filelines['col1']
        file_ew = filelines['col2']
        file_e_ew = np.maximum(filelines['col4'], filelines['col5'])

        for p in names:
            ab = np.array([fe.split()[6] for fe in abund['lines'][p]], dtype=float)
            if p == 'FeI':

                a_list = np.array([list(map(fe.split().__getitem__, [2, 5]))\
                                   for fe in abund['lines'][p]], dtype=float).T
                ep_list = a_list[0]
                rw_list = a_list[1]

                wave = np.array([fe.split()[0] for fe in abund['lines'][p]], dtype=float)
                ew = np.array([file_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
                ew_err = np.array([file_e_ew[int(np.where(file_wave == wv)[0])] for wv in wave])

                weights = 1./(ew_err/ew)
                weights = weights/np.sum(weights)

                ep = runODR(ep_list, ab, weights=weights)
                rw = runODR(rw_list, ab, weights=weights)

                abund['ab'][p] = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew))).n

                del wave, ew_err, a_list, ep_list, rw_list, weights, ew

            else:
                wave = np.array([fe.split()[0] for fe in abund['lines'][p]], dtype=float)
                ew = np.array([file_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
                ew_err = np.array([file_e_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
                abund['ab'][p] = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew))).n

                del wave, ew_err, ew

            del ab
        del filelines, file_wave, file_ew, file_e_ew

    if w:
        filename = './EW/%s.txt' % starname
        filelines = ascii.read(filename, include_names=('col1', 'col2', 'col3', 'col4', 'col5'))
        file_wave = filelines['col1']
        file_ew = filelines['col2']
        file_e_ew = np.maximum(filelines['col4'], filelines['col5'])

        for p in names:
            a_list = np.array([list(map(fe.split().__getitem__, [0, 6]))\
                               for fe in abund['lines'][p]], dtype=float).T
            wave = a_list[0]
            ab = a_list[1]

            w = np.array([1./file_e_ew[int(np.where(file_wave == wv)[0])] \
                         for wv in wave])

            if sum(w) != 0.:
                abund['ab'][p] = round(np.average(ab, weights=w), 3)
            else:
                abund['ab'][p] = np.mean(ab)

        del filename, filelines, file_wave, file_ew, file_e_ew, wave, ab, w

    dif = abund['ab']['FeI'] - abund['ab']['FeII']
    final_Fe = abund['ab']['FeI']

    del abund, names, names_file

    return final_Fe, ep, dif, rw, nfailed


#******************************************************************************


def mult(x, base, level):
    """
    Finds the multiple of 'base' closer to the number x.
    Input: x: number
           base: base for which we want to compute the closer value
           level: 'upper' or 'floor'
                  the higher or lower multiple of 'base' of x.
    Return: closer multiple of 'base' to the number x.
    """

    num = math.floor(old_div(x, base))
    if level == 'upper':
        final = (num + 1)*base
    else:
        final = num*base
        if final == x:
            final = final - base
    return final


def calc_params(star, hold, init_vals, debug, log_f=None,\
                set_boundaries=False, vals_boundaries=None,\
                in_errors=False, alias='test',
                minimization='per_parameter',
                one_round=False, read_mode='linearregression'):
    """
    Computes the stellar parameters for a certain star.
    Uses the methods: Metallicity, Temperature, Pressure and Velocity.

    Input: star: name of the star.
           hold : array with the values that are not to be changed.
           When empty, all the values will be changed.
           init_vals : initial parameters.
                       When empty, default values will be used.
           debug : True/False,
                   if you want to turn on or off the debugging option.
           file_debug : File where the output from the debugging
                        will be stored.

    Returns: T: temperature of the model.
             logg: surface gravity.
             xmetal: metallicity.
             micro: microturbulence velocity.
             exception: 1 or 2. If 1, then it founds the correct parameters.
                        If 2, it encountered a problem and couldn't converge to a solution.
    """
    # Creates the atmosphere object, where all the data concerning the
    # computation of the atmospheric parameters will be stored.
    a = atmos(star, hold, init_vals, debug, log_f, in_errors, \
              set_boundaries, vals_boundaries, alias=alias,\
              one_round=one_round, read_mode=read_mode)

    h_array = a.show_hold
    a.write_log('hold_m=%s, hold_t=%s, hold_p=%s, hold_v=%s' % \
                          (h_array[0], h_array[1], h_array[2], h_array[3]))
    line_hold_param = ', '.join(['%s = (%.2f, %.2f)' % (n, b[0], b[1])\
                                 for n, b in zip(a.parameters, a.boundaries)])
    a.write_log('Boundaries are: %s' % line_hold_param)
    del line_hold_param

    a.write_log('Initial values are: feh=%.2f, T=%.0f, logg=%.2f, vt=%.2f' % \
                          (init_vals[0], init_vals[1], init_vals[2], init_vals[3]))

    # Modifies the output from MOOG if hold metallicity is 'yes'
    xmetal_i = a.metallicity.value
    a.check_hold(xmetal_i)

    if minimization == 'per_parameter':

        # First iteration with MOOG.
        ab, ep, dif, rw, nfailed = runMOOG(a.starname, a.values, a.alias, a.read_mode)

        a.moog_output([ab, ep, dif, rw], nfailed)
        a.write_log('ab=%.3f, ep=%.3f, dif=%.3f, rw=%.3f, nfailed=%d' % \
                              (ab, ep, dif, rw, nfailed))
        a.metallicity.ranges[2] = ab

        a.write_log('change=%s' % a.change)
        i = 0
        while True:
            # Values that need to be reset each iteration
            a.new_iteration(xmetal_i)

            # Found the right values
            if a.check_correct_vals() == -1:
                break

            # If all the parameters are out of range, break the calculation
            if a.check_nout() == -1:
                break

            # Parameters still in the permitted ranges
            change = a.change

            if change == 'metallicity':
                a = Metallicity(a)
                ab = a.moog[0]
                a.metallicity.ranges = [-999., -999., ab]
                a.temperature.ranges = [-999., -999.]
                a.gravity.ranges = [-999., -999.]
                a.velocity.ranges = [-999., -999.]
                i = 0
                a.check_nfailed_it(change)

            elif change == 'temperature':
                a, i = Temperature(a, i)
                a.check_nfailed_it(change)

            elif change == 'pressure':
                a, i = Pressure(a, i)
                a.check_nfailed_it(change)

            else:
                a, i = Velocity(a, i)
                a.check_nfailed_it(change)

            # If an iteration failed, change the input parameters
            # according to a normal distribution
            a.check_nfailed()

            # If the iteration has failed more than 5 times,
            # break the calculation
            if a.check_nbreak() == -1:
                break

            # If the parameters for an iteration are the same
            # as a previous one, save them
            a.check_params_rep()

            # If the parameters are repeated more than
            # 500 times, break the calculation
            if a.check_nrepeat() == -1:
                break

            # If mnore than 1 million iterations have been performed
            # and the values have not converge, stop the calculation
            if a.check_nit_total() == -1:
                break
            a.write_log('change is %s' % change)


        xmetal, T, logg, micro = list(a.values)
        exception = a.exception
        a.write_log('Final parameters for %s: feh=%.3f, T=%.0f, logg=%.3f, micro=%.3f' %\
                              (star, xmetal, T, logg, micro))

        del a, i, xmetal_i, h_array

    elif minimization == 'downhill_simplex':

        # First iteration with MOOG.
        ab, ep, dif, rw, nfailed = runMOOG(a.starname, a.values, a.alias, a.read_mode)
        ab, ep, dif, rw = a.update_moog([ab, ep, dif, rw], nfailed)

        if a.check_correct_vals() == -1:
            vals = a.values
            del a
            return vals[1], vals[2], vals[0], vals[3], 1

        [T, logg, micro], xmetal, exception = a.nelder_optimizer(60, 10)
        del a

    return T, logg, xmetal, micro, exception


def Metallicity(atm):
    """
    Runs MOOG with a model atmosphere in which the
    metallicity is the value that changes.
    It stops running when the abundances derived by MOOG
    are the same as the input value.
    """

    c = 0
    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atm.write_log('Iteration in Metallicity was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atm.new_change()
            break

        xmetal = atm.metallicity.value
        xmetal_antes = xmetal
        ab, ep, dif, rw = atm.moog
        if abs(ab - xmetal) <= atm.tol.ab:
            atm.new_change()
            break

        else:
            if c > 50:
                if abs(ep) <= atm.tol.ep:
                    atm.new_change('temperature')
                else:
                    atm.new_change()
                break

            atm.change_metallicity(ab)
            xmetal = atm.metallicity.value
            bound_min, bound_max = atm.metallicity.bounds
            if xmetal < bound_min or xmetal > bound_max:
                atm.write_log('Not possible to compute parameters '\
                                          'for [Fe/H] < %.1f or [Fe/H] > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atm.new_change()
                atm.change_metallicity(xmetal_antes)
                break
            ab, ep, dif, rw, nfailed = runMOOG(atm.starname, atm.values,
                                               atm.alias, atm.read_mode)
            atm.moog_output([ab, ep, dif, rw], nfailed)
            atm.write_debug_moog()
            del ab, ep, dif, rw, xmetal, xmetal_antes
            if nfailed > 0:
                atm.new_change()
                break
            c += 1

    del c, nit_total, nfailed
    return atm


#******************************************************************************


def Temperature(atm, i):
    """
    Runs MOOG with a model atmosphere in which the
    temperature is the value that changes.
    It stops running when the correlation between Ab(FeI)
    and excitation potential is less than 0.002
    """

    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atm.write_log('Iteration in Temperature was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atm.new_change()
            break

        Tantes = atm.temperature.value
        ab, ep, dif, rw = atm.moog

        if abs(ep) <= atm.tol.ep:
            if ab == atm.metallicity.value:
                atm.new_change()
            else:
                atm.new_change('velocity')
            break

        else:
            atm.change_parameter('temperature', 1, 250., 0)
            T = atm.temperature.value
            atm.write_log('T is %.0f, Tantes is %.0f' % (T, Tantes))
            bound_min, bound_max = atm.temperature.bounds
            if T > bound_max or T < bound_min:
                atm.write_log('Not possible to compute parameters '\
                                          'for T < %d or T > %d. '\
                                          'Check the boundaries of your parameter.' % \
                                          (int(bound_min), int(bound_max)))
                atm.new_change()
                atm.temperature.value = Tantes #???????????????????
                if T < 3500.:
                    atm.temperature.value = 3500.
                break

            ab, ep, dif, rw, nfailed = runMOOG(atm.starname, atm.values,
                                               atm.alias, atm.read_mode)
            atm.moog_output([ab, ep, dif, rw], nfailed)

            atm.write_debug_moog()
            atm.metallicity.ranges[i] = ab
            del ab, ep, dif, rw, T, Tantes
            i += 1
            if i == 3:
                i = 0

            if nfailed > 0:
                atm.new_change()
                break

            if atm.check_met():
                atm.new_change('velocity')
                break

    del nit_total, nfailed
    return atm, i


#******************************************************************************


def Pressure(atm, i):
    """
    Runs MOOG with a model atmosphere in which the
    surface gravity is the value that changes.
    It stops running when the difference between the abundances
    derived for FeI and FeII is less than 0.002.
    """

    nfailed = 0
    nit_total = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atm.write_log('Iteration in Pressure was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atm.new_change()
            break

        logg_antes = atm.gravity.value
        ab, ep, dif, rw = atm.moog

        if abs(dif) <= atm.tol.dif:
            if ab == atm.metallicity.value:
                atm.new_change()
            else:
                atm.new_change('velocity')
            break

        else:
            atm.change_parameter('gravity', 2, 0.25, 5)
            logg = atm.gravity.value
            bound_min, bound_max = atm.gravity.bounds
            if logg < bound_min or logg > bound_max:
                atm.write_log('Not possible to compute parameters '\
                                          'for log g < %.1f or log g > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atm.new_change()
                atm.gravity.value = logg_antes
                break

            ab, ep, dif, rw, nfailed = runMOOG(atm.starname, atm.values,
                                               atm.alias, atm.read_mode)
            atm.moog_output([ab, ep, dif, rw], nfailed)
            atm.write_debug_moog()
            atm.metallicity.ranges[i] = ab
            del ab, ep, dif, rw, logg, logg_antes
            i += 1
            if i == 3:
                i = 0

            if nfailed > 0:
                atm.new_change()
                break

            if atm.check_met():
                atm.new_change('velocity')
                break

    del nit_total, nfailed

    return atm, i


#******************************************************************************


def Velocity(atm, i):
    """
    Runs MOOG with a model atmosphere in which the
    microturbulence velocity is the value that changes.
    It stops running when the correlation between Ab(FeI)
    and reduced EW is less than 0.002
    """

    nfailed = 0
    nit_total = 0
    v = 0

    while True:
        nit_total += 1
        if nit_total > 100000:
            atm.write_log('Iteration in Velocity was completed '\
                                      'more than 100000 times. '\
                                      'Stopping the computation.')
            atm.new_change()


        micro_antes = atm.velocity.value
        ab, ep, dif, rw = atm.moog
        v += 1

        if v == 50:
            atm.new_change()
            break

        if abs(rw) <= atm.tol.rw:
            if ab == atm.metallicity.value:
                atm.new_change('metallicity')
            else:
                atm.new_change('velocity')
            break

        else:
            atm.change_parameter('velocity', 3, 0.25, 5)
            micro = atm.velocity.value
            bound_min, bound_max = atm.velocity.bounds
            if micro < bound_min or micro > bound_max:
                atm.write_log('Not possible to compute parameters '\
                                          'for micro < %.1f or micro > %.1f. '\
                                          'Check the boundaries of your parameter.' % \
                                          (bound_min, bound_max))
                atm.new_change()
                atm.velocity.value = micro_antes
                break

            ab, ep, dif, rw, nfailed = runMOOG(atm.starname, atm.values,
                                               atm.alias, atm.read_mode)
            atm.moog_output([ab, ep, dif, rw], nfailed)
            atm.write_debug_moog()
            atm.metallicity.ranges[i] = ab
            del ab, ep, dif, rw, micro_antes, micro
            i += 1
            if i == 3:
                i = 0
            if nfailed > 0:
                atm.new_change()
                break
            if atm.check_met():
                atm.new_change('velocity')
                break

    del nit_total, nfailed, v

    return atm, i


#******************************************************************************


def runMOOG(starname, atmos_values, alias='test', read_mode='linearregression', fesun=7.50):
    m, T, g, vt = atmos_values
    interpol(starname, T, g, m, vt, alias, fesun=fesun)
    cmd = str("MOOGSILENT > temp.log 2>&1 <<EOF\nMOOGFEB2017/ab_%s.par\n\nEOF" % alias)
    os.system(cmd)
    ab, ep, dif, rw, nfailed = compute_average_abundance(starname, w=False,
                                                         alias=alias, mode=read_mode)
    ab = ab - fesun
    return ab, ep, dif, rw, nfailed


elements = np.array(['H','He','Li','Be','B','C','N','O','F','Ne',
                     'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca',
                     'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                     'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',
                     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
                     'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
                     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
                     'Lu','Hf','Ta','Wl','Re','Os','Ir','Pt','Au','Hg',
                     'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
                     'Pa','U','Np','Pu','Am'])
            
# From Asplund et al. (2009, Ann. Rev. Ast. Ap., 47, 481)
solarabund = np.array([12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
                        6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
                        3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
                        3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
                        1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
                        1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
                       -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
                        0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
                        0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
                       -5.00,-0.54,-5.00,-5.00,-5.00])

def runMOOG_ab(starname, T, g, m, vt, use_w=False, alias='test', nions=14,\
               nameions=['NaI', 'MgI', 'AlI', 'SiI', 'CaI', 'TiI', 'TiII', 'CrI', 'MnI',\
                         'FeI', 'FeII', 'NiI', 'CuI', 'ZnI']):
    use_w = True
    interpol(starname, T, g, m, vt, alias)
    cmd = 'MOOGSILENT > temp.log 2>&1 << EOF\nMOOGFEB2017/ab_%s.par%sEOF'\
           % (alias, '\n'*nions)
    os.system(cmd)

    abund = {}
    lines = {}
    wave = {}
    ab = {}
    dev = {}
    
    names_file = [n[:n.index('I')]+' '+n[n.index('I'):] for n in nameions]
    element_name = [n[:n.index('I')] for n in nameions]
    
    for n in nameions:
        lines[n] = np.array([])
        wave[n] = np.array([])
        ab[n] = np.nan
        dev[n] = np.nan
        
    abund['lines'] = lines
    abund['wave'] = wave
    abund['ab'] = ab
    abund['dev'] = dev
    
    del lines, wave, ab, dev

    flagab = 0
    with open('./output/%s_out.test' % alias) as output:
        for linea in output:
            linea = linea.strip()
            for p in range(len(nameions)):
                m = re.search(r'Abundance Results for Species ' + names_file[p] + '\.*', linea)
                if m:
                    flagab = p+1

            m = re.search(r'[a-z]', linea)
            if m is None:
                m = re.search(r'[\d]+   [\d]+.+', linea)
                if m:
                    linea = linea.split()
                    wave_linea = float(linea[0])
                    ab_linea = float(linea[6])
                    dev_linea = float(linea[7])
                    if abs(dev_linea) < 10.0:
                        abund['lines'][nameions[flagab-1]] = np.append(abund['lines'][nameions[flagab-1]],\
                                                                       ab_linea)
                        abund['wave'][nameions[flagab-1]] = np.append(abund['wave'][nameions[flagab-1]],\
                                                                      wave_linea)

    del names_file

    for p in nameions:
        c, low, upp = sigmaclip(abund['lines'][p], 1.5, 1.5)
        index = np.where(np.in1d(abund['lines'][p], c))[0]
        if index.size > 0:
            abund['lines'][p] = abund['lines'][p][index]
            abund['wave'][p] = abund['wave'][p][index]
        del c, low, upp, index

    if use_w:
        filename = './EW/' + starname + '.txt'
        filelines = ascii.read(filename, include_names=('col1', 'col2', 'col4', 'col5'))
        file_wave = np.array([np.rint(v*100.0)/100.0 for v in filelines['col1']])
        file_ew = filelines['col2']
        file_e_ew = np.maximum(filelines['col4'], filelines['col5'])
        #print(file_wave)

        w = {}
        for p in nameions:
            #print(p)
            #print(np.std(abund['lines'][p]))
            #w[p] = np.array([1./file_e_ew[np.argmin(np.abs(abund['wave'][p][i]-file_wave))]\
            #                 for i in range(len(abund['lines'][p]))])

            #if sum(w[p]) != 0:
            #    abund['ab'][p] = round(np.average(abund['lines'][p], weights=w[p]), 3)
            #    abund['dev'][p] = round(np.sqrt(np.average((abund['lines'][p] - \
            #                                                np.mean(abund['lines'][p]))**2., \
            #                                                weights=w[p])), 3)
            #else:
            if len(abund['lines'][p]) > 0:
                m = unumpy.uarray(abund['lines'][p], 0.2*np.ones(abund['lines'][p].size))
                #abund['ab'][p] = np.mean(abund['lines'][p])
                #abund['dev'][p] = np.std(abund['lines'][p])
                #print(np.mean(m))
                abund['ab'][p] = np.mean(m).n
                abund['dev'][p] = np.mean(m).s
            else:
                abund['ab'][p] = np.nan
                abund['dev'][p] = np.nan

        del w, filelines, file_wave, file_ew, file_e_ew

    else:
        for p in nameions:
            if abund['lines'][p]:
                abund['ab'][p] = np.mean(abund['lines'][p])
                abund['dev'][p] = np.std(abund['lines'][p])
                
    abund_dict = {}
    for ii, n in enumerate(nameions):
        solar = 0.0
        if element_name[ii] in elements:
            solar = solarabund[np.where(elements == element_name[ii])[0]]
        if np.isnan(abund['ab'][n]):
            abund_dict[n] = [abund['ab'][n], abund['dev'][n], len(abund['lines'][n])]
        else:
            abund_dict[n] = [abund['ab'][n]-solar, abund['dev'][n], len(abund['lines'][n])]

    del abund
    #print(abund_dict)
    return abund_dict


def plot_output_file(starname):
    # Read the result from the moog output file
    nfailed = 0
    flag = -1
    ep_moog = rw_moog = -99

    abund = {'ab':{'FeI':-99, 'FeII': -999},\
             'lines':{'FeI':[], 'FeII':[]}}

    names = ['FeI', 'FeII']
    names_file = ['Fe I', 'Fe II']

    failed1 = re.compile(r'OH NO! ANOTHER FAILED ITERATION!')
    failed2 = re.compile(r'CANNOT DECIDE ON A LINE WAVELENGTH STEP FOR')
    line1 = re.compile(r'[a-z]')
    line2 = re.compile(r'[\d]+   [\d]+.+')
    epline = re.compile(r'E.P. correlation')
    rwline = re.compile(r'R.W. correlation')
    abline = re.compile(r'average abundance')

    with open('output/MOOG_output_files/%s_out.test' % starname, 'r') as filemoog:
        for line in filemoog:
            line = line.strip()
            if failed1.search(line) or failed2.search(line):
                nfailed += 1
            for p in range(2):
                m = re.search(r'Abundance Results for Species (%s\s)\.*' % names_file[p], line)
                if m:
                    flag = p

            m = line1.search(line)
            if m is None:
                m = line2.search(line)
                if m:
                    abund['lines'][names[flag]].append(line)

            m = epline.search(line)
            if m and flag == 0:
                ep_moog = float(line.split()[4])

            m = rwline.search(line)
            if m and flag == 0:
                rw_moog = float(line.split()[4])

            m = abline.search(line)
            if m:
                abund['ab'][names[flag]] = float(line.split()[3])

    del failed1, failed2, line1, line2, epline, rwline, abline, filemoog

    ab_FeI_moog = abund['ab']['FeI']
    ab_FeII_moog = abund['ab']['FeII']

    # Read the EW information
    filelines = ascii.read('EW/%s.txt' % starname, include_names=('col1', 'col2', 'col4', 'col5'))
    file_wave = filelines['col1']
    file_ew = filelines['col2']
    file_e_ew = np.maximum(filelines['col4'], filelines['col5'])

    fig, ax = plt.subplots(3, 1, figsize=(10, 7))

    for p in names:
        ab = np.array([fe.split()[6] for fe in abund['lines'][p]], dtype=float)
        if p == 'FeI':
            # Using sigma clipping and linear regression
            iclip = sigma_clip(ab, maxiters=1)
            a_list = np.array([list(map(fe.split().__getitem__, [2, 5]))\
                               for fe in abund['lines'][p]], dtype=float).T
            ep_list = a_list[0]
            rw_list = a_list[1]
            isort_ep = np.argsort(ep_list[~iclip.mask])
            isort_rw = np.argsort(rw_list[~iclip.mask])
            wave = np.array([fe.split()[0] for fe in abund['lines'][p]], dtype=float)
            ew = np.array([file_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
            ew_err = np.array([file_e_ew[int(np.where(file_wave == wv)[0])] for wv in wave])

            ep, ep_intercept, _, _, ep_err = linregress(ep_list[~iclip.mask][isort_ep],
                                                        ab[~iclip.mask][isort_ep])
            rw, rw_intercept, _, _, rw_err = linregress(rw_list[~iclip.mask][isort_rw],
                                                        ab[~iclip.mask][isort_rw])
            abund['ab'][p] = np.mean(ab[~iclip.mask])

            # Using ODR
            isort_ep_odr = np.argsort(ep_list)
            isort_rw_odr = np.argsort(rw_list)
            weights = 1./(ew_err/ew)
            weights = weights/np.sum(weights)
            func = ODR.polynomial(1)
            mydata = ODR.Data(x=ep_list[isort_ep_odr], y=ab[isort_ep_odr],
                              we=weights[isort_ep_odr])
            myodr = ODR.ODR(mydata, func)
            myoutput = myodr.run()
            ep_odr = myoutput.beta[1]
            ep_intercept_odr = myoutput.beta[0]

            mydata = ODR.Data(x=rw_list[isort_rw_odr], y=ab[isort_rw_odr],
                              we=weights[isort_rw_odr])
            myodr = ODR.ODR(mydata, func)
            myoutput = myodr.run()
            rw_odr = myoutput.beta[1]
            rw_intercept_odr = myoutput.beta[0]

            del weights, func, mydata, myodr, myoutput

            ab_FeI_w = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew)))

            ax[0].plot(ep_list, ab, ls='None', marker='o', color='steelblue')
            ax[0].plot(ep_list[~iclip.mask], ab[~iclip.mask], ls='None', marker='o', color='orange')
            ax[0].plot(ep_list[~iclip.mask][isort_ep],
                       ep*ep_list[~iclip.mask][isort_ep]+ep_intercept,
                       color='red', label='ep slope = %.4f $\pm$ %.4f' % (ep, ep_err))
            ax[0].plot(ep_list[~iclip.mask][isort_ep],
                       ep_moog*ep_list[~iclip.mask][isort_ep]+ep_intercept,
                       ls='--', color='red', label='ep slope moog = %.4f' % ep_moog)
            ax[0].plot(ep_list[~iclip.mask][isort_ep],
                       ep_odr*ep_list[~iclip.mask][isort_ep]+ep_intercept_odr,
                       ls='--', color='green', label='ep slope ODR = %.4f' % ep_odr)

            ax[1].plot(rw_list, ab, ls='None', marker='o', color='steelblue')
            ax[1].plot(rw_list[~iclip.mask], ab[~iclip.mask], ls='None', marker='o', color='orange')
            ax[1].plot(rw_list[~iclip.mask][isort_rw],
                       rw*rw_list[~iclip.mask][isort_rw]+rw_intercept,
                       color='red', label='rw slope = %.4f $\pm$ %.4f' % (rw, rw_err))
            ax[1].plot(rw_list[~iclip.mask][isort_rw],
                       rw_moog*rw_list[~iclip.mask][isort_rw]+rw_intercept,
                       ls='--', color='red', label='rw slope moog = %.4f' % rw_moog)
            ax[1].plot(rw_list[~iclip.mask][isort_rw],
                       rw_odr*rw_list[~iclip.mask][isort_rw]+rw_intercept_odr,
                       ls='--', color='green', label='rw slope odr = %.4f' % rw_odr)

            ax[2].plot(wave[iclip.mask], ab[iclip.mask], ls='None', marker='x', color='steelblue')
            ax[2].plot(wave[~iclip.mask], ab[~iclip.mask], ls='None', marker='o', color='steelblue')
            ax[2].axhline(np.mean(ab[~iclip.mask]), color='steelblue',
                          label='FeI = %.4f' % np.mean(ab[~iclip.mask]))
            ax[2].axhline(ab_FeI_w.n, color='steelblue', ls=':',
                          label='FeI ODR = %.4f' % ab_FeI_w.n)
            ax[2].axhline(ab_FeI_moog, ls='--', color='steelblue',
                          label='FeI moog = %.4f' % ab_FeI_moog)

            del iclip, ep_list, rw_list, isort_ep, isort_rw, a_list, isort_ep_odr, isort_rw_odr,\
                wave, ew_err, ew

        else:
            abund['ab'][p] = np.median(ab)
            wave = np.array([fe.split()[0] for fe in abund['lines'][p]], dtype=float)
            ew = np.array([file_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
            ew_err = np.array([file_e_ew[int(np.where(file_wave == wv)[0])] for wv in wave])
            ab_FeII_w = np.mean(unumpy.uarray(ab, (ew_err/ew)/np.sum(ew_err/ew)))

            ax[2].plot(wave, ab, ls='None', marker='o', color='orange')
            ax[2].axhline(np.mean(ab), color='orange',
                          label='FeII = %.4f\ndif = %.4f' %\
                                (np.median(ab), abund['ab']['FeI'] - abund['ab']['FeII']))
            ax[2].axhline(ab_FeII_w.n, color='orange', ls=':',
                          label='FeII ODR = %.4f\ndif ODR = %.4f' % (ab_FeII_w.n,
                                                                     ab_FeI_w.n - ab_FeII_w.n))
            ax[2].axhline(ab_FeII_moog, color='orange', ls='--',
                          label='FeII moog = %.4f\ndif moog = %.4f' % (ab_FeII_moog,
                                                                       ab_FeI_moog - ab_FeII_moog))

            del wave, ew_err, ew

    ax[0].set_xlabel('Excitation Potential')
    ax[0].set_ylabel('FeI abundance')
    ax[1].set_xlabel('Reduced Equivalent Width')
    ax[1].set_ylabel('FeI abundance')
    ax[2].set_xlabel('Wavelength')
    ax[2].set_ylabel('Abundance')

    ax[0].legend(loc='upper left', fontsize='x-small')
    ax[1].legend(loc='upper left', fontsize='x-small')
    ax[2].legend(loc='upper left', ncol=3, fontsize='x-small')

    fig.subplots_adjust(hspace=0.35, left=0.08, right=0.95, top=0.98, bottom=0.1)
    fig.savefig('output/MOOG_output_files/%s.pdf' % starname)
    plt.close('all')

    del filelines, file_wave, file_e_ew, abund, fig, ab
