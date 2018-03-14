import os, math
import numpy as np

#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


def create_atmos_params(star, hold, init_vals, debug, log_f,\
                        set_boundaries, vals_boundaries):
    """
    Creates de dictionary where all the values relating the
    atmospheric parameters will be store.
    """
    atmos_params = {}
    atmos_params['star'] = star
    atmos_params['ini'] = init_vals
    atmos_params['debug'] = debug
    atmos_params['file_debug'] = log_f
    atmos_params['change'] = 'metallicity'
    moog_output = [0., 0., 0., 0.]

    dic_boundaries = ['']
    if set_boundaries == True and vals_boundaries != None:
        dic_boundaries = vals_boundaries.keys()

    if len(init_vals) < 4:
        init_vals = [5500., 4.36, 0.0, 1.23]

    temperature = {}
    temperature['value'] = init_vals[0]
    if 'temperature' in hold:
        temperature['hold'] = 'yes'
        moog_output[1] = 0.0
    else:
        temperature['hold'] = 'no'
    temperature['ranges'] = [-999., -999.]
    temperature['name'] = 'temperature'
    if 'temperature' in dic_boundaries:
        bound_min = vals_boundaries['temperature'][0]
        bound_max = min(vals_boundaries['temperature'][1], 9500.)
        temperature['boundaries'] = (bound_min, bound_max)
    else:
        temperature['boundaries'] = (3500., 9000.)


    gravity = {}
    gravity['value'] = init_vals[1]
    if 'pressure' in hold:
        gravity['hold'] = 'yes'
        moog_output[2] = 0.0
    else:
        gravity['hold'] = 'no'
    gravity['ranges'] = [-999., -999.]
    gravity['name'] = 'pressure'
    if 'gravity' in dic_boundaries:
        bound_min = max(vals_boundaries['gravity'][0], 1.0)
        bound_max = min(vals_boundaries['gravity'][1], 5.0)
        gravity['boundaries'] = (bound_min, bound_max)
    else:
        gravity['boundaries'] = (0.5, 4.9)



    metallicity = {}
    metallicity['value'] = init_vals[2]
    if 'metallicity' in hold:
        metallicity['hold'] = 'yes'
        moog_output[0] = metallicity['value']
    else:
        metallicity['hold'] = 'no'
    metallicity['ranges'] = [-999., -999., -999.]
    metallicity['name'] = 'metallicity'
    if 'metallicity' in dic_boundaries:
        bound_min = max(vals_boundaries['metallicity'][0], -3.0)
        bound_max = min(vals_boundaries['metallicity'][1], 1.0)
        metallicity['boundaries'] = (bound_min, bound_max)
    else:
        metallicity['boundaries'] = (-3.0, 1.0)



    velocity = {}
    velocity['value'] = init_vals[3]
    if 'velocity' in hold:
        velocity['hold'] = 'yes'
        moog_output[3] = 0.0
    else:
        velocity['hold'] = 'no'
    velocity['ranges'] = [-999., -999.]
    velocity['name'] = 'velocity'
    if 'velocity' in dic_boundaries:
        bound_min = max(vals_boundaries['velocity'][0], 0.0)
        bound_max = min(vals_boundaries['velocity'][1], 5.0)
        velocity['boundaries'] = (bound_min, bound_max)
    else:
        velocity['boundaries'] = (0.0, 5.0)



    atmos_params['temperature'] = temperature
    atmos_params['gravity'] = gravity
    atmos_params['metallicity'] = metallicity
    atmos_params['velocity'] = velocity


    atmos_params['moog'] = moog_output
    atmos_params['nfailed'] = 0
    atmos_params['exception'] = 1
    atmos_params['change_antes'] = atmos_params['change']

    return atmos_params


#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************


class atmos:

    def __init__(self, star, hold, init_vals, debug, file_debug, in_errors, set_boundaries, vals_boundaries):
        self.atmos_params = create_atmos_params(star, hold, init_vals, debug, file_debug, set_boundaries, vals_boundaries)
        if in_errors:
            self.n_repeat = 100
        else:
            self.n_repeat = 200

        self.params = []
        self.nit = 0
        self.nbreak = 0
        self.nit_total = 0
        self.nout = 0

    def data(self):
        return self.atmos_params

    def values(self):
        T = self.atmos_params['temperature']['value']
        logg = self.atmos_params['gravity']['value']
        met = self.atmos_params['metallicity']['value']
        micro = self.atmos_params['velocity']['value']
        return (T, logg, met, micro)

    def boundaries(self):
        T_r = self.atmos_params['temperature']['boundaries']
        logg_r = self.atmos_params['gravity']['boundaries']
        met_r = self.atmos_params['metallicity']['boundaries']
        micro_r = self.atmos_params['velocity']['boundaries']
        return (T_r, logg_r, met_r, micro_r)

    def write_debug_moog(self):
        if self.atmos_params['debug']:
            f = self.atmos_params['file_debug']
            name = self.atmos_params['change']
            vals = self.values()
            output = self.atmos_params['moog']
            nfailed = self.atmos_params['nfailed']
            f.debug('Ran %s with: T=%.2f, logg=%.2f, xmetal=%.2f, micro=%.2f\n'\
                    '\t\t Obtained: ab=%.3f, ep=%.3f, dif=%.3f, rw=%.3f, nfailed=%d.',\
                    name, vals[0], vals[1], vals[2], vals[3], \
                    output[0], output[1], output[2], output[3], nfailed)
            self.atmos_params['file_debug'] = f

            del f, name, vals, output, nfailed

    def write_debug_message(self, message):
        if self.atmos_params['debug']:
            f = self.atmos_params['file_debug']
            f.debug(message)
            self.atmos_params['file_debug'] = f
            del f

    def add_param(self):
        self.params.append(self.values())

    def check_correct_vals(self):
        output = self.atmos_params['moog']
        if abs(output[0] - self.values()[2]) <= 0.02 and \
            abs(output[1]) <= 0.002 and \
            abs(output[2]) <= 0.002 and \
            abs(output[3]) <= 0.002:

            self.write_debug_message('Found right parameters')

            del output

            return -1

        del output

        return 0



    def check_nout(self):
        nout = self.nout
        vals = self.values()
        boundaries = self.boundaries()

        if vals[2] < boundaries[2][0] or vals[2] > boundaries[2][1]: nout += 1
        if vals[1] < boundaries[1][0] or vals[1] > boundaries[1][1]: nout += 1
        if vals[0] > boundaries[0][1] or vals[0] < boundaries[0][0]: nout += 1
        if nout == 3: # All parameters all out of ranges.
            self.write_debug_message('[Fe/H], T and log g are out of the possible ranges. Can not find final parameters.')

            self.atmos_params['exception'] = 2
            del vals, nout

            return -1

        del nout, vals
        return 0


    def check_nfailed(self):
        if self.atmos_params['nfailed'] > 0:
            T = self.atmos_params['temperature']
            if T['hold'] == 'no':
                new_T = np.random.normal(T['value'], 50.)
                if new_T > T['boundaries'][1]: new_T = T['boundaries'][1]
                if new_T < T['boundaries'][0]: new_T = T['boundaries'][0]
                self.atmos_params['temperature']['value'] = new_T
            del T, new_T

            logg = self.atmos_params['gravity']
            if logg['hold'] == 'no':
                new_logg = np.random.normal(logg['value'], 0.25)
                if new_logg > logg['boundaries'][1]: new_logg = logg['boundaries'][1]
                if new_logg < logg['boundaries'][0]: new_logg = logg['boundaries'][0]
                self.atmos_params['gravity']['value'] = new_logg
            del T, new_T

            xmetal = self.atmos_params['metallicity']
            if xmetal['hold'] == 'no':
                new_xmetal = np.random.normal(xmetal['value'], 0.25)
                if new_xmetal < xmetal['boundaries'][0]: new_xmetal = xmetal['boundaries'][0]
                if new_xmetal > xmetal['boundaries'][1]: new_xmetal = xmetal['boundaries'][1]
                self.atmos_params['metallicity']['value'] = new_xmetal
            del xmetal, new_xmetal

            micro = self.atmos_params['velocity']
            if micro['hold'] == 'no':
                new_micro = np.random.normal(micro['value'], 0.25)
                if new_micro < micro['boundaries'][0]: new_micro = micro['boundaries'][0]
                if new_micro > micro['boundaries'][1]: new_micro = micro['boundaries'][1]
                self.atmos_params['velocity']['value'] = new_micro
            del micro, new_micro

            self.atmos_params['metallicity']['ranges'] = [-999., -999., self.atmos_params['metallicity']['value']]
            self.atmos_params['temperature']['ranges'] = [-999., -999.]
            self.atmos_params['gravity']['ranges'] = [-999., -999.]
            self.atmos_params['velocity']['ranges'] = [-999., -999.]


    def check_nbreak(self):
        if self.nbreak > 5:
            self.atmos_params['exception'] = 2
            self.write_debug_message('Failed more than 5 times in the models.')

            return -1

        return 0


    def check_params_rep(self):
        params = self.params
        vals = self.values()
        for q in range(len(params)):
            if params[q][0] == vals[0] and params[q][1] == vals[1] and params[q][2] == vals[2] and params[q][3] == vals[3]:
                self.write_debug_message('xmetal, T, logg and micro have already been used in another iteration.')
                if self.atmos_params['change_antes'] == 'metallicity':
                    n = 'temperature'
                elif self.atmos_params['change_antes'] == 'temperature':
                    n = 'pressure'
                elif self.atmos_params['change_antes'] == 'pressure':
                    n = 'velocity'
                else:
                    n = 'metallicity'
                self.atmos_params['change'] = n
                self.nit += 1
                break

        del params, vals

    def check_nrepeat(self):
        if self.nit >= self.n_repeat:
            self.write_debug_message('Parameters were repeated more than %d times.' % (self.n_repeat))
            self.atmos_params['exception'] = 2
            return -1
        return 0


    def check_nit_total(self):
        if self.nit_total >= 500000:
            self.write_debug_message('More than 500000 iterations for the same star. Stopping the computation.')
            self.atmos_params['exception'] = 2
            return -1
        return 0


    def check_hold(self, xmetal):
        if self.atmos_params['temperature']['hold'] == 'yes':
            self.atmos_params['moog'][1] = 0.0
        if self.atmos_params['gravity']['hold'] == 'yes':
            self.atmos_params['moog'][2] = 0.0
        if self.atmos_params['metallicity']['hold'] == 'yes':
            self.atmos_params['moog'][0] = xmetal
        if self.atmos_params['velocity']['hold'] == 'yes':
            self.atmos_params['moog'][3] = 0.0

    def show_hold(self):
        hold_array = []
        for t in ['temperature', 'gravity', 'metallicity', 'velocity']:
            hold_array.append(self.atmos_params[t]['hold'])
        return hold_array


    def new_iteration(self, xmetal):
        self.nit_total += 1
        self.check_hold(xmetal)
        self.nout = 0
        self.add_param()
        self.atmos_params['change_antes'] = self.atmos_params['change']


    def moog_output(self, output, nfail):
        self.atmos_params['moog'] = output
        self.atmos_params['nfailed'] = nfail


    def new_values(self, new_vals):
        self.atmos_params['temperature']['value'] = new_vals[0]
        self.atmos_params['gravity']['value'] = new_vals[1]
        self.atmos_params['metallicity']['value'] = new_vals[2]
        self.atmos_params['velocity']['value'] = new_vals[3]
        #return self.atmos_params

    def next_change(self, change_ini):
        c = ['metallicity', 'temperature', 'pressure', 'velocity']
        i = c.index(change_ini)
        if i == 3: i = -1
        return c[i+1]

    def new_change(self, change_ini = None):
        if change_ini == None:
            change_ini = self.atmos_params['change']
        self.atmos_params['change'] = self.next_change(change_ini)

    def check_nfailed_it(self, change_ini):
        if self.atmos_params['nfailed'] > 0:
            self.new_change(change_ini)
            self.write_debug_message('Failed in metallicity. Change=%s' % (self.atmos_params['change']))
            self.nbreak += 1

    def change_metallicity(self, new_val):
        self.atmos_params['metallicity']['value'] = new_val

    def check_met(self):
        met = self.atmos_params['metallicity']['ranges']
        if (met[0] == met[1]) and (met[1] == met[2]) and (met[0] != self.atmos_params['metallicity']['value']):
            return True
        else:
            return False


    def change_parameter(self, name_par, moog_output, range_m, decimals):
        m_val = self.atmos_params['moog'][moog_output]
        ext = self.atmos_params[name_par]['ranges']
        val = self.atmos_params[name_par]['value']

        if m_val > 0.002:
            self.atmos_params[name_par]['ranges'][0] = val
            if val < ext[1]:
                self.atmos_params[name_par]['value'] = round(np.mean(ext), decimals)
            else:
                self.atmos_params[name_par]['value'] = round(mult(val, range_m, 'upper'), decimals)

        else:
            self.atmos_params[name_par]['ranges'][1] = val
            if ext[0] != -999. and val > ext[0]:
                self.atmos_params[name_par]['value'] = round(np.mean(ext), decimals)
            else:
                self.atmos_params[name_par]['value'] = round(mult(val, range_m, 'floor'), decimals)

        del m_val, ext, val


#******************************************************************************
#******************************************************************************
#******************************************************************************
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

    num = math.floor(x/base)
    if level == 'upper':
        final = (num + 1)*base
    else:
        final = num*base
        if final == x:
            final = final - base
    return final
