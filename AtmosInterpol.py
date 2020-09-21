#interpol_new.py
from scipy.interpolate import griddata
import numpy as np
import os
import time
from AtmosGrid import grid

Tgg = grid['tgrid']
Ggg = grid['ggrid']
Mgg = grid['mgrid']
col0 = grid['col0']
col1 = grid['col1']
col2 = grid['col2']
col3 = grid['col3']
col4 = grid['col4']
col5 = grid['col5']
col6 = grid['col6']

def interpol(starname, xteff, xlogg, xmetal, micro, alias, fesun=7.50):
    ii = np.where((np.abs(Tgg-xteff) <= 301.) & (np.abs(Mgg-xmetal) <= 0.501) & (np.abs(Ggg-xlogg) <= 0.501))
    lT=np.shape(Tgg[ii])[0]
    lm=np.shape(Mgg[ii])[0]
    lg=np.shape(Ggg[ii])[0]
    test = False
    if all([lT > 1, lm > 1, lg > 1, micro <= 5., micro >= 0.]):
        dT = max(Tgg[ii])-min(Tgg[ii])
        dm = max(Mgg[ii])-min(Mgg[ii])
        dg = max(Ggg[ii])-min(Ggg[ii])
        tt = all([min(Tgg[ii]) <= xteff, max(Tgg[ii]) >= xteff, dT > 0.])
        tm = all([min(Mgg[ii]) <= xmetal, max(Mgg[ii]) >= xmetal, dm > 0.])
        tg = all([min(Ggg[ii]) <= xlogg, max(Ggg[ii]) >= xlogg, dg > 0.])
        if tt and tm and tg:
            test = True

    if not test:
        del ii
        return False

    gg = (Tgg[ii], Ggg[ii], Mgg[ii])
    xpoint = (xteff, xlogg, xmetal)
    col1g = griddata(gg, col1[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
    if np.isfinite(col1g[0]):
        col0g = griddata(gg, col0[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        col2g = griddata(gg, col2[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        col3g = griddata(gg, col3[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        col4g = griddata(gg, col4[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        col5g = griddata(gg, col5[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        col6g = griddata(gg, col6[ii], xpoint, method='linear', fill_value=np.nan, rescale=True)
        nlayers = len(col1g)
        with open('./atm_models/%s.atm' % alias, 'w') as model:
            model.write('KURUCZ\n')
            model.write('          Teff = {0:5.0f} log g = {1:4.2f}\n'.format(xteff, xlogg))
            model.write('NTAU        {0:3.0f}\n'.format(nlayers))
            for i in range(nlayers):
                model.write(' %14.8E  %7.1f %9.3E %9.3E %9.3E %9.3E %9.3E\n' % \
                                     (col0g[i], col1g[i], col2g[i], col3g[i],\
                                      col4g[i], col5g[i], col6g[i]))
            model.write('    {0:5.3f}e+05\n'.format(micro))
            model.write('NATOMS     1  {0:5.2f}\n'.format(xmetal))
            model.write('      26.0   {0:4.2f}\n'.format(xmetal+fesun))
            model.write('NMOL      19\n')
            model.write('      606.0    106.0    607.0    608.0    107.0    108.0    112.0    707.0\n')
            model.write('      708.0    808.0     12.1  60808.0  10108.0    101.0      6.1      7.1\n')
            model.write('        8.1    822.0     22.1\n')

        del gg, xpoint, col0g, col1g, col2g, col3g, col4g, col5g, col6g, ii, model
        return True
    del ii, lT, lm, lg, gg, xpoint, col1g
    return False
