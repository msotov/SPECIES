from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
import os
import glob
import subprocess
from past.utils import old_div
from astropy.io import fits
import numpy as np
from PyAstronomy import pyasl
from scipy import interpolate
from scipy.io.idl import readsav
from Ccf import restframe


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
    for r in ranges:
        data_range = data[np.where((x >= r[0]) & (x <= r[1]))[0]]
        if data_range.size > 4:
            m = np.median(data_range)
            #s1, _ = beq.betaSigma(data_range, 1, 1, returnMAD=True)
            #n = len(data_range)
            #s2 = 0.6052697*np.median(np.abs(2.0*data_range[2:n-2]-data_range[0:n-4]-data_range[4:n]))
            s = np.std(data_range)
            #print(m/s1, m/s2, m/s)
            sn.append(old_div(m, s))

            del m, s

        del data_range

    del ranges
    inonan = np.where(~np.isnan(sn))[0]
    return min(np.median(np.array(sn)[inonan]), 300)


#*****************************************************************************


def values(h, j):# funcion que extrae las longitudes de onda del header
    N = h['NAXIS' + str(j)]
    CRPIX = float(h['CRPIX' + str(j)])
    CDELT = float(h['CDELT' + str(j)])
    CRVAL = float(h['CRVAL' + str(j)])
    val = np.zeros(N)
    for i in range(N):
        val[i] = (i + 1 - CRPIX)*CDELT + CRVAL

    del N, CRPIX, CDELT, CRVAL
    return val


#*****************************************************************************

def create_new_wave_flux(xmin, xmax, N, rangos_w, deltas, tcks):
    xnew = np.linspace(xmin, xmax, N)
    ynew = np.zeros(N)

    for p, _ in enumerate(deltas):
        xmin_p0 = rangos_w[p][0]
        xmax_p0 = rangos_w[p][1]

        indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
        tck = tcks[p]
        ynew[indices] = tck(xnew[indices])
        del indices, xmin_p0, xmax_p0, tck

    return xnew, ynew

def combine_orders(wave, flux, reverse=False):
    rangos_w = []
    tcks = []
    deltas = []
    rangos = []

    for x, y in zip(wave, flux):
        inan = np.where((np.isnan(y) == False) & (x != 0.0))[0]
        x = x[inan]
        y = y[inan]
        del inan
        if x.size == 0:
            continue
        rangos_w.append([x[0], x[-1]])
        deltas.append(np.mean(x[1:] - x[:-1]))
        tck_flux = interpolate.UnivariateSpline(x, y, k=5, s=0)
        tcks.append(tck_flux)
        rangos.append(x[0])
        rangos.append(x[-1])

        del tck_flux

    Nnew = old_div(abs(min(rangos) - max(rangos)), max(deltas))
    Nnew = int(Nnew) if round(Nnew) >= Nnew else int(Nnew+1)

    if reverse:
        deltas.reverse()
        rangos_w.reverse()
        tcks.reverse()

    xnew, ynew = create_new_wave_flux(min(rangos), max(rangos), Nnew, rangos_w, deltas, tcks)
    del Nnew, rangos_w, deltas, tcks
    return xnew, ynew


#*****************************************************************************

def other_instrument(starname, do_restframe=True, new_res=False, make_1d=False):
    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    hdu = fits.open('%s.fits' % starname)
    header0 = hdu[0].header

    if 'SNR' in header0.keys():
        snr = header0['SNR']
    
    # Assuming there is only one header in the data
    if 'CRPIX1' in header0.keys() and 'CDELT1' in header0.keys() and 'CRVAL1' in header0.keys():
        wave, flux = pyasl.read1dFitsSpec('%s.fits' % starname)
    else:
        if header0['NAXIS'] == 2:
            wave = hdu[0].data[0]
            flux = hdu[0].data[1]

            inan = np.where(~np.isnan(flux))[0]
            wave = wave[inan]
            flux = flux[inan]
            del inan
        elif header0['NAXIS'] == 3 and header0['NAXIS1'] == 2:
            dataT = hdu[0].data.T
            wave = dataT[0].T
            flux = dataT[1].T
            del dataT

            #inan = np.where((~np.isnan(flux)) & (wave != 0.0))[0]
            #wave = wave[inan]
            #flux = flux[inan]

            wave_1d, flux_1d = combine_orders(wave, flux, reverse=True)
            #del inan

        elif header0['NAXIS'] == 3 and header0['NAXIS1'] > 2:
            wave = hdu[0].data[0]
            flux = hdu[0].data[1]

            #inan = np.where(~np.isnan(flux))[0]
            #wave = wave[inan]
            #flux = flux[inan]

            wave_1d, flux_1d = combine_orders(wave, flux, reverse=True)
            #del inan



    hdu.close()
    del hdu

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0
    if os.path.isfile(starname + '_res.fits'):
        return snr
    return 0.0



def harps(starname, do_restframe=True, new_res=False, make_1d=False):
    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile(starname + '_res.fits') and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    if os.path.isfile(starname + '.fits'):
        hdu = fits.open(starname + '.fits')
        header0 = hdu[0].header

        if hdu[0].data is None:
            try:
                wave = hdu[1].data['WAVE'][0]
                flux = hdu[1].data['FLUX'][0]
            except (ValueError, KeyError):
                wave = hdu[1].data['wav']
                flux = hdu[1].data['flux']

            if 'HIERARCH ESO DRS SPE EXT SN0' in header0:
                snr = np.median([header0['HIERARCH ESO DRS SPE EXT SN%d' % i] for i in range(72)])

            hdu.close()
            del hdu

        else:
            if len(hdu[0].data.shape) == 1:
                if 'ESO DRS SPE EXT SN0' in header0:
                    snr = np.median([header0['ESO DRS SPE EXT SN%d' % i] for i in range(72)])
                elif 'ESO DRS CAL EXT SN0' in header0:
                    snr = np.median([header0['ESO DRS CAL EXT SN%d' % i] for i in range(72)])
                elif 'SNR' in header0:
                    snr = header0['SNR']
                hdu.close()
                del hdu
                wave, flux = pyasl.read1dFitsSpec('%s.fits' % starname)

            elif (len(hdu[0].data.shape) == 2) and (hdu[0].data.shape[0] == 2):
                wave = hdu[0].data[0]
                flux = hdu[0].data[1]

                inan = np.where(~np.isnan(flux))[0]
                wave = wave[inan]
                flux = flux[inan]

                hdu.close()
                del hdu

                if 'HIERARCH ESO DRS SPE EXT SN0' in header0:
                    snr = np.median([header0['HIERARCH ESO DRS SPE EXT SN%d' % i]
                                     for i in range(72)])
                elif 'SNR' in header0:
                    snr = header0['SNR']

            elif (len(hdu[0].data.shape) == 2) and (hdu[0].data.shape[1] == 2):
                wave = hdu[0].data[:, 0]
                flux = hdu[0].data[:, 1]

                inan = np.where(~np.isnan(flux))[0]
                wave = wave[inan]
                flux = flux[inan]

                hdu.close()
                del hdu

                if 'HIERARCH ESO DRS SPE EXT SN0' in header0:
                    snr = np.median([header0['HIERARCH ESO DRS SPE EXT SN%d' % i]
                                     for i in range(72)])
                elif 'SNR' in header0:
                    snr = header0['SNR']

            else:
                flux = hdu[0].data[5]
                wave = hdu[0].data[0]
                sn = hdu[0].data[8]

                wave_1d, flux_1d = combine_orders(wave, flux, reverse=True)

                hdu.close()

                snr = [np.median(sno[np.where(sno > 0.0)[0]]) for sno in sn\
                        if np.where(sno > 0.0)[0].size > 0]
                del sn, hdu

    else:
        header0 = None

        wave, flux = np.loadtxt('%s.dat' % starname).T
        isort = np.argsort(wave)
        wave = wave[isort]
        flux = flux[isort]
        inan = np.where(~np.isnan(flux))[0]
        wave = wave[inan]
        flux = flux[inan]
        
        del isort, inan


    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0
    if os.path.isfile('%s_res.fits' % starname):
        return snr
    return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def feros(starname, do_restframe=True, new_res=False, make_1d=False):

    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named ' + starname + '_res.fits')
        return 100.0

    hdu = fits.open(starname + '.fits')
    header0 = hdu[0].header

    if hdu[0].data is None:

        wave = hdu[1].data['WAVE'][0]
        flux = hdu[1].data['FLUX'][0]

        inan = np.where(~np.isnan(flux))
        wave = wave[inan]
        flux = flux[inan]

        if 'SNR' in header0:
            snr = header0['SNR']
        hdu.close()
        del hdu

    else:
        if len(hdu[0].data.shape) >= 2:
            rangos_w = []
            tcks = []
            deltas = []

            rangos = []
            
            if ('NAXIS3' not in header0) and (header0['NAXIS2'] == 2):
                wave = hdu[0].data[0]
                flux = hdu[0].data[1]
                hdu.close()
            
            else:
                if ('PIPELINE' in header0) or ('NAXIS3' in header0):
                    if header0['NAXIS3'] > 2:
                        wave = hdu[0].data[0]
                        flux = hdu[0].data[5]
                        sn = hdu[0].data[8]

                        data = hdu[0].data[5]
                        x = hdu[0].data[0]
                        #print(np.unique(flux))
                        if len(np.unique(flux)) <= 1:
                            flux = hdu[0].data[1]
                            data = hdu[0].data[1]

                        snr = np.array([np.median(sno[np.where(sno > 0.0)[0]]) for sno in sn\
                                                    if np.where(sno > 0.0)[0].size > 0])

                    else:
                        wave = hdu[0].data[0]
                        flux = hdu[0].data[1]
    
                        data = hdu[0].data[1]
                        x = hdu[0].data[0]

                    for i in range(x.shape[0]):
                        data_i = data[i]
                        x_i = x[i]

                        rangos_w.append([x_i[0], x_i[-1]])
                        deltas.append(x_i[1] - x_i[0])
                        tck_flux = interpolate.InterpolatedUnivariateSpline(x_i, data_i, k=5)
                        tcks.append(tck_flux)

                        rangos.append(x_i[0])
                        rangos.append(x_i[-1])

                        del data_i, x_i, tck_flux

                    del data, x

                else:
                    wave = np.array([])
                    flux = np.array([])

                    for i in range(39):
                        header = hdu[i].header
                        data = hdu[i].data
                        x = values(header, 1)

                        rangos_w.append([x[0], x[-1]])
                        deltas.append(x[1] - x[0])
                        tck_flux = interpolate.InterpolatedUnivariateSpline(x, data, k=5)
                        tcks.append(tck_flux)

                        rangos.append(x[0])
                        rangos.append(x[-1])

                        wave = np.append(wave, x)
                        flux = np.append(flux, data)

                        del header, data, x, tck_flux

                hdu.close()
                del hdu

                delta_x = max(deltas)
                xmin = min(rangos)
                xmax = max(rangos)
                Nnew = old_div(abs(xmin - xmax), delta_x)
                Nnew = int(Nnew) if round(Nnew) >= Nnew else int(Nnew+1)

                deltas.reverse()
                rangos_w.reverse()
                tcks.reverse()

                wave_1d, flux_1d = create_new_wave_flux(xmin, xmax, Nnew, rangos_w, deltas, tcks)

        else:
            try:
                wave, flux = pyasl.read1dFitsSpec('%s.fits' % starname)
            except TypeError:
                crpix1 = float(header0['CRPIX1'])
                cdelt1 = float(header0['CDELT1'])
                crval1 = float(header0['CRVAL1'])

                N = int(header0['NAXIS1'])
                wave = ((np.arange(N) + 1.0) - crpix1) * cdelt1 + crval1
                flux = np.array(hdu[0].data)

                header0['CRPIX1'] = crpix1
                header0['CDELT1'] = cdelt1
                header0['CRVAL1'] = crval1

            hdu.close()
            del hdu

    if (snr is None) and ('SNR' in header0):
        snr = header0['SNR']

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************


def aat(starname, do_restframe=True, new_res=False, make_1d=False):

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    hdulist = fits.open('%s.fits' % starname)
    header0 = hdulist[0].header

    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    data = hdulist[0].data
    wave = data[:52]
    flux = data[52:]

    hdulist.close()
    del hdulist

    wave_1d, flux_1d = combine_orders(wave, flux)

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, data, header0

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************

def lconres(starname, do_restframe=True, new_res=False, make_1d=False):
    
    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    hdulist = fits.open('%s.fits' % starname)
    header0 = hdulist[0].header

    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    wave = hdulist['WAVESPEC'].data
    flux = hdulist['SPECBLAZE'].data

    hdulist.close()
    del hdulist

    wave_1d, flux_1d = combine_orders(wave, flux)#, reverse=True)

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0



#*****************************************************************************


def uves(starname, do_restframe=True, new_res=False, make_1d=False):

    snr = None
    wave, flux, wave_1d, flux_1d = None, None, None, None
    header0 = None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0
        
    if os.path.isfile('%s.fits' % starname):
        archivo = starname + '.fits'
        if fits.getval(archivo, 'NAXIS', 0) > 1:
            snr = fits.getval(archivo, 'SNR', 0)
            header0 = fits.getheader(archivo, 0)
            data = fits.getdata(archivo, 1)
            wave = data[0][0]
            flux = data[0][3]
            wave_1d = np.copy(wave)
            flux_1d = np.copy(flux)
            del data
        else:
            d = fits.open(archivo)
            if len(d) == 2:
                snr = fits.getval(archivo, 'SNR', 0)
                header0 = fits.getheader(archivo, 0)
                data = fits.getdata(archivo, 1)
                wave = data['WAVE'][0]
                flux = data['FLUX'][0]
                wave_1d = np.copy(wave)
                flux_1d = np.copy(flux)
                del data
            else:
                wave, flux = pyasl.read1dFitsSpec(archivo)
                header0 = fits.getheader(archivo, 0)
                wave_1d = np.copy(wave)
                flux_1d = np.copy(flux)
            d.close()
            del d

        snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr,
                        do_restframe=do_restframe, make_1d=make_1d)
        del wave, flux, wave_1d, flux_1d, header0

        if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
            return snr
        return 0.0

    # Creates the 1D spectra fits file
    # Combining the red and the blue parts
    # of the star

    archivo1 = '%s_red.fits' % starname
    archivo2 = '%s_blue.fits' % starname

    # Defines the red part
    #-------------------------------------------
    hdulist1 = fits.open(archivo1)

    header1 = hdulist1[1].header
    data1 = hdulist1[1].data
    x1 = np.array(data1.field('WAVE')[0])
    if 'FLUX' in data1.columns.names:
        flux1 = data1.field('FLUX')[0]
    else:
        flux1 = data1.field('FLUX_REDUCED')[0]
    xmax1 = max(x1)
    xmin1 = min(x1)

    header1_0 = hdulist1[0].header
    sn1 = header1_0['SNR']

    hdulist1.close()

    #--------------------------------------------
    # Defines the blue part
    #--------------------------------------------
    hdulist2 = fits.open(archivo2)

    header2 = hdulist2[1].header
    data2 = hdulist2[1].data
    x2 = np.array(data2.field('WAVE')[0])
    if 'FLUX' in data2.columns.names:
        flux2 = data2.field('FLUX')[0]
    else:
        flux2 = data2.field('FLUX_REDUCED')[0]
    xmax2 = max(x2)
    xmin2 = min(x2)

    header2_0 = hdulist2[0].header
    sn2 = header2_0['SNR']

    hdulist2.close()

    #--------------------------------------------
    # Combines the red and blue part
    #--------------------------------------------
    x = []
    flux = []
    for i, x2_ in enumerate(x2):
        if ~(xmin1 <= x2_ <= xmax1):
            x.append(x2_)
            flux.append(flux2[i])
        else:
            if sn2 > sn1:
                x.append(x2_)
                flux.append(flux2[i])
    for i, x1_ in enumerate(x1):
        if ~(xmin2 <= x1_ <= xmax2):
            x.append(x1_)
            flux.append(flux1[i])
        else:
            if sn1 > sn2:
                x.append(x1_)
                flux.append(flux1[i])

    xmax = max(x)
    xmin = min(x)

    del header1, data1, header1_0, header2, data2, header2_0, hdulist1, hdulist2

    #--------------------------------------------

    isort = np.argsort(np.array(x))
    x = np.array(x)[isort]
    flux = np.array(flux)[isort]
    del isort

    # Create new wavelength and flux arrays
    #---------------------------------------------
    mean2 = np.mean(x2[1:]-x2[:-1])
    mean1 = np.mean(x1[1:]-x1[:-1])

    mean2 = round(mean2, 11)
    mean1 = round(mean1, 11)

    delta_x = mean2 if mean2 > mean1 else mean1
    Nnew = (xmax - xmin)/ delta_x
    Nnew = int(Nnew) if round(Nnew) >= Nnew else int(Nnew+1)

    xnew = np.linspace(xmin, xmax, Nnew)
    ynew = np.zeros(Nnew)

    #----------------------------------------------
    # Defines the cero_flux array
    #----------------------------------------------
    x_flux0 = x[np.where(flux == 0.)]

    xini = x_flux0[0]
    cero_flux = []

    for j in range(1, len(x_flux0)-1):
        xfinal = x_flux0[j]
        if (x_flux0[j+1] - xfinal) > 10.:
            cero_flux.append([xini, xfinal])
            xini = x_flux0[j+1]
            p = int(np.where(x == xini)[0])
            if (x[p] - x[p-1]) > 10.:
                cero_flux.append([xfinal, xini])

        del xfinal

    cero_flux = np.array(sorted(cero_flux))

    #----------------------------------------------
    # Defines the non_cero_flux array
    #----------------------------------------------

    non_cero_flux = []
    if cero_flux[0][0] != x[0]:
        indice = int(np.where(x == cero_flux[1][0])[0])
        non_cero_flux.append([x[0], x[indice-1]])
        del indice
    for p0, p1 in zip(cero_flux[:-1], cero_flux[1:]):
        indice_min = int(np.where(x == p0[1])[0])
        indice_max = int(np.where(x == p1[0])[0])
        if indice_min != indice_max:
            non_cero_flux.append([x[indice_min+1], x[indice_max-1]])
        del indice_min, indice_max
    if x[-1] != cero_flux[-1][1]:
        indice = int(np.where(x == cero_flux[-1][1])[0])
        non_cero_flux.append([x[indice+1], x[-1]])
        del indice
    non_cero_flux = np.array(non_cero_flux)

    #----------------------------------------------
    # Replace the values given by the interpolation
    # in the ranges at which the flux is different
    # from cero.
    #----------------------------------------------

    for non_cero in non_cero_flux:
        x_1 = non_cero[0]
        x_2 = non_cero[1]
        i1 = int(np.where(x == x_1)[0])
        i2 = int(np.where(x == x_2)[0])
        rango_x = x[i1: i2+1]
        rango_flux = flux[i1: i2+1]
        tck = interpolate.UnivariateSpline(rango_x, rango_flux, k=5, s=0)
        ii = np.where((xnew >= x_1) & (xnew <= x_2))[0]
        xnew_rango = xnew[ii]
        i1_rango = int(np.where(xnew == xnew_rango[0])[0])
        i2_rango = int(np.where(xnew == xnew_rango[-1])[0])
        ynew[ii[:i2_rango-i1_rango+1]] = tck(xnew_rango[:i2_rango-i1_rango+1])

        del x_1, x_2, i1, i2, rango_x, rango_flux, tck, xnew_rango, i1_rango, i2_rango, ii

    del mean1, mean2, delta_x, Nnew,\
        x_flux0, cero_flux, non_cero_flux

    #--------------------------------------------
    # Create an array containing the red and blue
    # parts, as well as the two different SNR
    #--------------------------------------------

    isize = max(len(x1), len(x2))
    wave = np.zeros((2, isize))
    flux = np.zeros((2, isize))

    wave[0, :len(x1)] = x1
    wave[1, :len(x2)] = x2
    flux[0, :len(x1)] = flux1
    flux[1, :len(x2)] = flux2

    wave_1d = xnew
    flux_1d = ynew
    snr = np.array([sn1, sn2])

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0
    del x1, x2, flux1, flux2, xnew, ynew

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************


def hires(starname, do_restframe=True, new_res=False, make_1d=False):

    """
    Combines the values for the tables of each ccd
    to create a 1D image.
    """

    snr = None
    header0 = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    path = '%s/HIRES/extracted/tbl' % starname
    folders = glob.glob(path + '/*')
    folders = sorted(folders)

    # Looks for the .tbl files in each ccd, and saves
    # the names in files.

    files = []
    borders_ccds = []
    j = 0
    for f in folders:
        path1 = '%s/flux' % f
        files1 = glob.glob('%s/*.tbl' % path1)
        if not files1:
            files1_un = glob.glob('%s/*.gz' % path1)
            if not files1_un:
                print('\t\tNo tables files in ccd %s' % f[-1])
            else:
                for f_ in files1_un:
                    subprocess.call(['gunzip', f_])
                files1 = glob.glob('%s/*.tbl' % path1)
            del files1_un
        files1 = sorted(files1)
        for f_ in files1:
            files.append(f_)

        borders_ccds.append([j, j + len(files1) - 1])
        j += len(files1)

        del files1

    # Open each .tbl file, and saves the interpolation
    # of the flux and s/n for each range of wavelength.
    # It also saves the ranges of wavelength and the
    # delta between them.

    tcks_flux = [] # Interpolation params for the flux
    tcks_sn = [] # Interpolation params for the S/N
    rangos_w = [] # Ranges of wavelength
    deltas = [] # Difference between two wavelengths

    rangos_ccds = []
    valores_x = []
    valores_flux = []
    valores_sn = []

    if not files:
        print('\t\tNo .tbl files for this star.')

    for j, archivo in enumerate(files):
        tabla = np.genfromtxt(archivo, dtype=None, skip_header=1,\
                              usecols=(4, 5, 8), names=('wave', 'flux', 'sn'))
        x = tabla['wave']
        flux = tabla['flux']
        sn = tabla['sn']


        valores_x.append(x)
        valores_flux.append(flux)
        valores_sn.append(sn)

        # Creates the interpolation parameters and saves them
        # in tcks_flux and tcks_sn

        if ~np.all(x == -1.):
            tck_flux = interpolate.InterpolatedUnivariateSpline(x, flux, k=5)
            tck_sn = interpolate.InterpolatedUnivariateSpline(x, sn, k=5)

            tcks_flux.append(tck_flux)
            tcks_sn.append(tck_sn)
            rangos_w.append([x[0], x[-1]])

            d = x[1:]-x[:-1]
            deltas.append(np.mean(d))
            if np.mean(d) == 0.:
                print('\t\tWave values are the same for every row.')
                break

            for border in borders_ccds:
                if j == border[0]:
                    rangos_ccds.append(x[0])
                if j == border[1]:
                    rangos_ccds.append(x[-1])

            del tck_flux, tck_sn, d

        del x, flux, sn, tabla

    # Do not create the image if the wavelength
    # values are wrong.
    deltas = np.array(deltas)
    if (np.mean(deltas) == 0.) or (not files) or deltas.size == 0:
        snr = 0.
        print('\t\tNo valid wavelength range. Skipping this image.')

        del borders_ccds, tcks_flux, tcks_sn, rangos_w, deltas,\
            rangos_ccds, valores_x, valores_flux, folders, files

        return snr

    # If the wavelength values are correct,
    # continue with creating the image.

    # Set the delta, xmin, xmax and number of points
    # the new image will have

    delta_x = np.mean(deltas)
    xmin = rangos_w[0][0]
    xmax = rangos_w[-1][1]
    Nnew = (xmax - xmin)/ delta_x
    Nnew = int(Nnew) if round(Nnew) >= Nnew else int(Nnew+1)

    xnew = np.linspace(xmin, xmax, Nnew)
    ynew = np.zeros(len(xnew))
    snr = np.zeros((xnew.size, 2))

    # For each wavelength that corresponds to more than
    # one order, choose the values that have a
    # higher S/N.

    for i, _ in enumerate(xnew):
        for p, _ in enumerate(deltas):
            # x can be in any order but the final one
            if p < (len(deltas) - 1):
                xmin_p0 = rangos_w[p][0]
                xmax_p0 = rangos_w[p][1]
                xmin_p1 = rangos_w[p + 1][0]
                xmax_p1 = rangos_w[p + 1][1]

                # x is in only in one order.
                if (xmin_p0 <= xnew[i] < xmax_p0) and ~(xmin_p1 <= xnew[i] < xmax_p1):
                    ynew[i] = tcks_flux[p](xnew[i])
                    snr[i][1] = tcks_sn[p](xnew[i])
                    snr[i][0] = xnew[i]
                    break
                # x is in more than one order
                if (xmin_p0 <= xnew[i] < xmax_p0) and (xmin_p1 <= xnew[i] < xmax_p1):
                    sn1 = tcks_sn[p](xnew[i])
                    sn2 = tcks_sn[p+1](xnew[i])

                    if sn1 >= sn2:
                        ynew[i] = tcks_flux[p](xnew[i])
                        snr[i][1] = sn1
                        snr[i][0] = xnew[i]
                    else:
                        ynew[i] = tcks_flux[p+1](xnew[i])
                        snr[i][1] = sn2
                        snr[i][0] = xnew[i]
                    break
            # x is in the last order only
            else:
                if (rangos_w[p][0] <= xnew[i] <= rangos_w[p][1]) and\
                    ~(rangos_w[p - 1][0] <= xnew[i] < rangos_w[p - 1][1]):
                    ynew[i] = tcks_flux[p](xnew[i])
                    snr[i][1] = tcks_sn[p](xnew[i])
                    snr[i][0] = xnew[i]

    del borders_ccds, tcks_flux, tcks_sn, rangos_w, deltas,\
        rangos_ccds, folders, files

    wave = np.copy(valores_x)
    flux = np.copy(valores_flux)

    snr = np.array([np.median(sno[np.where(sno > 0.0)[0]]) for sno in valores_sn])

    wave_1d = xnew[:]
    flux_1d = ynew[:]

    del valores_x, valores_flux, valores_sn, xnew, ynew

    fits.writeto('%s.fits' % starname, data=np.vstack((wave, flux)), overwrite=True)

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************

def spectra_sav(starname, do_restframe=True, new_res=False, make_1d=False):

    snr = None
    header0 = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    data = readsav('%s.sav' % starname)
    wave = data.w
    flux = data.star

    wave_1d, flux_1d = combine_orders(wave, flux)

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0, data

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************


def pfs(starname, do_restframe=True, new_res=False, make_1d=False):

    snr = None
    header0 = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    if ~os.path.isfile('%s.fits' % starname):

        data = readsav('%s.sav' % starname)
        wave = data.w
        flux = data.star

        wave_1d, flux_1d = combine_orders(wave, flux)
        del data

    else:
        wave, flux = pyasl.read1dFitsSpec('%s.fits' % starname)

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0


#*****************************************************************************


def coralie(starname, do_restframe=True, new_res=False, make_1d=False):

    snr = None
    header0 = None
    wave, flux, wave_1d, flux_1d = None, None, None, None

    if os.path.isfile('%s_res.fits' % starname) and not new_res:
        print('\t\tThere is already a file named %s_res.fits' % starname)
        return 100.0

    hdulist = fits.open('%s.fits' % starname)
    header0 = hdulist[0].header

    if header0['NAXIS'] != 1:
        if header0['NAXIS'] == 3:
            data = hdulist[0].data
            wave = data[0]
            if header0['NAXIS3'] == 2:
                flux = data[1]
            else:
                flux = data[5]

        else:
            data = hdulist[0].data
            wave = data[0]
            try:
                flux = data[3]
            except IndexError:
                flux = data[1]

        hdulist.close()

        if len(wave.shape) > 1:
            rangos_w = []
            tcks = []
            deltas = []
            rangos = []

            for i in range(len(wave))[::-1]:
                x = wave[i]
                y = flux[i]

                rangos_w.append([x[0], x[-1]])
                deltas.append(np.mean(x[1:]-x[:-1]))
                tck_flux = interpolate.InterpolatedUnivariateSpline(x, y, k=5)
                tcks.append(tck_flux)

                rangos.append(x[0])
                rangos.append(x[-1])

            delta_x = max(deltas)
            xmin = min(rangos)
            xmax = max(rangos)
            Nnew = old_div(abs(xmin - xmax), delta_x)
            Nnew = int(Nnew) if round(Nnew) >= Nnew else int(Nnew+1)

            xnew, ynew = create_new_wave_flux(xmin, xmax, Nnew, rangos_w, deltas, tcks)
            del rangos_w, tcks, deltas, rangos

        else:
            xnew = wave[:]
            ynew = flux[:]

        wave_1d = xnew[:]
        flux_1d = ynew[:]

        del data, xnew, ynew

    else:
        wave, flux = pyasl.read1dFitsSpec('%s.fits' % starname)
        hdulist.close()

    snr = restframe(starname, wave, flux, wave_1d, flux_1d, header0, snr, do_restframe=do_restframe,\
                    make_1d=make_1d)
    del wave, flux, wave_1d, flux_1d, header0, hdulist

    if os.path.isfile('%s_res.fits' % starname) and (snr is not None):
        return snr
    return 0.0
