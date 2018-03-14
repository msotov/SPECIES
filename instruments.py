import pyfits, os, glob
import numpy as np
from ccf import restframe
from PyAstronomy import pyasl
from scipy import interpolate
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import collections


def compute_snr(x, data):
    """
    Computes the S/N for spectra using a set of ranges where there shouldn't be
    any lines, only continuum.
    """
    ranges = [[5840.22, 5844.01], [5979.25, 5983.13], [6066.49, 6075.55], \
              [6195.95, 6198.74], [6386.36, 6392.14], [5257.83, 5258.58], \
              [5256.03, 5256.7]]
    sn = []
    for r in range(len(ranges)):
        data_range = data[np.where((x >= ranges[r][0]) & (x <= ranges[r][1]))[0]]
        if len(data_range) != 0:
            m = np.median(data_range)
            s = np.std(data_range)
            sn.append(m/s)

            del m, s

        del data_range

    del ranges

    return np.median(sn)


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def values(h, j):# funcion que extrae las longitudes de onda del header
    N = h['NAXIS' + str(j)]
    CRPIX = float(h['CRPIX' + str(j)])
    CDELT = float(h['CDELT' + str(j)])
    CRVAL = float(h['CRVAL' + str(j)])
    val = np.zeros(N);
    for i in range(N):
        val[i] = (i + 1 - CRPIX)*CDELT + CRVAL

    del N, CRPIX, CDELT, CRVAL
    return val


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************

def harps(starname, abundances = False):

    hdu = pyfits.open(starname + '.fits')
    header0 = hdu[0].header

    snr = 0.0

    if hdu[0].data is None:
        xnew = hdu[1].data['WAVE'][0]
        ynew = hdu[1].data['FLUX'][0]
        dif = [xnew[i+1] - xnew[i] for i in range(len(xnew)-1)]
        delta_x = np.mean(dif)

        if ('HIERARCH ESO DRS SPE EXT SN0' in header0):
            sn = [header0['HIERARCH ESO DRS SPE EXT SN' + str(i)] for i in range(72)]
            snr = np.median(sn)
            del sn
        hdu.close()


        # Create the new image

        os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
        os.system('rm -f ' + starname + '.fits')

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

        if 'RA' and 'DEC' in header0.keys():
            header_new['RA'] = (header0['RA'], '')
            header_new['DEC'] = (header0['DEC'], '')

        header_new.add_comment('Image created after combining orders')

        hdu.flush()
        hdu.close()

        del hdu, header0, xnew, ynew, dif, delta_x, header_new


    else:


        if len(hdu[0].data.shape) == 1:

            if ('HIERARCH ESO DRS SPE EXT SN0' in header0):
                sn = [header0['HIERARCH ESO DRS SPE EXT SN' + str(i)] for i in range(72)]
                snr = np.median(sn)
                del sn
            elif ('SNR' in header0):
                snr = header0['SNR']
            hdu.close()

            del hdu, header0

        elif len(hdu[0].data.shape) == 2 and hdu[0].data.shape[0] == 2:

            xnew_0 = hdu[0].data[0]
            ynew_0 = hdu[0].data[1]

            hdu.close()

            inan = np.where(np.isnan(ynew_0) == False)[0]

            tck = interpolate.InterpolatedUnivariateSpline(xnew_0[inan], ynew_0[inan], k = 5)

            del inan

            delta_av = np.median([xnew_0[i+1] - xnew_0[i] for i in range(len(xnew_0)-1)])

            xnew = np.arange(xnew_0[0], xnew_0[-1], delta_av)
            ynew = tck.__call__(xnew)

            os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
            os.system('rm -f ' + starname + '.fits')

            pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

            hdu = pyfits.open(starname + '.fits', mode = 'update')
            header_new = hdu[0].header

            header_new['CRPIX1'] = (1., 'Reference pixel')
            header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
            header_new['CDELT1'] = (delta_av, 'Coordinate increment per pixel')

            if 'RA' and 'DEC' in header0.keys():
                header_new['RA'] = (header0['RA'], '')
                header_new['DEC'] = (header0['DEC'], '')

            header_new.add_comment('Image created after combining orders')

            hdu.flush()
            hdu.close()

            if ('HIERARCH ESO DRS SPE EXT SN0' in header0):
                sn = [header0['HIERARCH ESO DRS SPE EXT SN' + str(i)] for i in range(72)]
                snr = np.median(sn)
                del sn
            elif ('SNR' in header0):
                snr = header0['SNR']

            del hdu, header0

        else:

            rangos_w = []
            tcks = []
            deltas = []

            rangos = []

            data = hdu[0].data[3]
            x = hdu[0].data[0]

            for i in range(x.shape[0]):
                data_i = data[i]
                x_i = x[i]

                rangos_w.append([x_i[0], x_i[-1]])
                deltas.append(x_i[1] - x_i[0])
                tck_flux = interpolate.InterpolatedUnivariateSpline(x_i, data_i, k = 5)
                tcks.append(tck_flux)

                rangos.append(x_i[0])
                rangos.append(x_i[-1])

                del data_i, x_i, tck_flux

            del data, x

            hdu.close()

            delta_x = max(deltas)
            xmin = min(rangos)
            xmax = max(rangos)
            Nnew = abs(xmin - xmax)/delta_x
            if round(Nnew) < Nnew: Nnew = Nnew + 1
            Nnew = int(Nnew)

            xnew = np.linspace(xmin, xmax, Nnew)
            ynew = np.zeros(len(xnew))

            deltas.reverse()
            rangos_w.reverse()
            tcks.reverse()

            for p in range(len(deltas)):
                xmin_p0 = rangos_w[p][0]
                xmax_p0 = rangos_w[p][1]

                indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
                tck = tcks[p]
                ynew[indices] = tck.__call__(xnew[indices])

                del xmin_p0, xmax_p0, indices, tck


            # Creates the fits file with the 1D spectra

            os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
            os.system('rm -f ' + starname + '.fits')

            pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

            hdu = pyfits.open(starname + '.fits', mode = 'update')
            header_new = hdu[0].header

            header_new['CRPIX1'] = (1., 'Reference pixel')
            header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
            header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

            if 'RA' and 'DEC' in header0.keys():
                header_new['RA'] = (header0['RA'], '')
                header_new['DEC'] = (header0['DEC'], '')

            header_new.add_comment('Image created after combining orders')

            hdu.flush()
            hdu.close()

            del header0, rangos_w, tcks, deltas, rangos, xnew, ynew,\
                delta_x, xmin, xmax, Nnew, header_new, hdu


    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')

    ##############################################
    # Plot absorption lines
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        if snr == 0.0:
            snr = compute_snr(x, data)

        ab = abundances
        plot_lines(x, data, starname, abundances = ab, save_fig = True)

        del x, data, ab

        return snr

    else:
        return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def feros(starname, abundances = False):

    hdu = pyfits.open(starname + '.fits')
    header0 = hdu[0].header

    snr = 0.0

    if hdu[0].data is None:
        xnew = hdu[1].data['WAVE'][0]
        ynew = hdu[1].data['FLUX'][0]
        dif = [xnew[i+1] - xnew[i] for i in range(len(xnew)-1)]
        delta_x = np.mean(dif)

        if 'SNR' in header0.keys():
            snr = header0['SNR']
        hdu.close()

        # Create the new image

        os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
        os.system('rm -f ' + starname + '.fits')

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

        if 'RA' and 'DEC' in header0.keys():
            header_new['RA'] = (header0['RA'], '')
            header_new['DEC'] = (header0['DEC'], '')

        header_new.add_comment('Image created after combining orders')

        hdu.flush()
        hdu.close()

        del xnew, ynew, dif, delta_x, header_new

    else:
        if len(hdu[0].data.shape) > 2:
            snr = feros_orders(starname, abundances = abundances)
            hdu.close()
            del hdu, header0
            return snr
        hdu.close()

    del hdu, header0


    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')

    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        if snr == 0.0:
            snr = compute_snr(x, data)

        ab = abundances
        plot_lines(x, data, starname, abundances = ab)

        del x, data, ab

        return snr

    else:
        return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def feros_orders(starname, abundances = False, do_restframe = True):

    hdulist = pyfits.open(starname + '.fits')

    header0 = hdulist[0].header

    is_raw_image = True
    if 'COMMENT' in header0.keys():
        if header0['COMMENT'] == 'Image created after combining orders':
            is_raw_image = False

    if is_raw_image:

        rangos_w = []
        tcks = []
        deltas = []

        rangos = []

        if ('PIPELINE' in header0.keys()) or ('NAXIS3' in header0.keys()):

            if header0['NAXIS3'] > 2:

                data = hdulist[0].data[5]
                x = hdulist[0].data[0]
                print 'hola'

            else:
                data = hdulist[0].data[1]
                x = hdulist[0].data[0]

            for i in range(x.shape[0]):
                data_i = data[i]
                x_i = x[i]

                rangos_w.append([x_i[0], x_i[-1]])
                deltas.append(x_i[1] - x_i[0])
                tck_flux = interpolate.InterpolatedUnivariateSpline(x_i, data_i, k = 5)
                tcks.append(tck_flux)

                rangos.append(x_i[0])
                rangos.append(x_i[-1])

                del data_i, x_i, tck_flux

            del data, x

        else:

            for i in range(39):
                header = hdulist[i].header
                data = hdulist[i].data
                x = values(header, 1)

                rangos_w.append([x[0], x[-1]])
                deltas.append(x[1] - x[0])
                tck_flux = interpolate.InterpolatedUnivariateSpline(x, data, k = 5)
                tcks.append(tck_flux)

                rangos.append(x[0])
                rangos.append(x[-1])

                del header, data, x, tck_flux

        hdulist.close()

        delta_x = max(deltas)
        xmin = min(rangos)
        xmax = max(rangos)
        Nnew = abs(xmin - xmax)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew + 1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))

        deltas.reverse()
        rangos_w.reverse()
        tcks.reverse()

        for p in range(len(deltas)):
            xmin_p0 = rangos_w[p][0]
            xmax_p0 = rangos_w[p][1]

            indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
            tck = tcks[p]
            ynew[indices] = tck.__call__(xnew[indices])

            del xmin_p0, xmax_p0, indices, tck


        # Creates the fits file with the 1D spectra

        os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
        os.system('rm -f ' + starname + '.fits')

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

        if 'RA' and 'DEC' in header0.keys():
            header_new['RA'] = (header0['RA'], '')
            header_new['DEC'] = (header0['DEC'], '')

        header_new.add_comment('Image created after combining orders')

        hdu.flush()
        hdu.close()

        del header0, rangos_w, tcks, deltas, rangos, xnew, ynew,\
            delta_x, xmin, xmax, Nnew, header_new, hdu

    else:
        hdulist.close()
        del header0

    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        if do_restframe:
            restframe(starname + '.fits')
        else:
            os.system('cp ' + starname + '.fits ' + starname + '_res.fits')
            print '\t\tNo restframe correction will be performed.'

    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        snr = compute_snr(x, data)

        ab = abundances

        plot_lines(x, data, starname, abundances = ab)
        del x, data, ab

        return snr

    else:
        return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def aat(starname, abundances = False):

    hdulist = pyfits.open(starname + '.fits')
    header0 = hdulist[0].header

    if os.path.isfile(starname + '_original.fits') == False:

        data = hdulist[0].data

        wave = data[:52]
        flux = data[52:]
        hdulist.close()

        rangos_w = []
        tcks = []
        deltas = []
        rangos = []

        for i in range(len(wave)):
            x = wave[i]
            y = flux[i]

            rangos_w.append([x[0], x[-1]])
            dif = np.array([x[i+1]-x[i] for i in range(len(x)-1)])
            deltas.append(np.mean(dif))
            tck_flux = interpolate.InterpolatedUnivariateSpline(x, y, k = 5)
            tcks.append(tck_flux)

            rangos.append(x[0])
            rangos.append(x[-1])

            del x, y, dif, tck_flux

        delta_x = max(deltas)
        xmin = min(rangos)
        xmax = max(rangos)
        Nnew = abs(xmin - xmax)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew + 1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))

        for p in range(len(deltas)):
            xmin_p0 = rangos_w[p][0]
            xmax_p0 = rangos_w[p][1]

            indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
            tck = tcks[p]
            ynew[indices] = tck.__call__(xnew[indices])

            del indices, tck


        # Creates the fits file with the 1D spectra
        #------------------------------------------

        os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
        os.system('rm -f ' + starname + '.fits')

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')
        header_new.add_comment('Image created after combining orders')

        hdu.flush()
        hdu.close()

        del data, wave, flux, deltas, rangos, tcks, delta_x, xmin, xmax, Nnew,\
            xnew, ynew, header_new, rangos_w

    else:
        hdulist.close()


    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')


    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        ##############################################
        # Plots the lines
        ##############################################

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        snr = compute_snr(x, data)

        ab = abundances

        plot_lines(x, data, starname, abundances = ab)

        del x, data, ab

        return snr

    else:
        return 0.


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def uves(starname, abundances = False):

    ##############################################
    # Creates the 1D spectra fits file
    # Combining the red and the blue parts
    # of the star
    ##############################################

    archivo1 = starname + '_red.fits'
    archivo2 = starname + '_blue.fits'

    # Defines the red part
    #-------------------------------------------
    hdulist1 = pyfits.open(archivo1)

    header1 = hdulist1[1].header
    data1 = hdulist1[1].data
    x1 = data1.field('WAVE')[0]
    if ('FLUX' in data1.columns.names):
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
    hdulist2 = pyfits.open(archivo2)

    header2 = hdulist2[1].header
    data2 = hdulist2[1].data
    x2 = data2.field('WAVE')[0]
    if ('FLUX' in data2.columns.names):
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
    for i in range(len(x2)):
        if (xmin1 <= x2[i] <= xmax1) == False:
            x.append(x2[i])
            flux.append(flux2[i])
        else:
            if sn2 > sn1:
                x.append(x2[i])
                flux.append(flux2[i])
    for i in range(len(x1)):
        if (xmin2 <= x1[i] <= xmax2) == False:
            x.append(x1[i])
            flux.append(flux1[i])
        else:
            if sn1 > sn2:
                x.append(x1[i])
                flux.append(flux1[i])

    xmax = max(x)
    xmin = min(x)

    del header1, data1, flux1, header1_0,\
        header2, data2, flux2, header2_0,\
        hdulist1, hdulist2

    #--------------------------------------------

    snr = [[xmin2, xmax2, sn2], [xmin1, xmax1, sn1]]

    if os.path.isfile(starname + '.fits'):
        print '\t\tImage ' + starname + '.fits has already been created.'
        print '\t\tSkipping the next step.\n'

    else:

        arr = np.zeros((len(x), 2))
        for i in range(len(x)):
            arr[i][0] = x[i]
            arr[i][1] = flux[i]

        arr_sort = arr[np.argsort(arr[:,0])]
        x_sort = [arr[i][0] for i in range(len(x))]
        flux_sort = [arr[i][1] for i in range(len(x))]


        # Create new wavelength and flux arrays
        #---------------------------------------------
        deltas = []
        for i in range(len(x2) - 1):
            deltas.append(x2[i+1] - x2[i])
        mean2 = np.mean(deltas)

        deltas2 = []
        for i in range(len(x1) - 1):
            deltas2.append(x1[i+1] - x1[i])
        mean1 = np.mean(deltas2)

        mean2 = round(mean2, 11)
        mean1 = round(mean1, 11)

        delta_x = 0
        if mean2 > mean1: delta_x = mean2
        else: delta_x = mean1



        Nnew = (xmax - xmin)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew+1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))


        #----------------------------------------------

        # Defines the cero_flux array
        #----------------------------------------------
        flux = np.array(flux_sort)
        x = np.array(x_sort)
        x_flux0 = x[np.where(flux == 0.)]

        xini = x_flux0[0]
        cero_flux = []


        for j in range(1, len(x_flux0)-1):
            xfinal = x_flux0[j]
            if (x_flux0[j+1] - xfinal) > 10.:
                cero_flux.append([xini, xfinal])
                xini = x_flux0[j+1]
                p = int(np.where(x == xini)[0])
                if (x[p] - x[p-1])>10.:
                    cero_flux.append([xfinal, xini])

                del p, xini
            del xfinal

        cero_flux = sorted(cero_flux)

        #----------------------------------------------

        # Defines the non_cero_flux array
        #----------------------------------------------

        non_cero_flux = []
        if cero_flux[0][0] != x[0]:
            indice=int(np.where(x == cero_flux[1][0])[0])
            non_cero_flux.append([x[0], x[indice-1]])
            del indice
        for p in range(len(cero_flux)-1):
            indice_min = int(np.where(x == cero_flux[p][1])[0])
            indice_max = int(np.where(x == cero_flux[p+1][0])[0])
            if indice_min != indice_max:
                non_cero_flux.append([x[indice_min+1], x[indice_max-1]])
            del indice_min, indice_max
        if x[-1]!=cero_flux[-1][1]:
            indice = int(np.where(x == cero_flux[-1][1])[0])
            non_cero_flux.append([x[indice+1], x[-1]])
            del indice

        #----------------------------------------------

        # Replace the values given by the interpolation
        # in the ranges at which the flux is different
        # from cero.
        #----------------------------------------------

        for i in range(len(non_cero_flux)):
            x_1 = non_cero_flux[i][0]
            x_2 = non_cero_flux[i][1]
            indice1 = int(np.where(x == x_1)[0])
            indice2 = int(np.where(x == x_2)[0])
            rango_x = [x[j] for j in range(indice1,indice2+1)]
            rango_flux = [flux[j] for j in range(indice1,indice2+1)]
            tck = interpolate.InterpolatedUnivariateSpline(rango_x, rango_flux, k = 5)
            xnew_rango = []
            for p in range(len(xnew)):
                if (x_1 <= xnew[p] <= x_2) == True:
                    xnew_rango.append(xnew[p])
            indice1_rango = int(np.where(xnew == xnew_rango[0])[0])
            indice2_rango = int(np.where(xnew == xnew_rango[-1])[0])
            ynew_rango = interpolate.UnivariateSpline.__call__(tck, xnew_rango)
            for p in range(indice1_rango, indice2_rango+1):
                ynew[p] = ynew_rango[p - indice1_rango]

            del x_1, x_2, indice1, indice2, rango_x, rango_flux, tck,\
                xnew_rango, indice1_rango, indice2_rango, ynew_rango

        #-----------------------------------------------


        # Creates the fits file with the 1D spectra
        # including in the the header the
        # CRVAL1, CDELT1 and CRPIX1 cards
        # Also, creates the 'starname_full_res.fits'
        # file, which is a copy of the
        # 'starname_full.fits' file, but will be use
        # to edit with restframe correction.
        #-----------------------------------------------

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu_new = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu_new[0].header
        data = hdu_new[0].data

        header_new['CRPIX1'] = (1.,'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x,'Coordinate increment par pixel')
        hdu_new.flush()
        hdu_new.close()

        del arr, arr_sort, x_sort, flux_sort, deltas,\
            deltas2, mean1, mean2, delta_x, Nnew, xnew, ynew,\
            x_flux0, cero_flux, non_cero_flux

    del files, x, flux, x1, x2
    #-----------------------------------------------

    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')

    ##############################################
    # Checks the number of lines present in the spectra
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        ab = abundances
        plot_lines(x, data, starname, abundances = ab)

        del x, data, ab

        return snr

    else:
        snr[0][2] = 0.0
        snr[1][2] = 0.0

        return snr


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def hires(starname, abundances = False):
    """
    Combines the values for the tables of each ccd
    to create a 1D image. It then corrects the spectra
    to restframe, which will be the image to
    be used by ARES.
    """

    ###############################################
    # Checks if the 1D image has been created.
    ###############################################

    if os.path.isfile(starname + '.fits'):
        print '\t\tImage ' + starname + '.fits has already been created.'
        print '\t\tSkipping the next step.\n'

    else:
        path = starname + '/HIRES/extracted/tbl'
        folders = glob.glob(path + '/*')
        folders = sorted(folders)

        ###############################################
        # Looks for the .tbl files in each ccd, and saves
        # the names in files.
        ###############################################

        files = []
        borders_ccds = []
        j = 0
        for i in range(len(folders)):
            path1 = folders[i] + '/flux'
            files1 = glob.glob(path1 + '/' + '*.tbl')
            if len(files1) == 0:
                files1_un = glob.glob(path1 + '/' + '*.gz')
                if len(files1_un) == 0:
                    print '\t\tNo tables files in ccd' + folders[i][len(folders[i]) - 1]
                else:
                    for f in files1_un:
                        subprocess.call(['gunzip', f])
                    files1 = glob.glob(path1 + '/' + '*.tbl')
                del files1_un
            files1 = sorted(files1)
            for f in files1:
                files.append(f)

            borders_ccds.append([j, j + len(files1) - 1])
            j += len(files1)

            del files1

        ###############################################
        # Open each .tbl file, and saves the interpolation
        # of the flux and s/n for each range of wavelength.
        # It also saves the ranges of wavelength and the
        # delta between them.
        ###############################################

        j = 0

        tcks_flux = [] # Interpolation params for the flux
        tcks_sn = [] # Interpolation params for the S/N
        rangos_w = [] # Ranges of wavelength
        deltas = [] # Difference between two wavelengths


        rangos_ccds = []

        valores_x = np.array([])
        valores_flux = np.array([])

        if len(files) == 0:
            print '\t\tNo .tbl files for this star.'

        for archivo in files:
            j += 1
            tabla = np.genfromtxt(archivo, dtype = None, skip_header = 1,\
                                  usecols = (4, 5, 8), names = ('wave', 'flux', 'sn'))
            x = tabla['wave']
            flux = tabla['flux']
            sn = tabla['sn']

            valores_x = np.hstack([valores_x, x])
            valores_flux = np.hstack([valores_flux, flux])

            # Creates the interpolation parameters and saves them
            # in tcks_flux and tcks_sn

            if np.all(x == -1.) == False:


                tck_flux = interpolate.InterpolatedUnivariateSpline(x, flux, k = 5)
                tck_sn = interpolate.InterpolatedUnivariateSpline(x, sn, k = 5)

                tcks_flux.append(tck_flux)
                tcks_sn.append(tck_sn)
                rangos_w.append([x[0], x[-1]])

                d = np.array([x[i + 1] - x[i] for i in range(len(x) - 1)])
                deltas.append(np.mean(d))
                if np.mean(d) == 0.:
                    print '\t\tWave values are the same for every row.'
                    break


                for p in range(len(borders_ccds)):
                    if (j - 1) == borders_ccds[p][0]:
                        rangos_ccds.append(x[0])
                    if (j - 1) == borders_ccds[p][1]:
                        rangos_ccds.append(x[-1])

                del tck_flux, tck_sn, d


            del x, flux, sn, tabla

        ###############################################
        # Do not create the image if the wavelength
        # values are wrong.
        ###############################################
        deltas = np.array(deltas)
        if (np.mean(deltas) == 0.) or (len(files) == 0) or (len(deltas) == 0):
            snr = 0.
            print '\t\tWas not able to generate the fits file.'

            del borders_ccds, tcks_flux, tcks_sn, rangos_w, deltas,\
                rangos_ccds, valores_x, valores_flux, folders, files

            return snr

        ###############################################
        # If the wavelength values are correct,
        # continue with creating the image.
        ###############################################

        else:
            ccds = [[rangos_ccds[i], rangos_ccds[i + 1]] for i in range(len(rangos_ccds)) if i % 2 == 0]

            ###############################################
            # Set the delta, xmin, xmax and number of points
            # the new image will have
            ###############################################

            delta_x = np.mean(deltas)
            xmin = rangos_w[0][0]
            xmax = rangos_w[-1][1]
            Nnew = abs(xmin - xmax)/delta_x
            if round(Nnew) < Nnew: Nnew = Nnew + 1
            Nnew = int(Nnew)

            xnew = np.linspace(xmin, xmax, Nnew)
            ynew = np.zeros(len(xnew))
            snr = np.array([[0.,0.] for i in range(len(xnew))])

            ###############################################
            # For each wavelength that corresponds to more than
            # one order, choose the values that have a
            # higher S/N.
            ###############################################

            for i in range(len(xnew)):
                for p in range(len(deltas)):
                    # x can be in any order but the final one
                    if p < (len(deltas) - 1):
                        xmin_p0 = rangos_w[p][0]
                        xmax_p0 = rangos_w[p][1]
                        xmin_p1 = rangos_w[p + 1][0]
                        xmax_p1 = rangos_w[p + 1][1]

                        # x is in only in one order.
                        if (xmin_p0 <= xnew[i] < xmax_p0) == True and (xmin_p1 <= xnew[i] < xmax_p1) == False:
                            tck = tcks_flux[p]
                            ynew[i] = interpolate.UnivariateSpline.__call__(tck, xnew[i])
                            tck_sn1 = tcks_sn[p]
                            snr[i][1] = interpolate.UnivariateSpline.__call__(tck_sn1, xnew[i])
                            snr[i][0] = xnew[i]
                            break
                        # x is in more than one order
                        if (xmin_p0 <= xnew[i] < xmax_p0) == True and (xmin_p1 <= xnew[i] < xmax_p1) == True:
                            tck_sn1 = tcks_sn[p]
                            tck_sn2 = tcks_sn[p + 1]

                            sn1 = interpolate.UnivariateSpline.__call__(tck_sn1, xnew[i])
                            sn2 = interpolate.UnivariateSpline.__call__(tck_sn2, xnew[i])

                            if sn1 >= sn2:
                                tck = tcks_flux[p]
                                ynew[i] = interpolate.UnivariateSpline.__call__(tck, xnew[i])
                                snr[i][1] = sn1
                                snr[i][0] = xnew[i]
                            else:
                                tck = tcks_flux[p + 1]
                                ynew[i] = interpolate.UnivariateSpline.__call__(tck, xnew[i])
                                snr[i][1] = sn2
                                snr[i][0] = xnew[i]
                            break
                    # x is in the last order only
                    else:
                        if (rangos_w[p][0] <= xnew[i] <= rangos_w[p][1]) == True and (rangos_w[p - 1][0] <= xnew[i] < rangos_w[p - 1][1]) == False:
                            tck = tcks_flux[p]
                            ynew[i] = interpolate.UnivariateSpline.__call__(tck, xnew[i])
                            tck_sn1 = tcks_sn[p]
                            snr[i][1] = interpolate.UnivariateSpline.__call__(tck_sn1, xnew[i])
                            snr[i][0] = xnew[i]

            del borders_ccds, tcks_flux, tcks_sn, rangos_w, deltas,\
                rangos_ccds, valores_x, valores_flux, folders, files

            ###############################################
            # Creates the fits file with the 1D spectra
            ###############################################

            pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

            hdu = pyfits.open(starname + '.fits', mode = 'update')
            header_new = hdu[0].header

            header_new['CRPIX1'] = (1., 'Reference pixel')
            header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
            header_new['CDELT1'] = (delta_x, 'Coordinate increment par pixel')

            hdu.flush()
            hdu.close()

            del xnew, ynew, snr


    ###############################################
    # Correct to restframe using a ccf between the
    # spectra and a binary mask.
    ###############################################

    # Check that the .fits image was created

    if os.path.isfile(starname + '.fits') == False:
        print '\t\t' + starname + '.fits has not yet been created.'
        print '\t\tThere must have been a problem when creating the image.'
        print '\t\tPlease check the tables files for ' + starname + '\n'
        snr = 0.
        return snr

    else:
        if os.path.isfile(starname + '_res.fits'):
            print '\t\tThere is already a file named ' + starname + '_res.fits'
        else:
            restframe(starname + '.fits')

        ###############################################
        # Checks the number of lines present in the spectra
        ###############################################

        if os.path.isfile(starname + '_res.fits'):

            x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

            ab = abundances
            plot_lines(x, data, starname, abundances = ab)

            range1 = '5302.94,5306.69'
            range2 = '5546.95,5553.64'
            range3 = '5721.31,5726.53'

            snr = []
            indice1 = range1.index(',')
            snr.append(float(range1[:indice1]))
            snr.append(float(range1[indice1 + 1:]))

            indice2 = range2.index(',')
            snr.append(float(range2[:indice2]))
            snr.append(float(range2[indice2 + 1:]))

            indice3 = range3.index(',')
            snr.append(float(range3[:indice3]))
            snr.append(float(range3[indice3 + 1:]))

            del x, data, ab, range1, range2, range3, \
                indice1, indice2, indice3

            return snr

        else:

            return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def spectra_sav(starname, abundances = False):

    if os.path.isfile(starname + '.fits') == False:

        data = readsav(starname + '.sav')

        wave = data.w
        flux = data.star

        rangos_w = []
        tcks = []
        deltas = []
        rangos = []


        for i in range(len(wave)):
            x = wave[i]
            y = flux[i]

            rangos_w.append([x[0], x[-1]])
            dif = np.array([x[i+1]-x[i] for i in range(len(x)-1)])
            deltas.append(np.mean(dif))
            tck_flux = interpolate.InterpolatedUnivariateSpline(x, y, k = 5)
            tcks.append(tck_flux)

            rangos.append(x[0])
            rangos.append(x[-1])

            del x, y, dif, tck_flux

        delta_x = max(deltas)
        xmin = min(rangos)
        xmax = max(rangos)
        Nnew = abs(xmin - xmax)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew + 1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))

        for p in range(len(deltas)):
            xmin_p0 = rangos_w[p][0]
            xmax_p0 = rangos_w[p][1]

            indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
            tck = tcks[p]
            ynew[indices] = tck.__call__(xnew[indices])

            del indices, xmin_p0, xmax_p0, tck


        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)
        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

        hdu.flush()
        hdu.close()

        del data, wave, flux, rangos_w, tcks, deltas, rangos,\
            xnew, ynew, delta_x, xmin, xmax, Nnew, hdu, header_new



    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')


    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        snr = compute_snr(x, data)

        ab = abundances

        plot_lines(x, data, starname, abundances = ab)

        del x, data, ab

        return snr

    else:
        return 0.0



#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def pfs(starname, abundances = False):

    if os.path.isfile(starname + '.fits') == False:

        data = readsav(starname + '.sav')

        wave = data.w
        flux = data.star

        rangos_w = []
        tcks = []
        deltas = []
        rangos = []


        for i in range(len(wave)):
            x = wave[i]
            y = flux[i]

            rangos_w.append([x[0], x[-1]])
            dif = np.array([x[i+1]-x[i] for i in range(len(x)-1)])
            deltas.append(np.mean(dif))
            tck_flux = interpolate.InterpolatedUnivariateSpline(x, y, k = 5)
            tcks.append(tck_flux)

            rangos.append(x[0])
            rangos.append(x[-1])

            del x, y, dif, tck_flux

        delta_x = max(deltas)
        xmin = min(rangos)
        xmax = max(rangos)
        Nnew = abs(xmin - xmax)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew + 1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))

        for p in range(len(deltas)):
            xmin_p0 = rangos_w[p][0]
            xmax_p0 = rangos_w[p][1]

            indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
            tck = tcks[p]
            ynew[indices] = tck.__call__(xnew[indices])

            del indices, xmin_p0, xmax_p0, tck


        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)
        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')

        hdu.flush()
        hdu.close()

        del data, wave, flux, rangos_w, tcks, deltas, rangos,\
            xnew, ynew, delta_x, xmin, xmax, Nnew, hdu, header_new



    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')


    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')

        snr = compute_snr(x, data)

        ab = abundances

        plot_lines(x, data, starname, abundances = ab)

        del x, data, ab

        return snr

    else:
        return 0.0



#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def coralie(starname, abundances = False):

    hdulist = pyfits.open(starname + '.fits')

    header0 = hdulist[0].header

    is_raw_image = True
    if 'COMMENT' in header0.keys():
        if header0['COMMENT'] == 'Image created after combining orders':
            is_raw_image = False


    if (is_raw_image == True) and (header0['NAXIS'] != 1):

        if header0['NAXIS'] == 3:
            if (header0['NAXIS3'] == 2):
                data = hdulist[0].data
                wave = data[0]
                flux = data[1]

            else:
                data = hdulist[0].data
                wave = data[0]
                flux = data[3]

        else:

            data = hdulist[0].data
            wave = data[0]
            flux = data[3]

        hdulist.close()

        rangos_w = []
        tcks = []
        deltas = []
        rangos = []

        for i in range(len(wave))[::-1]:
            x = wave[i]
            y = flux[i]

            rangos_w.append([x[0], x[-1]])
            dif = np.array([x[i+1]-x[i] for i in range(len(x)-1)])
            deltas.append(np.mean(dif))
            tck_flux = interpolate.InterpolatedUnivariateSpline(x, y, k = 5)
            tcks.append(tck_flux)

            rangos.append(x[0])
            rangos.append(x[-1])

        delta_x = max(deltas)
        xmin = min(rangos)
        xmax = max(rangos)
        Nnew = abs(xmin - xmax)/delta_x
        if round(Nnew) < Nnew: Nnew = Nnew + 1
        Nnew = int(Nnew)

        xnew = np.linspace(xmin, xmax, Nnew)
        ynew = np.zeros(len(xnew))

        for p in range(len(deltas)):
            xmin_p0 = rangos_w[p][0]
            xmax_p0 = rangos_w[p][1]

            indices = np.where((xnew >= xmin_p0) & (xnew < xmax_p0))[0]
            tck = tcks[p]
            ynew[indices] = tck.__call__(xnew[indices])

        del data, wave, flux, header0, rangos_w, tcks, deltas, rangos

        # Creates the fits file with the 1D spectra
        #------------------------------------------

        os.system('cp ' + starname + '.fits ' + starname + '_original.fits')
        os.system('rm -f ' + starname + '.fits')

        pyfits.writeto(starname + '.fits', data = ynew, clobber = True)

        hdu = pyfits.open(starname + '.fits', mode = 'update')
        header_new = hdu[0].header

        header_new['CRPIX1'] = (1., 'Reference pixel')
        header_new['CRVAL1'] = (xnew[0], 'Coordinate at reference pixel')
        header_new['CDELT1'] = (delta_x, 'Coordinate increment per pixel')
        header_new.add_comment('Image created after combining orders')

        hdu.flush()
        hdu.close()

        del xnew, ynew

    else:
        hdulist.close()
        del header0

    ##############################################
    # Correct to restframe
    ##############################################

    if os.path.isfile(starname + '_res.fits'):
        print '\t\tThere is already a file named ' + starname + '_res.fits'
    else:
        restframe(starname + '.fits')

    ##############################################
    # Computes the S/N
    ##############################################

    if os.path.isfile(starname + '_res.fits'):

        ##############################################
        # Plots the lines
        ##############################################

        x, data = pyasl.read1dFitsSpec(starname + '_res.fits')
        ab = abundances

        plot_lines(x, data, starname, abundances = ab)
        snr = compute_snr(x, data)

        del x, data, ab

        return snr

    else:
        return 0.0


#*****************************************************************************
#*****************************************************************************
#*****************************************************************************


def plot_lines(x, data, nombre, hold_plot = False, abundances = False, save_fig = True):
    if abundances == True:
        archivo = 'lines_ab.dat'
    else:
        archivo = 'linelist.dat'
    lines = np.genfromtxt('./Spectra/' + archivo, dtype = None, usecols = (0), skip_header = 2)

    fig, ax = plt.subplots()
    ax.plot(x, data)

    for l in lines:
        ax.axvline(x = l, color = 'red', linestyle = '--')

    ax.axvline(x = 6562.81, color = 'green', linestyle = '--')
    ax.axvline(x = (3933.664 + 0.545), color = 'green', linestyle = '--')
    ax.axvline(x = (3933.664 - 0.545), color = 'green', linestyle = '--')
    ax.axvline(x = (3968.470 + 0.545), color = 'green', linestyle = '--')
    ax.axvline(x = (3968.470 - 0.545), color = 'green', linestyle = '--')

    axins = inset_axes(ax, width = "25%", height = "25%", loc = 2)
    i_range  = np.where((x > 6550.) & (x < 6575.))[0]
    x_range = x[i_range]
    data_range = data[i_range]

    axins.plot(x_range, data_range)
    axins.axvline(x = 6562.81, color = 'green')
    plt.setp(axins.get_xticklabels(), visible = False)
    plt.setp(axins.get_yticklabels(), visible = False)

    plt.draw()

    if save_fig == True:
        if not os.path.exists('./EW/plots_spectra'):
            os.makedirs('./EW/plots_spectra')
        repeats = collections.Counter(nombre)['/']
        for p in range(repeats):
            indice = nombre.index('/')
            new_name = nombre[indice + 1:]
            nombre = new_name
        if abundances == False:
            fig.savefig('./EW/plots_spectra/' + nombre + '.pdf')
        else:
            fig.savefig('./EW/plots_spectra/' + nombre + '_ab.pdf')


    if hold_plot == True:
        press = raw_input('\t\tpress enter to continue')

    plt.close('all')

    del lines, fig, ax, axins, i_range, x_range, data_range
