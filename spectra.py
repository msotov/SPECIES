from __future__ import print_function
import os
import glob
import numpy as np
from astropy.io import fits
from Instruments import feros, harps,\
                         uves, hires, spectra_sav,\
                         coralie, aat, pfs, other_instrument,\
                         lconres



def find_spectra(starname,\
                 use_HARPS=True, use_FEROS=True, use_UVES=True, use_HIRES=True, \
                 use_CORALIE=True, use_AAT=True, use_PFS=True, use_lconres=True,\
                 PATH_SPECTRA='./Spectra', restframe=True,
                 use_other=False, new_res=False, make_1d=False):

    print('\n\t%s' % starname)
    nfiles = 0
    name_instrument = []

    #############################################
    # Checks for spectra for the star
    #############################################

    spectra_feros = os.path.join(PATH_SPECTRA, '%s_feros.fits' % starname)
    spectra_feros_raw = os.path.join(PATH_SPECTRA, '%s_feros_o.fits' % starname)
    spectra_harps = os.path.join(PATH_SPECTRA, '%s_harps.fits' % starname)
    spectra_harps_dat = os.path.join(PATH_SPECTRA, '%s_harps.dat' % starname)
    spectra_uves = os.path.join(PATH_SPECTRA, '%s_uves.fits' % starname)
    spectra_uves_red = os.path.join(PATH_SPECTRA, '%s_uves_red.fits' % starname)
    spectra_uves_blue = os.path.join(PATH_SPECTRA, '%s_uves_blue.fits' % starname)
    spectra_hires = os.path.join(PATH_SPECTRA, '%s_hires/' % starname)
    spectra_HIRES = os.path.join(PATH_SPECTRA, '%s_HIRES.sav' % starname)
    spectra_coralie = os.path.join(PATH_SPECTRA, '%s_coralie.fits' % starname)
    spectra_aat = os.path.join(PATH_SPECTRA, '%s_aat.fits' % starname)
    spectra_pfs = os.path.join(PATH_SPECTRA, '%s_pfs.fits' % starname)
    spectra_pfs_sav = os.path.join(PATH_SPECTRA, '%s_pfs.sav' % starname)
    spectra_lconres = os.path.join(PATH_SPECTRA, '%s_LCONRES.fits' % starname)


    #############################################
    # FEROS spectra
    #############################################

    if use_FEROS:

        if os.path.isfile(spectra_feros):
            snr_feros = 0.0
            nfiles += 1
            print('\t\tFound FEROS spectra')

            starname_feros = '%s_feros' % starname
            append_inst = 0

            snr_feros = feros(os.path.join(PATH_SPECTRA, starname_feros),\
                              do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_feros, "__len__"):
                append_inst += 1
            elif snr_feros != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_feros)


        elif os.path.isfile(spectra_feros_raw):
            snr_feros = 0.0
            nfiles += 1
            print('\t\tFound FEROS spectra')

            starname_feros = '%s_feros_o' % starname
            append_inst = 0

            snr_feros = feros(os.path.join(PATH_SPECTRA, starname_feros),\
                              do_restframe=restframe, new_res=new_res, make_1d=make_1d)

            if hasattr(snr_feros, "__len__"):
                append_inst += 1
            elif snr_feros != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_feros)


    #############################################
    # HARPS spectra
    #############################################

    if use_HARPS:

        if os.path.isfile(spectra_harps) or os.path.isfile(spectra_harps_dat):
            snr_harps = 0.0
            nfiles += 1
            print('\t\tFound HARPS spectra')

            starname_harps = '%s_harps' % starname
            append_inst = 0

            snr_harps = harps(os.path.join(PATH_SPECTRA, starname_harps),\
                              do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_harps, "__len__"):
                append_inst += 1
            elif snr_harps != 0.:
                append_inst += 1
            if append_inst > 0:
                name_instrument.append(starname_harps)


    #############################################
    # UVES spectra
    #############################################

    if use_UVES:

        if os.path.isfile(spectra_uves_blue) and os.path.isfile(spectra_uves_red):
            snr_uves = np.zeros((2, 3))
            nfiles += 1
            print('\t\tFound UVES spectra')
            starname_uves = '%s_uves' % starname
            append_inst = 0

            snr_uves = uves(os.path.join(PATH_SPECTRA, starname_uves),\
                            do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_uves, "__len__"):
                append_inst += 1
            elif snr_uves != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_uves)

        elif os.path.isfile(spectra_uves):
            snr_uves = 0.0
            nfiles += 1
            print('\t\tFound UVES spectra')
            starname_uves = '%s_uves' % starname
            append_inst = 0

            snr_uves = uves(os.path.join(PATH_SPECTRA, starname_uves),\
                            do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if snr_uves > 0.0:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_uves)


    #############################################
    # HIRES spectra
    #############################################

    if use_HIRES:

        if os.path.exists(spectra_hires):
            snr_hires = 0.0
            nfiles += 1
            print('\t\tFound HIRES spectra')
            starname_hires = '%s_hires' % starname
            append_inst = 0

            snr_hires = hires(os.path.join(PATH_SPECTRA, starname_hires),\
                              do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_hires, "__len__"):
                append_inst += 1
            else:
                if snr_hires == 0.:
                    nfiles = nfiles - 1
                else:
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_hires)

        #############################################
        # .sav files
        #############################################

        if os.path.isfile(spectra_HIRES):
            snr_HIRES = 0.0
            nfiles += 1
            print('\t\tFound .sav files')

            starname_HIRES = '%s_HIRES' % starname
            append_inst = 0

            snr_HIRES = spectra_sav(os.path.join(PATH_SPECTRA, starname_HIRES),\
                                    do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_HIRES, "__len__"):
                append_inst += 1
            elif snr_HIRES != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_HIRES)


    #############################################
    # PFS spectra
    #############################################

    if use_PFS:

        if os.path.exists(spectra_pfs) or os.path.isfile(spectra_pfs_sav):
            snr_pfs = 0.0
            nfiles += 1
            print('\t\tFound PFS spectra')
            starname_pfs = '%s_pfs' % starname
            append_inst = 0

            snr_pfs = pfs(os.path.join(PATH_SPECTRA, starname_pfs),\
                          do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_pfs, "__len__"):
                append_inst += 1
            else:
                if snr_pfs == 0.:
                    nfiles = nfiles - 1
                else:
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_pfs)



    #############################################
    # CORALIE spectra
    #############################################

    if use_CORALIE:

        if os.path.isfile(spectra_coralie):
            snr_coralie = 0.0
            nfiles += 1
            print('\t\tFound CORALIE spectra')

            starname_coralie = '%s_coralie' % starname
            append_inst = 0

            snr_coralie = coralie(os.path.join(PATH_SPECTRA, starname_coralie),\
                                  do_restframe=restframe, new_res=new_res, make_1d=make_1d)

            if hasattr(snr_coralie, "__len__"):
                append_inst += 1
            elif snr_coralie != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_coralie)


    #############################################
    # AAT spectra
    #############################################

    if use_AAT:

        if os.path.isfile(spectra_aat):
            snr_aat = 0.0
            nfiles += 1
            print('\t\tFound AAT spectra')

            starname_aat = '%s_aat' % starname
            append_inst = 0

            snr_aat = aat(os.path.join(PATH_SPECTRA, starname_aat),\
                          do_restframe=restframe, new_res=new_res, make_1d=make_1d)

            if hasattr(snr_aat, "__len__"):
                append_inst += 1
            elif snr_aat != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_aat)


    #############################################
    # LCO-NRES spectra
    #############################################

    if use_lconres:

        if os.path.isfile(spectra_lconres):
            snr_lconres = 0.0
            nfiles += 1
            print('\t\tFound LCO NRES spectra')

            starname_lconres = '%s_LCONRES' % starname
            append_inst = 0

            snr_lconres = lconres(os.path.join(PATH_SPECTRA, starname_lconres),\
                                  do_restframe=restframe, new_res=new_res, make_1d=make_1d)

            if hasattr(snr_lconres, "__len__"):
                append_inst += 1
            elif snr_lconres != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_lconres)

    #############################################
    # Spectra from other instruments
    #############################################

    if use_other:

        list_files = glob.glob(os.path.join(PATH_SPECTRA, '%s_*.fits' % starname))
        to_remove = []
        for l in list_files:
            if any([x in l for x in ['feros', 'harps', 'uves', 'hires', 'coralie',
                                     'aat', 'pfs', '_res']]):
                to_remove.append(l)
        _ = [list_files.remove(l) for l in to_remove]

        for spectra_other in list_files:
            snr_other = 0.0
            nfiles += 1
            try:
                inst = fits.getval(spectra_other, 'INSTRUME', 0)
            except KeyError:
                try:
                    inst = fits.getval(spectra_other, 'INSTRUMENT', 0)
                except KeyError:
                    inst = 'other'
            print('\t\tFound spectra from %s' % inst)

            starname_other = spectra_other[spectra_other.index(starname):-5]
            if '_res' in starname_other:
                starname_other = starname_other[:starname_other.index('_res')]
            append_inst = 0

            snr_other = other_instrument(os.path.join(PATH_SPECTRA, starname_other),\
                                         do_restframe=restframe, new_res=new_res, make_1d=make_1d)
            if hasattr(snr_other, "__len__"):
                append_inst += 1
            elif snr_other != 0.:
                append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_other)


    if nfiles == 0:
        print('\n')
        print('\t\tNo spectra found for %s' % starname)
        print('\n')

    return name_instrument


#######################################################################
#######################################################################


def spectra(lista_estrellas,\
            HARPS=True, FEROS=True, UVES=True, HIRES=True, \
            CORALIE=True, AAT=True, PFS=True, LCONRES=True,\
            s_path='./Spectra', restframe=True, use_other=False, new_res=False,\
            make_1d=False):

    total_spectra = []

    for starname in lista_estrellas:
        name_instrument = find_spectra(starname, \
                                       use_HARPS=HARPS, \
                                       use_FEROS=FEROS, \
                                       use_UVES=UVES, \
                                       use_HIRES=HIRES, \
                                       use_CORALIE=CORALIE, \
                                       use_AAT=AAT, \
                                       use_PFS=PFS, \
                                       use_lconres=LCONRES, \
                                       PATH_SPECTRA=s_path, \
                                       restframe=restframe,\
                                       use_other=use_other,\
                                       new_res=new_res,
                                       make_1d=make_1d)

        total_spectra = total_spectra + name_instrument#_final

    return total_spectra
