import os, sys
import numpy as np
from instruments import feros, feros_orders, harps,\
                        uves, hires, spectra_sav,\
                        coralie, aat, pfs


#######################################################################
#######################################################################

def modify_file2(starname, instrument, snr, abundances = False, PATH_SPECTRA = './Spectra'):
    PATH_EW = '../EW/'

    archivo_ares = open('./Spectra/mine.opt', 'r')
    archivo_ares_output = open('./Spectra/mine2.opt', 'w')
    for linea in archivo_ares:
        if linea[:8] == 'specfits':
            linea = "specfits='%s_res.fits'\n" % (os.path.join(os.path.relpath(PATH_SPECTRA, './Spectra'), starname))
        if linea[:11] == 'readlinedat':
            if abundances == False:
                linea = "readlinedat='linelist.dat'\n"
            else:
                linea = "readlinedat='lines_ab.dat'\n"
        if linea[:7] == 'fileout':
            if abundances == False:
                linea = "fileout='%s'\n" % (os.path.join(PATH_EW, starname + '.ares'))
            else:
                linea = "fileout='%s'\n" % (os.path.join(PATH_EW, starname + '_ab.ares'))
        if linea[:4] == 'rejt':
            if instrument == 'uves':
                rej1 = 1.-1./(snr[0][2])
                rej2 = 1.-1./(snr[1][2])
                linea = "rejt=%f\n" % np.mean([rej1,rej2])
            elif instrument == 'hires':
                linea = "rejt=3;%f,%f,%f,%f,%f,%f\n" % (snr[0], snr[1], snr[2], snr[3], snr[4], snr[5])
            else:
                linea = "rejt=%f\n" % (1.-1./snr)
        archivo_ares_output.writelines(linea)
    archivo_ares.close()
    archivo_ares_output.close()

    os.system('cp ./Spectra/mine2.opt ./Spectra/mine.opt')
    os.system('rm -f ./Spectra/mine2.opt')

#######################################################################
#######################################################################


def modify_file(starname, instrument, snr, abundances = False, PATH_SPECTRA = './Spectra'):
    PATH_EW = '../EW/'
    output = open('./Spectra/mine2.opt', 'w')

    output.writelines("specfits='%s_res.fits'\n" % (os.path.join(os.path.relpath(PATH_SPECTRA, './Spectra'), starname)))
    if abundances == False:
        output.writelines("readlinedat='linelist.dat'\n")
        output.writelines("fileout='%s'\n" % (os.path.join(PATH_EW, starname + '.ares')))
    else:
        output.writelines("readlinedat='lines_ab.dat'\n")
        output.writelines("fileout='%s'\n" % (os.path.join(PATH_EW, starname + '_ab.ares')))
    output.writelines("lambdai=3600.\n")
    output.writelines("lambdaf=9000.\n")
    output.writelines("smoothder=4\n")
    output.writelines("space=3.0\n")
    if instrument == 'uves':
        rej1 = 1.-1./(snr[0][2])
        rej2 = 1.-1./(snr[1][2])
        output.writelines("rejt=%f\n" % np.mean([rej1,rej2]))
    elif instrument == 'hires':
        output.writelines("rejt=3;%f,%f,%f,%f,%f,%f\n" % (snr[0], snr[1], snr[2], snr[3], snr[4], snr[5]))
    else:
        output.writelines("rejt=%f\n" % (1.-1./snr))
    output.writelines("lineresol=0.1\n")
    output.writelines("miniline=2\n")
    output.writelines("plots_flag=0\n")

    output.close()
    os.system('cp ./Spectra/mine2.opt ./Spectra/mine.opt')
    os.system('rm -f ./Spectra/mine2.opt')

    del output


#######################################################################
#######################################################################

def modify_file_vsini(starname, instrument, snr, PATH_SPECTRA = './Spectra'):
    PATH_EW = '../EW/'
    output = open('./Spectra/mine2.opt', 'w')

    output.writelines("specfits='%s_res.fits'\n" % (os.path.join(os.path.relpath(PATH_SPECTRA, './Spectra'), starname)))
    output.writelines("readlinedat='linelist_vsini.dat'\n")
    output.writelines("fileout='%s'\n" % (os.path.join(PATH_EW, starname + '_vsini.ares')))
    output.writelines("lambdai=3600.\n")
    output.writelines("lambdaf=9000.\n")
    output.writelines("smoothder=4\n")
    output.writelines("space=3.0\n")
    if instrument == 'uves':
        rej1 = 1.-1./(snr[0][2])
        rej2 = 1.-1./(snr[1][2])
        output.writelines("rejt=%f\n" % np.mean([rej1,rej2]))
    elif instrument == 'hires':
        output.writelines("rejt=3;%f,%f,%f,%f,%f,%f\n" % (snr[0], snr[1], snr[2], snr[3], snr[4], snr[5]))
    else:
        output.writelines("rejt=%f\n" % (1.-1./snr))
    output.writelines("lineresol=0.1\n")
    output.writelines("miniline=2\n")
    output.writelines("plots_flag=0\n")

    output.close()
    os.system('cp ./Spectra/mine2.opt ./Spectra/mine.opt')
    os.system('rm -f ./Spectra/mine2.opt')

    del output


#######################################################################
#######################################################################


def compute_ew(starname, skip = True, abundances = False, \
                use_HARPS = True, use_FEROS = True, use_UVES = True, use_HIRES = True, \
                use_CORALIE = True, use_AAT = True, use_PFS = True,\
                vsini = False, PATH_SPECTRA = './Spectra'):

    print '\n\t' + starname

    nfiles = 0

    name_instrument = []

    #############################################
    # Checks for spectra for the star
    #############################################

    #PATH_SPECTRA = './Spectra/'
    PATH_EW = './EW/'

    spectra_feros = os.path.join(PATH_SPECTRA, starname + '_feros.fits')
    spectra_feros_raw = os.path.join(PATH_SPECTRA, starname + '_feros_o.fits')
    spectra_harps = os.path.join(PATH_SPECTRA, starname + '_harps.fits')
    spectra_harps_dat = os.path.join(PATH_SPECTRA, starname + '_harps.dat')
    spectra_uves_red = os.path.join(PATH_SPECTRA, starname + '_uves_red.fits')
    spectra_uves_blue = os.path.join(PATH_SPECTRA, starname + '_uves_blue.fits')
    spectra_hires = os.path.join(PATH_SPECTRA, starname + '_hires/')
    spectra_HIRES = os.path.join(PATH_SPECTRA, starname + '_HIRES.sav')
    spectra_coralie = os.path.join(PATH_SPECTRA, starname + '_coralie.fits')
    spectra_aat = os.path.join(PATH_SPECTRA, starname + '_aat.fits')
    spectra_pfs = os.path.join(PATH_SPECTRA, starname + '_pfs.fits')
    spectra_pfs_sav = os.path.join(PATH_SPECTRA, starname + '_pfs.sav')

    #print spectra_harps


    #############################################
    # FEROS spectra
    #############################################

    if use_FEROS == True:

        if os.path.isfile(spectra_feros):
            snr_feros = 0.0
            nfiles += 1
            print '\t\tFound FEROS spectra'

            starname_feros = starname + '_feros'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_feros + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_feros = feros(os.path.join(PATH_SPECTRA, starname_feros), abundances = abundances)

                if snr_feros != 0.:

                    modify_file(starname_feros, 'feros', snr_feros, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_feros)


            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_feros + '_vsini.ares'
                else:
                    if snr_feros != 0.:
                        modify_file_vsini(starname_feros, 'feros', snr_feros, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


        #############################################
        # FEROS spectra raw
        #############################################

        if (os.path.isfile(spectra_feros_raw)==True) and (os.path.isfile(spectra_feros)==False):
            snr_feros_orders = 0.0
            nfiles += 1
            print '\t\tFound FEROS spectra separated by orders'

            starname_feros_o = starname + '_feros_o'
            starname_feros = starname + '_feros'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros_o + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros_o + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_feros_o + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_feros_orders = feros_orders(os.path.join(PATH_SPECTRA, starname_feros_o), abundances = abundances)

                if snr_feros_orders != 0.:

                    modify_file(starname_feros_o, 'feros_o', snr_feros_orders, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_feros_o)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_feros_o + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_feros_o + '_vsini.ares'
                else:
                    if snr_feros_orders != 0.:
                        modify_file_vsini(starname_feros_o, 'feros_o', snr_feros_orders, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')

    #############################################
    # HARPS spectra
    #############################################

    if use_HARPS == True:

        if os.path.isfile(spectra_harps) or os.path.isfile(spectra_harps_dat):
            snr_harps = 0.0
            nfiles += 1
            print '\t\tFound HARPS spectra'

            starname_harps = starname + '_harps'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_harps + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_harps + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_harps + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_harps = harps(os.path.join(PATH_SPECTRA, starname_harps), abundances = abundances)

                if snr_harps != 0.:

                    modify_file(starname_harps, 'harps', snr_harps, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')

                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_harps)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_harps + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_harps + '_vsini.ares'
                else:
                    if snr_harps != 0.:
                        modify_file_vsini(starname_harps, 'harps', snr_harps, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')

                    if os.path.isfile(os.path.join(PATH_EW, starname_harps + '_vsini.ares')) == False:
                        print '\t\tCould not create file ' + starname_harps + '_vsini.ares'


    #############################################
    # UVES spectra
    #############################################

    if use_UVES == True:

        if os.path.isfile(spectra_uves_blue) and os.path.isfile(spectra_uves_red):
            snr_uves = np.zeros((2,3))
            nfiles += 1
            print '\t\tFound UVES spectra'
            starname_uves = starname + '_uves'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_uves + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_uves + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_uves + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_uves = uves(os.path.join(PATH_SPECTRA, starname_uves), abundances = abundances)

                if (snr_uves[0][2] != 0.0) and (snr_uves[1][2] != 0.0):
                    modify_file(starname_uves, 'uves', snr_uves, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_uves)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW + starname_uves + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_uves + '_vsini.ares'
                else:
                    if (snr_uves[0][2] != 0.) and (snr_uves[1][2] != 0.0):
                        modify_file_vsini(starname_uves, 'uves', snr_uves, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


    #############################################
    # HIRES spectra
    #############################################

    if use_HIRES == True:

        if os.path.exists(spectra_hires):
            snr_hires = 0.0
            nfiles += 1
            print '\t\tFound HIRES spectra'
            starname_hires = starname + '_hires'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_hires + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_hires + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_hires + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_hires = hires(os.path.join(PATH_SPECTRA, starname_hires), abundances = abundances)
                if snr_hires == 0.:
                    nfiles = nfiles - 1
                else:
                    modify_file(starname_hires, 'hires', snr_hires, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_hires)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_hires + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_hires + '_vsini.ares'
                else:
                    if snr_hires != 0.:
                        modify_file_vsini(starname_hires, 'hires', snr_hires, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')

        #############################################
        # .sav files
        #############################################

        if os.path.isfile(spectra_HIRES):
            snr_HIRES = 0.0
            nfiles += 1
            print '\t\tFound .sav files'

            starname_HIRES = starname + '_HIRES'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_HIRES + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_HIRES + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_HIRES + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_HIRES = spectra_sav(os.path.join(PATH_SPECTRA, starname_HIRES), abundances = abundances)

                if snr_HIRES != 0.:
                    modify_file(starname_HIRES, 'sav', snr_HIRES, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_HIRES)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_HIRES + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_HIRES + '_vsini.ares'
                else:
                    if snr_HIRES != 0.:
                        modify_file_vsini(starname_HIRES, 'sav', snr_HIRES, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


    #############################################
    # PFS spectra
    #############################################

    if use_PFS:

        if os.path.exists(spectra_pfs) or os.path.isfile(spectra_pfs_sav):
            snr_pfs = 0.0
            nfiles += 1
            print '\t\tFound PFS spectra'
            starname_pfs = starname + '_pfs'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_pfs + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_pfs + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_pfs + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_pfs = pfs(os.path.join(PATH_SPECTRA, starname_pfs), abundances = abundances)
                if snr_pfs == 0.:
                    nfiles = nfiles - 1
                else:
                    modify_file(starname_pfs, 'pfs', snr_pfs, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_pfs)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_pfs + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_pfs + '_vsini.ares'
                else:
                    if snr_pfs != 0.:
                        modify_file_vsini(starname_pfs, 'pfs', snr_pfs, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


    #############################################
    # CORALIE spectra
    #############################################

    if use_CORALIE == True:

        if os.path.isfile(spectra_coralie):
            snr_coralie = 0.0
            nfiles += 1
            print '\t\tFound CORALIE spectra'

            starname_coralie = starname + '_coralie'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_coralie + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_coralie + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_coralie + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_coralie = coralie(os.path.join(PATH_SPECTRA, starname_coralie), abundances = abundances)

                if snr_coralie != 0.:
                    modify_file(starname_coralie, 'coralie', snr_coralie, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_coralie)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_coralie + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_coralie + '_vsini.ares'
                else:
                    if snr_coralie != 0.:
                        modify_file_vsini(starname_coralie, 'coralie', snr_coralie, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


    #############################################
    # AAT spectra
    #############################################

    if use_AAT == True:

        if os.path.isfile(spectra_aat):
            snr_aat = 0.0
            nfiles += 1
            print '\t\tFound AAT spectra'

            starname_aat = starname + '_aat'
            append_inst = 0

            if abundances == False:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_aat + '.ares'))
            else:
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_aat + '_ab.ares'))
                print '\t\tComputing with abundances = True'

            if file_skip and skip:
                print '\t\tThere is already a file called ' + starname_aat + '.ares'
                append_inst += 1

            if (file_skip == False) or (skip == False):
                snr_aat = aat(os.path.join(PATH_SPECTRA, starname_aat), abundances = abundances)

                if snr_aat != 0.:
                    modify_file(starname_aat, 'aat', snr_aat, abundances = abundances, PATH_SPECTRA = PATH_SPECTRA)

                    os.system('bash run_ARES.bash')
                    append_inst += 1

            if append_inst > 0:
                name_instrument.append(starname_aat)

            if vsini:
                print '\t\tComputing for vsini lines'
                file_skip = os.path.isfile(os.path.join(PATH_EW, starname_aat + '_vsini.ares'))
                if file_skip and skip:
                    print '\t\tThere is already a file called ' + starname_aat + '_vsini.ares'
                else:
                    if snr_aat != 0.:
                        modify_file_vsini(starname_aat, 'aat', snr_aat, PATH_SPECTRA = PATH_SPECTRA)
                        os.system('bash run_ARES.bash')


    if nfiles == 0:
        print '\n'
        print '\t\tNo spectra found for ' + starname
        print '\n'

    return name_instrument


#######################################################################
#######################################################################


def spectra(lista_estrellas, skip = True,\
            HARPS = True, FEROS = True, UVES = True, HIRES = True, \
            CORALIE = True, AAT = True, PFS = True, abundance = True, \
            s_path = './Spectra'):

    skip = True
    total_spectra = []
    calc_abundance = []
    PATH_EW = './EW'

    for starname in lista_estrellas:
        name_instrument = compute_ew(starname, skip = skip, \
                                     abundances = False, \
                                     use_HARPS = HARPS, \
                                     use_FEROS = FEROS, \
                                     use_UVES = UVES, \
                                     use_HIRES = HIRES, \
                                     use_CORALIE = CORALIE, \
                                     use_AAT = AAT, \
                                     use_PFS = PFS, \
                                     vsini = True, \
                                     PATH_SPECTRA = s_path)
        if abundance:
            n = compute_ew(starname, skip = skip, abundances = True, \
                           use_HARPS = HARPS, \
                           use_FEROS = FEROS, \
                           use_UVES = UVES, \
                           use_HIRES = HIRES, \
                           use_CORALIE = CORALIE, \
                           use_AAT = AAT, \
                           use_PFS = PFS, \
                           PATH_SPECTRA = s_path)


        name_instrument_final = []
        inst_abundance = []

        for inst in name_instrument:
            use_ab = False
            if os.path.isfile(os.path.join(PATH_EW, inst + '.ares')):
                if (os.stat(os.path.join(PATH_EW, inst + '.ares')).st_size != 0):

                    name_instrument_final.append(inst)

                    if os.path.isfile(os.path.join(PATH_EW, inst + '_ab.ares')):
                        if (os.stat(os.path.join(PATH_EW, inst + '_ab.ares')).st_size != 0):
                            use_ab = True

                    inst_abundance.append(use_ab)

        total_spectra = total_spectra + name_instrument_final
        calc_abundance = calc_abundance + inst_abundance

    return total_spectra, calc_abundance
