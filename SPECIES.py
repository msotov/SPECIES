'''
Computes the stellar parameters for the EWs as given by ARES.

Please send your questions and comments to
maritsoto@ug.uchile.cl

Created in: 31/07/2015
Last modified in: 31/05/2017
'''

print '\n\n\t*****************************************'
print '\t WELCOME TO SPECIES\n'
print '\t Written by Maritza Soto'
print '\t Send your comments to maritsoto@ug.uchile.cl'
print '\t*****************************************\n'
print '\t Reading the input values...'


#########################################################################
# Read the input parameters from command line
#########################################################################

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-starlist', default = 'any_star',\
                    help = 'Name of the stars to use, separated by commas. '\
                           'If the list is too long, retrieve it inside '\
                           'the code, in the part that says '\
                           '# List of stars for which to compute the parameters.')
parser.add_argument('-inst', default = 'all',\
                    help = 'Names of the instruments you want to use '\
                           'separated by commas. If "all" is included, '\
                           'spectra from all accepted instruments will be used')
parser.add_argument('-ncores', default = 3, type = int,\
                    help = 'Number of cores to use.')
parser.add_argument('-output_format', default = 'Both',\
                    choices = ['Both', 'ASCII', 'FITS'],\
                    help = 'Format of the output file.')
parser.add_argument('-output', default = 'summary_output',\
                    help = 'Name of the output file, without extension.')
parser.add_argument('-hold_params', action = 'store_true', default = False,\
                    help = 'Option to set the atmospheric parameters to a '\
                           'certain value.')
parser.add_argument('-file_params', default = 'any_file.txt',\
                    help = 'Name of the file with the parameters to hold. ' \
                           'Will be used only if hold_params == True.')
parser.add_argument('-hold_mass', action = 'store_true', default = False,\
                    help = 'Option to set the mass, age, and photometric logg '\
                           'to a certain value.')
parser.add_argument('-file_mass', default = 'any_file.txt',\
                    help = 'Name of the file with the hold values for the mass, '\
                           'age and photometric logg. '\
                           'Will only be used if hold_mass == True.')
parser.add_argument('-set_boundaries', action = 'store_true', default = False,\
                    help = 'Option to change the boundaries of the '\
                           'of the atmospheric parameters.')
parser.add_argument('-file_with_boundaries', default = 'any_file.txt',\
                    help = 'Name of the file with the boundaries for the '\
                           'atmospheric parameters. '\
                           'Will only be used if set_boundaries == True.')
parser.add_argument('-no_abundance', action = 'store_true', default = False, \
                    help = 'Option to compute abundances for a set of elements.')
parser.add_argument('-no_error', action = 'store_true', default = False,\
                    help = 'Option to compute uncertainties to the '\
                           'atmospheric parameters.')
parser.add_argument('-no_debug', action = 'store_true', default = False,\
                    help = 'Option to create a debug file for each star, '\
                           'where the details of the computation will be stored.')
parser.add_argument('-no_check_t', action = 'store_true', default = False,\
                    help = 'Check the temperature with the one computed using '\
                           'the formulas from Casagrande et al. 2012.')
parser.add_argument('-file_colors',\
                    help = 'Name of the file with the colors to use ' \
                           'to compute the temperature from Casagrande et al. '\
                           'Used only if check_t == True')
parser.add_argument('-use_coords', action = 'store_true', default = False,\
                    help = 'Option to use the coordinates of the star, '\
                           'instead of its name, to look for the star colors. '\
                           'Used only if check_t == True')
parser.add_argument('-file_coords', default = 'any_file.txt',\
                    help = 'Name of the file wiith the coordinates '\
                           'of the stars (in deg). '\
                           'Used only if check_t == True and use_coords == True.')

parser.add_argument('-mass_from_file', action = 'store_true', default = False,\
                    help = 'Load the chains created during a previous mass calculation,\
                            stored in dotter_isochrones/*.h5')


args = parser.parse_args()

# -output summary_N021969 -use_coords -starlist N021969

#########################################################################
# List of stars for which to compute the parameters
#########################################################################

import numpy as np
import glob

starlist = args.starlist.split(',')
if 'any_star' in starlist:
    print '\t Starlist was not given in the command line. ' + \
          'Reading from inside SPECIES.py\n'


    #starlist = ['sun01']
    #starlist = ['ceres01','ceres02','ceres03','moon01','ganymede01','sun01','sun02','sun03','sun04','sun05']
    #starlist = ['ceres03']
    #starlist = ['HD10002']

    #starlist = ['GL397', 'HD9986', 'HD72490', 'HD111814', 'HD171918', 'HD202696', 'HD214823', 'HD26965']
    #starlist = ['GL397']

    #starlist = ['HD100942', 'HD11112']
    #starlist = ['HD26965']
    starlist = np.genfromtxt('./Spectra/stars_sousa.txt', dtype = None, unpack = True, delimiter = '\n')
    starlist = [str(s).replace(' ', '') for s in starlist]
    #starlist = ['HD11964A']

    #starlist = ['arcturus01', 'arcturus02', 'arcturus03']
    #starlist = ['arcturus03']
    #print starlist
    #press = raw_input('press')
    #starlist = ['HD28185']

    #starlist = ['NG021676']
    #starlist = ['N021969']

    #files = glob.glob('./mari_cheps/*.fits')
    #names = [f[30:-5] for f in files]
    #for i in range(len(names)):
    #    if names[i][:2] == 'HD' or names[i][:3] == 'HIP':
    #        names[i] = names[i].replace('-','')
    #starlist = np.hstack((starlist,names))

    #from astropy.io import ascii

    #data1 = ascii.read('./output/summary_stars_sousa_all_inst_test3.dat')
    #data2 = ascii.read('./output/summary_cheps_stars.dat')
    #data = np.hstack((data1, data2))

    #exception = data['exception']
    #use_Tc = data['use_Tc']
    #inst = data['Instrument']

    #i = np.where((exception == 1) & (inst == 'HARPS'))[0]

    #starlist = data['Starname'][i]

'''
    starlist = ['HD283', 'HD6348', 'HD9796', 'HD10647', 'HD12345', 'HD12617', 'HD14744',
 'HD17051', 'HD34688', 'HD38858', 'HD44573', 'HD65907A', 'HD69830', 'HD72673',
 'HD74014', 'HD83529', 'HD85512', 'HD86065', 'HD93380', 'HD98356', 'HD100508',
 'HD100777', 'HD101581', 'HD104263', 'HD114386', 'HD116858', 'HD116920',
 'HD123265', 'HD129642', 'HD134606', 'HD143295', 'HD146233', 'HD154363',
 'HD176157', 'HD177565', 'HD199933', 'HD207129', 'HD210277', 'HD210975',
 'HD212301', 'HD216777', 'HD221420', 'Gamma-Tau', 'HD126535', 'HD38467',
 'HD68402', 'HD72892', 'HD77338', 'HD9174', 'HIP4909', 'HIP6407', 'HIP8507',
 'HIP11915', 'HIP30037', 'HIP31831', 'HIP53084', 'HIP69724']
'''
#########################################################################
# Instruments to use
#########################################################################

all_inst, harps, feros, uves, hires, coralie, aat, pfs = [False, False, False, False,\
                                                     False, False, False, False]
inst = args.inst.split(',')
if 'all' in inst: all_inst = True
if 'harps' in inst: harps = True
if 'feros' in inst: feros = True
if 'uves' in inst: uves = True
if 'hires' in inst: hires = True
if 'coralie' in inst: coralie = True
if 'aat' in inst: aat = True
if 'pfs' in inst: pfs = True



#########################################################################
# Number of cores
#########################################################################

ncores = args.ncores

#########################################################################
# Output options
#########################################################################

output_format = args.output_format
file_output = args.output

#########################################################################
# hold atmospheric parameters option
#########################################################################

hold_params = args.hold_params
file_params = args.file_params

#########################################################################
# hold mass, age, and photometric logg option
#########################################################################

hold_mass = args.hold_mass
file_hold_mass = args.file_mass
mass_from_file = args.mass_from_file

#########################################################################
# Set boundaries for the atmospheric parameters
#########################################################################

set_boundaries = args.set_boundaries
file_with_boundaries = args.file_with_boundaries

#########################################################################
# Abundance option
#########################################################################

abundance = not args.no_abundance

#########################################################################
# Uncertainties option
#########################################################################

error = not args.no_error
if error:
    file_coefs_error = True
    plots_coefs_error = True

#########################################################################
# Debugging option
#########################################################################

set_debug = not args.no_debug

#########################################################################
# Check with the Casagrande et al. 2010 values for the temperature
#########################################################################

colors_from_file = False
check_T = not args.no_check_t
name_file_with_colors = args.file_colors
if name_file_with_colors != None:
    colors_from_file = True
use_coords = args.use_coords
file_with_coords = args.file_coords



#########################################################################
# Imports the python packages needed
# Make sure you've got all of them
#########################################################################

print '\t Loading the necessary packages...'

import os, time, warnings, logging
import calc_params as cp
import multiprocessing as mp
from spectra import spectra
from calc_params import multi_run_wrapper as multi_run_wrapper
from astropy.io import fits, ascii
from astropy.table import Table, Column


#########################################################################
# Initializes the time of execution
#########################################################################

start = time.time()

#########################################################################
# Ignore warnings
#########################################################################

warnings.simplefilter("ignore")

#########################################################################
# Corrects the spectra to restframe and computes the EW
#########################################################################

print '\t Creating the equivalent width files...'

if all_inst == True:
    starname, abundance_p = spectra(starlist, HARPS = True, FEROS = True, UVES = True, \
                       HIRES = True, CORALIE = True, AAT = True, PFS = True, \
                       abundance = abundance)
else:
    starname, abundance_p = spectra(starlist, HARPS = harps, FEROS = feros, UVES = uves, \
                       HIRES = hires, CORALIE = coralie, AAT = aat, PFS = pfs, \
                       abundance = abundance)

if len(starname) == 0:
    print '\t No star with valid .ares files was found. Stopping the computation.'

else:

    ############################################################################
    # Check that the list of stars is correct
    ############################################################################

    print '\n\n\t The computation will be carried out for the following stars: '

    starname = starname
    list_names = ''
    for s in starname:
        list_names += s + ', '
    print '\t' + list_names[:-2]

    #press = raw_input('\t Press enter to continue, or ctrl+C to abort.')


    ############################################################################
    # Creates the list of arguments to be used in the parallelization
    ############################################################################


    arg_list = []
    for i in range(len(starname)):
        arg_list.append((starname[i], set_debug, hold_params, abundance_p[i],\
                        error, file_params, file_coefs_error, plots_coefs_error,\
                        check_T, colors_from_file, name_file_with_colors,\
                        hold_mass, file_hold_mass,\
                        use_coords, file_with_coords,\
                        set_boundaries, file_with_boundaries, mass_from_file))


    ############################################################################
    # Run the computation of the parameters, with the stars split into
    # different cores.
    ############################################################################


    pool = mp.Pool(ncores, maxtasksperchild = 20)

    results = pool.map(multi_run_wrapper, arg_list, chunksize = 1)

    results = np.array(results)

    pool.close()


    i_none = np.where(results.T[0] == None)[0]
    results = np.delete(results, i_none, 0)


    del i_none


    ############################################################################
    # Writes the output in each variable
    ############################################################################


    star, T, logg, xmetal, micro, exception, nFeI, nFeII, ab_FeI, ab_FeII,\
        ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI, ab_AlI, dev_AlI, nAlI,\
        ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI, ab_TiI, dev_TiI, nTiI,\
        ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI, ab_MnI, dev_MnI, nMnI,\
        ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI, nCuI, ab_ZnI, dev_ZnI, nZnI,\
        exception_Fe, exception_Ti, err_T, err_T2, err_logg, err_met, err_vt,\
        err_vt2, vs, err_vs, vm, err_vm, mass, err_mass,\
        age, err_age, s_logg, err_s_logg, use_Tc = results.transpose()



    ############################################################################
    # Options for output file
    ############################################################################


    if output_format == 'FITS':
        save_as_ascii = False
        save_as_fits = True
    elif output_format == 'ASCII':
        save_as_ascii = True
        save_as_fits = False
    else:
        save_as_ascii = True
        save_as_fits = True


    ############################################################################
    # Sets the instrument column
    ############################################################################

    name_star = star
    instrument = np.array(['nofound' for i in range(len(star))])
    for i in range(len(star)):
        indices = [indx for indx, c in enumerate(star[i]) if '_' == c]
        if len(indices) > 0:
            name_inst = star[i][indices[-1]+1:]
            if name_inst == 'o':
                name_inst = star[i][indices[-2]+1:]
                name_star[i] = star[i][:indices[-2]]
            else:
                name_star[i] = star[i][:indices[-1]]
            if (name_inst == 'harps'):
                instrument[i] = 'HARPS'
            elif (name_inst == 'feros'):
                instrument[i] = 'FEROS'
            elif (name_inst == 'feros_o'):
                instrument[i] = 'FEROS'
            elif (name_inst == 'uves'):
                instrument[i] = 'UVES'
            elif (name_inst == 'hires') or (name_inst == 'HIRES'):
                instrument[i] = 'HIRES'
            elif (name_inst == 'coralie'):
                instrument[i] = 'CORALIE'
            elif (name_inst == 'aat'):
                instrument[i] = 'AAT'
            elif (name_inst == 'pfs'):
                instrument[i] = 'PFS'


    ############################################################################
    # Writes the results in the output files
    ############################################################################

    if save_as_ascii == True:

        data = Table([name_star, instrument, xmetal, T, logg, micro, \
            nFeI, nFeII, exception, vs, err_vs, vm, \
            err_vm, ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI,\
            ab_AlI, dev_AlI, nAlI, ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, \
            nCaI, ab_TiI, dev_TiI, nTiI, nTiII, ab_CrI, dev_CrI, nCrI, \
            ab_MnI, dev_MnI, nMnI, ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI,\
            nCuI, ab_ZnI, dev_ZnI, nZnI, exception_Fe, exception_Ti, \
            ab_FeI, ab_FeII, err_T, err_T2, err_logg, err_met, err_vt, err_vt2, \
            mass, err_mass, age, err_age, s_logg, err_s_logg, use_Tc],\
            names = ['Starname', 'Instrument', '[Fe/H]', 'Temperature', 'logg', \
            'vt', 'nFeI', 'nFeII', 'exception',\
            'vsini', 'err_vsini', 'vmac', 'err_vmac',\
            '[Na/H]', 'e_[Na/H]', 'nNaI', '[Mg/H]', 'e_[Mg/H]', 'nMgI', \
            '[Al/H]', 'e_[Al/H]', 'nAlI', '[Si/H]', 'e_[Si/H]', 'nSiI', \
            '[Ca/H]', 'e_[Ca/H]', 'nCaI', '[Ti/H]', 'e_[Ti/H]', 'nTiI', 'nTiII',\
            '[Cr/H]', 'e_[Cr/H]', 'nCrI', '[Mn/H]', 'e_[Mn/H]', 'nMnI',\
            '[Ni/H]', 'e_[Ni/H]', 'nNiI', '[Cu/H]', 'e_[Cu/H]', 'nCuI', \
            '[Zn/H]', 'e_[Zn/H]', 'nZnI', 'exception_Fe', 'exception_Ti', \
            '[FeI/H]', '[FeII/H]', 'err_T', 'err_T2', 'err_logg', 'err_[Fe/H]', \
            'err_vt', 'err_vt2', 'Mass', 'err_mass', 'age', 'err_age', \
            'Photo_logg', 'err_photo_logg', 'use_Tc'],\
            dtype = ['str', 'str', 'float', 'float', 'float', 'float', 'int', \
            'int', 'int', 'float', 'float', 'float', 'float', \
            'float', 'float', 'int', 'float', 'float', 'int', 'float', 'float', \
            'int', 'float', 'float', 'int', 'float', 'float', 'int', 'float', \
            'float', 'int', 'int', 'float', 'float', 'int', 'float', 'float', \
            'int', 'float', 'float', 'int', 'float', 'float', 'int', 'float', \
            'float', 'int', 'int', 'int', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', 'str'])
        data['[Fe/H]'].format = data['[Na/H]'].format = data['[Mg/H]'].format = \
                data['[Al/H]'].format = data['[Si/H]'].format = \
                data['[Ca/H]'].format = data['[Ti/H]'].format = '%8.2f'
        data['[Cr/H]'].format = data['[Mn/H]'].format = data['[Ni/H]'].format = \
                data['[Cu/H]'].format = data['[Zn/H]'].format = '%8.2f'
        data['Temperature'].format = '.1f'
        data['Mass'].format = data['age'].format = '%8.2f'
        data['logg'].format = data['vt'].format = \
                data['vsini'].format = data['err_vsini'].format = \
                data['vmac'].format = data['err_vmac'].format = \
                data['e_[Na/H]'].format = data['e_[Mg/H]'].format = \
                data['e_[Al/H]'].format = data['e_[Si/H]'].format = \
                data['e_[Ca/H]'].format = data['e_[Ti/H]'].format = \
                data['e_[Cr/H]'].format = data['e_[Mn/H]'].format = \
                data['e_[Ni/H]'].format = data['e_[Cu/H]'].format = \
                data['e_[Zn/H]'].format = data['[FeI/H]'].format = \
                data['[FeII/H]'].format = data['err_T'].format = \
                data['err_T2'].format = data['err_logg'].format = \
                data['err_[Fe/H]'].format = data['err_vt'].format = \
                data['err_vt2'].format = data['err_mass'].format = \
                data['err_age'].format = data['Photo_logg'].format = \
                data['err_photo_logg'].format = '.3f'

        ascii.write(data, './output/' + file_output + '.dat', \
                    format = 'fixed_width', delimiter = None)


    if save_as_fits == True:
        c1 = fits.Column(name = 'Starname', format = '20A', array = star)
        c2 = fits.Column(name = '[Fe/H]', format = 'E', array = xmetal)
        c3 = fits.Column(name = 'Temperature', format = 'E', array = T)
        c4 = fits.Column(name = 'logg', format = 'E', array = logg)
        c5 = fits.Column(name = 'vt', format = 'E', array = micro)
        c6 = fits.Column(name = 'nFeI', format = 'J', array = nFeI)
        c7 = fits.Column(name = 'nFeII', format = 'J', array = nFeII)
        c8 = fits.Column(name = 'exception', format = 'J', array = exception)
        c9 = fits.Column(name = '[FeI/H]', format = 'E', array = ab_FeI)
        c10 = fits.Column(name = '[FeII/H]', format = 'E', array = ab_FeII)

        c13 = fits.Column(name = 'err_T', format = 'E', array = err_T)
        c14 = fits.Column(name = 'err_T2', format = 'E', array = err_T2)
        c15 = fits.Column(name = 'err_logg', format = 'E', array = err_logg)
        c16 = fits.Column(name = 'err_[Fe/H]', format = 'E', array = err_met)
        c17 = fits.Column(name = 'err_vt', format = 'E', array = err_vt)
        c66 = fits.Column(name = 'err_vt2', format = 'E', array = err_vt2)

        c23 = fits.Column(name = '[Na/H]', format = 'E', array = ab_NaI)
        c24 = fits.Column(name = 'e_[Na/H]', format = 'E', array = dev_NaI)
        c25 = fits.Column(name = 'nNaI', format = 'J', array = nNaI)
        c26 = fits.Column(name = '[Mg/H]', format = 'E', array = ab_MgI)
        c27 = fits.Column(name = 'e_[Mg/H]', format = 'E', array = dev_MgI)
        c28 = fits.Column(name = 'nMgI', format = 'J', array = nMgI)
        c29 = fits.Column(name = '[Al/H]', format = 'E', array = ab_AlI)
        c30 = fits.Column(name = 'e_[Al/H]', format = 'E', array = dev_AlI)
        c31 = fits.Column(name = 'nAlI', format = 'J', array = nAlI)
        c32 = fits.Column(name = '[Si/H]', format = 'E', array = ab_SiI)
        c33 = fits.Column(name = 'e_[Si/H]', format = 'E', array = dev_SiI)
        c34 = fits.Column(name = 'nSiI', format = 'J', array = nSiI)
        c35 = fits.Column(name = '[Ca/H]', format = 'E', array = ab_CaI)
        c36 = fits.Column(name = 'e_[Ca/H]', format = 'E', array = dev_CaI)
        c37 = fits.Column(name = 'nCaI', format = 'J', array = nCaI)
        c38 = fits.Column(name = '[Ti/H]', format = 'E', array = ab_TiI)
        c39 = fits.Column(name = 'e_[Ti/H]', format = 'E', array = dev_TiI)
        c40 = fits.Column(name = 'nTiI', format = 'J', array = nTiI)
        c41 = fits.Column(name = 'nTiII', format = 'J', array = nTiII)
        c42 = fits.Column(name = '[Cr/H]', format = 'E', array = ab_CrI)
        c43 = fits.Column(name = 'e_[Cr/H]', format = 'E', array = dev_CrI)
        c44 = fits.Column(name = 'nCrI', format = 'J', array = nCrI)
        c45 = fits.Column(name = '[Mn/H]', format = 'E', array = ab_MnI)
        c46 = fits.Column(name = 'e_[Mn/H]', format = 'E', array = dev_MnI)
        c47 = fits.Column(name = 'nMnI', format = 'J', array = nMnI)
        c48 = fits.Column(name = '[Ni/H]', format = 'E', array = ab_NiI)
        c49 = fits.Column(name = 'e_[Ni/H]', format = 'E', array = dev_NiI)
        c50 = fits.Column(name = 'nNiI', format = 'J', array = nNiI)
        c51 = fits.Column(name = '[Cu/H]', format = 'E', array = ab_CuI)
        c52 = fits.Column(name = 'e_[Cu/H]', format = 'E', array = dev_CuI)
        c53 = fits.Column(name = 'nCuI', format = 'J', array = nCuI)
        c54 = fits.Column(name = '[Zn/H]', format = 'E', array = ab_ZnI)
        c55 = fits.Column(name = 'e_[Zn/H]', format = 'E', array = dev_ZnI)
        c56 = fits.Column(name = 'nZnI', format = 'J', array = nZnI)
        c57 = fits.Column(name = 'exception_Fe', format = 'J', array = exception_Fe)
        c58 = fits.Column(name = 'exception_Ti', format = 'J', array = exception_Ti)

        c59 = fits.Column(name = 'Mass', format = 'E', array = mass)
        c60 = fits.Column(name = 'err_mass', format = 'E', array = err_mass)
        c61 = fits.Column(name = 'Age', format = 'E', array = age)
        c62 = fits.Column(name = 'err_age', format = 'E', array = err_age)
        c63 = fits.Column(name = 'Photo_logg', format = 'E', array = s_logg)
        c64 = fits.Column(name = 'err_photo_logg', format = 'E', array = err_s_logg)
        #c65 = fits.Column(name = 'n_mass', format = 'J', array = n_mass)

        #c67 = fits.Column(name = 'SNR', format = 'E', array = sn)
        c68 = fits.Column(name = 'use_Tc', format = '5A', array = use_Tc)

        c69 = fits.Column(name = 'vsini', format = 'E', array = vs)
        c70 = fits.Column(name = 'err_vsini', format = 'E', array = err_vs)
        c71 = fits.Column(name = 'vmac', format = 'E', array = vm)
        c72 = fits.Column(name = 'err_vmac', format = 'E', array = err_vm)

        c73 = fits.Column(name = 'Instrument', format = '20A', array = instrument)

        coldefs = fits.ColDefs([c1, c73, c2, c3, c4, c5, c6, c7, c8,\
                    c69, c70, c71, c72, c23, c24, c25, c26, c27, c28, c29, c30, c31,\
                    c32, c33, c34, c35, c36, c37, c38, c39, c40, c41, c42, c43, c44,\
                    c45, c46, c47, c48, c49, c50, c51, c52, c53, c54, c55, c56, c57,\
                    c58, c9, c10, c13, c14, c15, c16, c17, c66, c59, c60, c61, c62,\
                    c63, c64, c68])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto('./output/' + file_output + '.fits', clobber = True)


time_seconds = (time.time() - start)

print '\n\n\t Time in seconds: %f' % time_seconds
print '\t Time in minutes: %f' % (time_seconds/60.)
print '\t Time in hours:   %f' % (time_seconds/3600.)

print '\n\t*****************************************\n'
print '\t Thank you for using SPECIES\n'
print '\t*****************************************\n'
