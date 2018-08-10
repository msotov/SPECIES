'''
Computes the stellar parameters for the EWs as given by ARES.

Please send your questions and comments to
maritsoto@ug.uchile.cl

Created in: 31/07/2015
Last modified in: 14/03/2018
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
                           'to a certain value. \
                           DEPRECATED.')
parser.add_argument('-file_mass', default = 'any_file.txt',\
                    help = 'Name of the file with the hold values for the mass, '\
                           'age and photometric logg. '\
                           'Will only be used if hold_mass == True. \
                           DEPRECATED.')
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

parser.add_argument('-vt_hold', action = 'store_true', default = False,\
                    help = 'Try the computation with vt = 1.2 km/s in the case it \
                            does not converge. Default is False.')

parser.add_argument('-no_r_iso_logg', action = 'store_true', default = False,\
                    help = 'Do not use the logg from isochrones interpolation \
                            to recompute the atmospheric parameters when there is \
                            disagreement with the spectroscopic one.')
parser.add_argument('-age_limits',\
                    help = 'Set limits to the parameters space for the stellar age. \
                            It has to be in the format log(age), with age in years \
                            (not Gyr). Default is None.')
parser.add_argument('-path_spectra', default = './Spectra',\
                    help = 'Path to the directory containing the spectra. It has \
                            to be relative to the path of SPECIES, or the \
                            absolute path.')

parser.add_argument('-columns', default = 'All',\
                    help = 'Columns to include in the output. The name of the star, \
                            instrument and exception will always be included.')

parser.add_argument('-giant', default = 'False', action = 'store_true',\
                    help = 'If the star, or list of stars, are considered to be giants.')


args = parser.parse_args()

# -output summary_N021969 -use_coords -starlist N021969

#########################################################################
# List of stars for which to compute the parameters
#########################################################################

import numpy as np
import glob
from astropy.io import fits, ascii
from astropy.table import vstack

starlist = args.starlist.split(',')
if 'any_star' in starlist:
    print '\t Starlist was not given in the command line. ' + \
          'Reading from inside SPECIES.py\n'


    #starlist = ['ceres01','ceres02','ceres03','moon01','ganymede01','sun01','sun02','sun03','sun04','sun05']


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
# Check with  temperature with the one from photometric relations
#########################################################################

colors_from_file = False
check_T = not args.no_check_t
name_file_with_colors = args.file_colors
if name_file_with_colors != None:
    colors_from_file = True
use_coords = args.use_coords
file_with_coords = args.file_coords

#########################################################################
# Try with vt = 1.2 km/s if the computation does not converge
#########################################################################

try_with_vt_hold = args.vt_hold

#########################################################################
# Recompute atm parameters with logg_p
#########################################################################

recompute_with_logg_p = not args.no_r_iso_logg


#########################################################################
# Age limits
#########################################################################

age_limits = args.age_limits
if age_limits != None:
    a = age_limits.split(',')
    age_limits = (float(a[0]), float(a[1]))

#########################################################################
# Path to spectra
#########################################################################

PATH_SPECTRA = args.path_spectra

#########################################################################
# Columns to use in the output
#########################################################################

out_columns = args.columns

#########################################################################
# If the star is giant
#########################################################################

is_giant = args.giant

#########################################################################
# Imports the python packages needed
# Make sure you've got all of them
#########################################################################

print '\t Loading the necessary packages...'

import os, time, warnings, logging

warnings.simplefilter("ignore")

import calc_params as cp
import multiprocessing as mp
from spectra import spectra
from calc_params import multi_run_wrapper as multi_run_wrapper
from astropy.table import Table, Column


#########################################################################
# Initializes the time of execution
#########################################################################

start = time.time()

#########################################################################
# Corrects the spectra to restframe and computes the EW
#########################################################################

print '\t Creating the equivalent width files...'

if all_inst == True:
    starname, abundance_p = spectra(starlist, HARPS = True, FEROS = True, UVES = True, \
                       HIRES = True, CORALIE = True, AAT = True, PFS = True, \
                       abundance = abundance, s_path = PATH_SPECTRA)
else:
    starname, abundance_p = spectra(starlist, HARPS = harps, FEROS = feros, UVES = uves, \
                       HIRES = hires, CORALIE = coralie, AAT = aat, PFS = pfs, \
                       abundance = abundance, s_path = PATH_SPECTRA)

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
                        set_boundaries, file_with_boundaries, mass_from_file,\
                        try_with_vt_hold, recompute_with_logg_p, age_limits,\
                        PATH_SPECTRA, is_giant))


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
        err_vt2, vs, err_vs, vm, err_vm, mass, err_mass1, err_mass2, age,\
        err_age1, err_age2, logg_iso, err_logg_iso1, err_logg_iso2, radius,\
        err_radius1, err_radius2, logL, err_logL1, err_logL2, use_Tc, use_vt, use_logg_p,\
        T_c, err_T_c, relation_temp = results.transpose()



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
    # Sets the instrument and rv columns
    ############################################################################

    name_star = star[:]
    instrument = np.array(['nofound' for i in range(len(star))])
    rv = np.zeros(len(star))
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

        h = fits.getheader('Spectra/' + name_star[i] + '_' + name_inst + '_res.fits', 0)
        try:
            rv[i] = h['RV']
        except KeyError:
            pass
        del h

        if os.path.samefile(PATH_SPECTRA, './Spectra') == False:
            if os.path.isfile(os.path.join('./Spectra', name_star[i] + '_' + name_inst + '_res.fits')):
                os.unlink(os.path.join('./Spectra', name_star[i] + '_' + name_inst + '_res.fits'))


    ############################################################################
    # Writes the results in the output files
    ############################################################################

    original_columns = np.array(['Starname', 'Instrument', 'RV', '[Fe/H]', 'err_[Fe/H]', \
    'Temperature', 'err_T', 'logg', 'err_logg', 'vt', 'err_vt', 'nFeI', 'nFeII',\
    'exception', 'vsini', 'err_vsini', 'vmac', 'err_vmac', '[Na/H]', 'e_[Na/H]',\
    'nNaI', '[Mg/H]', 'e_[Mg/H]', 'nMgI', '[Al/H]', 'e_[Al/H]', 'nAlI', '[Si/H]',\
    'e_[Si/H]', 'nSiI', '[Ca/H]', 'e_[Ca/H]', 'nCaI', '[TiI/H]', 'e_[TiI/H]', 'nTiI',\
    '[TiII/H]', 'e_[TiII/H]', 'nTiII', '[Cr/H]', 'e_[Cr/H]', 'nCrI', '[Mn/H]', \
    'e_[Mn/H]', 'nMnI','[Ni/H]', 'e_[Ni/H]', 'nNiI', '[Cu/H]', 'e_[Cu/H]', 'nCuI', \
    '[Zn/H]', 'e_[Zn/H]', 'nZnI', 'exception_Fe', 'exception_Ti', '[FeI/H]', '[FeII/H]',\
    'Mass', 'err_mass1', 'err_mass2', 'Age', 'err_age1', 'err_age2', 'iso_logg', \
    'err_iso_logg1', 'err_iso_logg2', 'Radius', 'err_radius1', 'err_radius2', \
    'logL', 'err_logL1', 'err_logL2', 'use_Tc', 'use_vt', 'use_iso_logg', 'err_vt2', 'err_T2',\
    'T_photo', 'err_T_photo', 'T_photo_relation'])

    data = Table([name_star, instrument, rv, xmetal, err_met, T, err_T,\
            logg, err_logg, micro, err_vt, nFeI, nFeII, exception, vs, err_vs,\
            vm, err_vm, ab_NaI, dev_NaI, nNaI, ab_MgI, dev_MgI, nMgI,\
            ab_AlI, dev_AlI, nAlI, ab_SiI, dev_SiI, nSiI, ab_CaI, dev_CaI, nCaI,\
            ab_TiI, dev_TiI, nTiI, ab_TiII, dev_TiII, nTiII, ab_CrI, dev_CrI, nCrI,\
            ab_MnI, dev_MnI, nMnI, ab_NiI, dev_NiI, nNiI, ab_CuI, dev_CuI, nCuI, \
            ab_ZnI, dev_ZnI, nZnI, exception_Fe, exception_Ti, ab_FeI, ab_FeII, \
            mass, err_mass1, err_mass2, age, err_age1, err_age2, logg_iso, err_logg_iso1,\
            err_logg_iso2, radius, err_radius1, err_radius2, logL, err_logL1, err_logL2,\
            use_Tc, use_vt, use_logg_p, err_vt2, err_T2, T_c, err_T_c, relation_temp],\
            names = ['Starname', 'Instrument', 'RV', '[Fe/H]', 'err_[Fe/H]', \
            'Temperature', 'err_T', 'logg', 'err_logg', 'vt', 'err_vt', 'nFeI', 'nFeII',\
            'exception', 'vsini', 'err_vsini', 'vmac', 'err_vmac', '[Na/H]', 'e_[Na/H]',\
            'nNaI', '[Mg/H]', 'e_[Mg/H]', 'nMgI', '[Al/H]', 'e_[Al/H]', 'nAlI', '[Si/H]',\
            'e_[Si/H]', 'nSiI', '[Ca/H]', 'e_[Ca/H]', 'nCaI', '[TiI/H]', 'e_[TiI/H]', 'nTiI',\
            '[TiII/H]', 'e_[TiII/H]', 'nTiII', '[Cr/H]', 'e_[Cr/H]', 'nCrI', '[Mn/H]', \
            'e_[Mn/H]', 'nMnI','[Ni/H]', 'e_[Ni/H]', 'nNiI', '[Cu/H]', 'e_[Cu/H]', 'nCuI', \
            '[Zn/H]', 'e_[Zn/H]', 'nZnI', 'exception_Fe', 'exception_Ti', '[FeI/H]', '[FeII/H]',\
            'Mass', 'err_mass1', 'err_mass2', 'Age', 'err_age1', 'err_age2', 'iso_logg', \
            'err_iso_logg1', 'err_iso_logg2', 'Radius', 'err_radius1', 'err_radius2', \
            'logL', 'err_logL1', 'err_logL2', 'use_Tc', 'use_vt', 'use_iso_logg', 'err_vt2', 'err_T2',\
            'T_photo', 'err_T_photo', 'T_photo_relation'],\
            dtype = ['str', 'str', 'float', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', 'float', 'int', \
            'int', 'int', 'float', 'float', 'float', 'float', \
            'float', 'float', 'int', 'float', 'float', 'int', 'float', 'float', \
            'int', 'float', 'float', 'int', 'float', 'float', 'int', 'float', \
            'float', 'int', 'float', 'float', 'int', 'float', 'float', 'int', 'float', 'float', \
            'int', 'float', 'float', 'int', 'float', 'float', 'int', 'float', \
            'float', 'int', 'int', 'int', 'float', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',\
            'float', 'float', 'float', 'str', 'str',\
            'str', 'float', 'float','float', 'float', 'str'])

    data['[Fe/H]'].format = data['[Na/H]'].format = data['[Mg/H]'].format = \
                data['[Al/H]'].format = data['[Si/H]'].format = \
                data['[Ca/H]'].format = data['[TiI/H]'].format = '%8.2f'
    data['[Cr/H]'].format = data['[Mn/H]'].format = data['[Ni/H]'].format = \
                data['[Cu/H]'].format = data['[Zn/H]'].format = data['[TiII/H]'].format = '%8.2f'
    data['Temperature'].format = '.1f'
    data['Mass'].format = data['Age'].format = data['Radius'].format = data['logL'].format = '%8.2f'
    data['logg'].format = data['vt'].format = \
                data['vsini'].format = data['err_vsini'].format = \
                data['vmac'].format = data['err_vmac'].format = \
                data['e_[Na/H]'].format = data['e_[Mg/H]'].format = \
                data['e_[Al/H]'].format = data['e_[Si/H]'].format = \
                data['e_[Ca/H]'].format = data['e_[TiI/H]'].format = data['e_[TiII/H]'].format = \
                data['e_[Cr/H]'].format = data['e_[Mn/H]'].format = \
                data['e_[Ni/H]'].format = data['e_[Cu/H]'].format = \
                data['e_[Zn/H]'].format = data['[FeI/H]'].format = \
                data['[FeII/H]'].format = data['err_T'].format = \
                data['err_T2'].format = data['err_logg'].format = \
                data['err_[Fe/H]'].format = data['err_vt'].format = \
                data['err_vt2'].format = '%.3f'
    data['err_mass1'].format = data['err_mass2'].format = \
                data['err_age1'].format = data['err_age2'].format = data['iso_logg'].format = \
                data['err_iso_logg1'].format = data['err_iso_logg2'].format = \
                data['err_radius1'].format = data['err_radius2'].format = \
                data['err_logL1'].format = data['err_logL2'].format = '%.3f'
    data['RV'].format = '%.3f'
    data['T_photo'].format = '%.1f'
    data['err_T_photo'].format = '%.3f'

    if out_columns != 'All':
        columns = np.array(['Starname', 'Instrument', 'exception'] + out_columns.split(','))
        not_included = original_columns[np.isin(original_columns, columns, invert = True)]
        data.remove_columns(not_included)


    if save_as_ascii == True:

        ascii.write(data, './output/' + file_output + '.dat', \
                    format = 'fixed_width', delimiter = None, overwrite = True)


    if save_as_fits == True:

        data.write('./output/' + file_output + '.fits', format='fits', overwrite = True)




time_seconds = (time.time() - start)

print '\n\n\t Time in seconds: %f' % time_seconds
print '\t Time in minutes: %f' % (time_seconds/60.)
print '\t Time in hours:   %f' % (time_seconds/3600.)

print '\n\t*****************************************\n'
print '\t Thank you for using SPECIES\n'
print '\t*****************************************\n'
