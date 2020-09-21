"""
Computes the stellar parameters from EWs of iron lines.

Created in: 31/07/2015
Last modified in: 21/09/2020

Example of use:

python SPECIES.py -output N021969 -use_coords -starlist N021969
"""

from __future__ import print_function

__author__ = "Maritza Soto (marisotov@gmail.com)"
__date__ = "2020-09-21"
__version__ = "3.0.0"

from builtins import range
import argparse
import os
import sys
import logging
import time
import warnings
warnings.simplefilter("ignore")
import multiprocessing as mp
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table, hstack, Column
from Spectra import spectra
from CalcParams import multi_run_wrapper


print('\n\n\t*****************************************')
print('\t WELCOME TO SPECIES\n')
print('\t Written by Maritza Soto')
print('\t Send your comments to marisotov@gmail.com')
print('\t*****************************************\n')
print('\t Reading the input values...')


#########################################################################
# Read the input parameters from command line
#########################################################################

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-starlist', default='any_star',\
                    help='Name of the stars to use, separated by commas. '\
                          'If the list is too long, retrieve it inside '\
                          'the code, in the part that says '\
                          '# List of stars for which to compute the parameters.')
PARSER.add_argument('-inst', default='all',\
                    help='Names of the instruments you want to use '\
                           'separated by commas. If "all" is included, '\
                           'spectra from all accepted instruments will be used')
PARSER.add_argument('-ncores', default=3, type=int,\
                    help='Number of cores to use.')
PARSER.add_argument('-output_format', default='Both',\
                    choices=['Both', 'ASCII', 'FITS'],\
                    help='Format of the output file.')
PARSER.add_argument('-output', default='summary_output',\
                    help='Name of the output file, without extension.')
PARSER.add_argument('-hold_params', action='store_true', default=False,\
                    help='Option to set the atmospheric parameters to a '\
                           'certain value.')
PARSER.add_argument('-file_params', default='any_file.txt',\
                    help='Name of the file with the parameters to hold. ' \
                           'Will be used only if hold_params == True.')
PARSER.add_argument('-hold_mass', action='store_true', default=False,\
                    help='Option to set the mass, age, and photometric logg '\
                           'to a certain value. \
                           DEPRECATED.')
PARSER.add_argument('-file_mass', default='any_file.txt',\
                    help='Name of the file with the hold values for the mass, '\
                           'age and photometric logg. '\
                           'Will only be used if hold_mass == True. \
                           DEPRECATED.')
PARSER.add_argument('-set_boundaries', action='store_true', default=False,\
                    help='Option to change the boundaries of the '\
                           'atmospheric parameters.')
PARSER.add_argument('-file_with_boundaries', default='any_file.txt',\
                    help='Name of the file with the boundaries for the '\
                           'atmospheric parameters. '\
                           'Will only be used if set_boundaries == True.')
PARSER.add_argument('-no_abundance', action='store_true', default=False, \
                    help='Option to compute abundances for a set of elements.')
PARSER.add_argument('-no_error', action='store_true', default=False,\
                    help='Option to compute uncertainties to the '\
                           'atmospheric parameters.')
PARSER.add_argument('-no_debug', action='store_true', default=False,\
                    help='Option to create a debug file for each star, '\
                           'where the details of the computation will be stored.')
PARSER.add_argument('-no_check_t', action='store_true', default=False,\
                    help='Check the temperature with the one computed using '\
                           'the formulas from Casagrande et al. 2012.')
PARSER.add_argument('-file_colors', default='any_file.txt',\
                    help='Name of the file with the colors to use ' \
                           'to compute the temperature from Casagrande et al. ')
PARSER.add_argument('-use_coords', action='store_true', default=False,\
                    help='Option to use the coordinates of the star, '\
                           'instead of its name, to look for the star colors. '\
                           'Used only if check_t == True')
PARSER.add_argument('-file_coords', default='any_file.txt',\
                    help='Name of the file wiith the coordinates '\
                           'of the stars (in deg). '\
                           'Used only if check_t == True and use_coords == True.')

PARSER.add_argument('-mass_from_file', action='store_true', default=False,\
                    help='Load the chains created during a previous mass calculation,\
                            stored in dotter_isochrones/*.h5')

PARSER.add_argument('-no_vt_hold', action='store_true', default=False,\
                    help='Do not try to run the computation with vt = 1.2 km/s in the case it \
                            does not converge. Default is False.')

PARSER.add_argument('-no_r_iso_logg', action='store_true', default=False,\
                    help='Do not use the logg from isochrones interpolation \
                            to recompute the atmospheric parameters when there is \
                            disagreement with the spectroscopic one.')
PARSER.add_argument('-age_limits',\
                    help='Set limits to the parameters space for the stellar age. \
                            It has to be in the format log(age), with age in years \
                            (not Gyr). Default is None.')
PARSER.add_argument('-path_spectra', default='./Spectra',\
                    help='Path to the directory containing the spectra. It has \
                            to be relative to the path of SPECIES, or the \
                            absolute path.')

PARSER.add_argument('-columns', default='All',\
                    help='Columns to include in the output. The name of the star, \
                            instrument and exception will always be included.')

PARSER.add_argument('-giant', default=False, action='store_true',\
                    help='If the star, or list of stars, are considered to be giants.')

PARSER.add_argument('-path_iso', default='./isochrones',\
                    help='Path to the directory where all the files related to \
                            the mass computation will be stored. It has \
                            to be relative to the path of SPECIES, or the \
                            absolute path.')

PARSER.add_argument('-path_ew', default='./EW/plots_EW',\
                    help='Path to the directory where the plots of the equivalent \
                          width for all the lines will be stored. It has \
                            to be relative to the path of SPECIES, or the \
                            absolute path.')

PARSER.add_argument('-skip_mass', default=False, action='store_true',\
                    help='Skip the mass, radius, age, photometric logg, and luminosity \
                            computation. Usually done to reduce computation time.')

PARSER.add_argument('-no_rv', default=False, action='store_true',\
                    help='Disable restframe correction. Default is False \
                          (restframe correction is performed).')

PARSER.add_argument('-skip_vsini', default=False, action='store_true',\
                    help='Skip the vsini and vmac computation.')

PARSER.add_argument('-only_EW', default=False, action='store_true',\
                    help='Compute only the equivalent widths, skipping the \
                          derivation of the stellar parameters.')

PARSER.add_argument('-linelist', default='linelist.dat')

PARSER.add_argument('-no_EW_recomp', default=False, action='store_true')

PARSER.add_argument('-linelist_ab', default='lines_ab.dat')

PARSER.add_argument('-no_EW_recomp_ab', default=False, action='store_true')

PARSER.add_argument('-refit', default=False, action='store_true')

PARSER.add_argument('-plot_mass', default=False, action='store_true')

PARSER.add_argument('-minimization', default='downhill_simplex')

PARSER.add_argument('-other_instrument', default=False, action='store_true')

PARSER.add_argument('-redo_restframe', default=False, action='store_true')

PARSER.add_argument('-read_mode', default='linearregression',\
                    choices=['linearregression', 'odr', 'moog'])

PARSER.add_argument('-make_1d', default=False, action='store_true')

PARSER.add_argument('-plotEW', default=False, action='store_true')


ARGS = PARSER.parse_args()


# List of stars for which to compute the parameters.
STARLIST = ARGS.starlist.split(',')
if 'any_star' in STARLIST:
    print('\t Starlist was not given in the command line. ' + \
          'Reading from inside SPECIES.py\n')


    STARLIST = ['ceres01','ceres02','ceres03','moon01','ganymede01',\
                'sun01','sun02','sun03','sun04','sun05']


INST = ARGS.inst.split(',')
if 'all' in INST:
    HARPS, FEROS, UVES, HIRES, CORALIE, AAT, PFS, LCONRES = [True, True, True, True, True, True, True, True]
else:
    HARPS = 'harps' in INST
    FEROS = 'feros' in INST
    UVES = 'uves' in INST
    HIRES = 'hires' in INST
    CORALIE = 'coralie' in INST
    AAT = 'aat' in INST
    PFS = 'pfs' in INST
    LCONRES = 'lconres' in INST

NCORES = ARGS.ncores
OUTPUT_FORMAT = ARGS.output_format
FILE_OUTPUT = ARGS.output
HOLD_PARAMS = ARGS.hold_params
FILE_PARAMS = ARGS.file_params
HOLD_MASS = ARGS.hold_mass
FILE_HOLD_MASS = ARGS.file_mass
MASS_FROM_FILE = ARGS.mass_from_file
SET_BOUNDARIES = ARGS.set_boundaries
FILE_WITH_BOUNDARIES = ARGS.file_with_boundaries
ABUNDANCE = not ARGS.no_abundance

ERROR = not ARGS.no_error
if ERROR:
    FILE_COEFS_ERROR = True
    PLOTS_COEFS_ERROR = True

SET_DEBUG = not ARGS.no_debug
COLORS_FROM_FILE = False
CHECK_T = not ARGS.no_check_t

NAME_FILE_WITH_COLORS = ARGS.file_colors
if (NAME_FILE_WITH_COLORS is not None) and (NAME_FILE_WITH_COLORS != 'any_file.txt'):
    COLORS_FROM_FILE = True

USE_COORDS = ARGS.use_coords
FILE_WITH_COORDS = ARGS.file_coords
TRY_WITH_VT_HOLD = not ARGS.no_vt_hold
RECOMPUTE_WITH_LOGG_P = not ARGS.no_r_iso_logg

AGE_LIMITS = ARGS.age_limits
if AGE_LIMITS is not None:
    a = AGE_LIMITS.split(',')
    AGE_LIMITS = (float(a[0]), float(a[1]))

PATH_SPECTRA = ARGS.path_spectra
PATH_ISO = ARGS.path_iso
PATH_EW = ARGS.path_ew
OUT_COLUMNS = ARGS.columns
IS_GIANT = ARGS.giant
SKIP_MASS = ARGS.skip_mass
RESTFRAME = not ARGS.no_rv
SKIP_VSINI = ARGS.skip_vsini
ONLY_EW = ARGS.only_EW

LINELIST = ARGS.linelist
EW_RECOMP = not ARGS.no_EW_recomp
LINELIST_AB = ARGS.linelist_ab
EW_RECOMP_AB = not ARGS.no_EW_recomp_ab
REFIT = ARGS.refit
PLOT_MASS = ARGS.plot_mass
MINIMIZATION = ARGS.minimization
OTHER_INSTRUMENT = ARGS.other_instrument
NEW_RES = ARGS.redo_restframe

if NEW_RES:
    REFIT = True
    
MAKE_1D = ARGS.make_1d
READ_MODE = ARGS.read_mode
PLOTEW = ARGS.plotEW

START = time.time()

#########################################################################
# Corrects the spectra to restframe
#########################################################################

print('\t Checking the spectra and correcting (or not) to restframe.')


STARNAME = spectra(STARLIST, HARPS=HARPS, FEROS=FEROS, UVES=UVES,\
                   HIRES=HIRES, CORALIE=CORALIE, AAT=AAT, PFS=PFS, LCONRES=LCONRES,\
                   s_path=PATH_SPECTRA,\
                   restframe=RESTFRAME,\
                   use_other=OTHER_INSTRUMENT, new_res=NEW_RES,\
                   make_1d=MAKE_1D)

if not STARNAME:
    print('\t No star with valid spectra were found. Stopping the computation.')

else:
    ntotal = len(STARNAME)
    print('\n\n\t The computation will be carried out for the following stars (%d): ' % ntotal)
    str_ = "\n\t\t\t"
    print('\t\t\t'+str_.join("{}".format(s) for s in STARNAME))
    
    if NCORES == 1 or ntotal == 1:
        ############################################################################
        # Will perform the computation without parallelization
        ############################################################################
        
        from CalcParams import run_iteration
        
        RESULTS = []
        for i_s, s in enumerate(STARNAME):
            try:
                result = run_iteration(s, SET_DEBUG, HOLD_PARAMS, ABUNDANCE,\
                                       ERROR, FILE_PARAMS, FILE_COEFS_ERROR, PLOTS_COEFS_ERROR,\
                                       CHECK_T, COLORS_FROM_FILE, NAME_FILE_WITH_COLORS,\
                                       USE_COORDS, FILE_WITH_COORDS,\
                                       SET_BOUNDARIES, FILE_WITH_BOUNDARIES, MASS_FROM_FILE,\
                                       TRY_WITH_VT_HOLD, RECOMPUTE_WITH_LOGG_P, AGE_LIMITS,\
                                       PATH_SPECTRA, IS_GIANT, PATH_ISO, SKIP_MASS, PATH_EW,\
                                       SKIP_VSINI, ONLY_EW, LINELIST, EW_RECOMP, REFIT, PLOT_MASS,\
                                       MINIMIZATION, READ_MODE, (i_s+1, ntotal), LINELIST_AB,\
                                       EW_RECOMP_AB, PLOTEW)
            except Exception as e:
                print('%s: ERROR IN CODE, STOPPING THE COMPUTATION' % s)
                logging.error('%s: ERROR IN CODE, STOPPING THE COMPUTATION', s)
                exc_type, exc_obj, exc_tb = sys.exc_info()
                logging.error('line %d: %s', exc_tb.tb_lineno, e)
                logging.error(exc_type)
                logging.error(exc_obj)
                logging.error(exc_tb)
                result = (s, 0.0, 0.0, 0.0, 0.0, 2, 0, 0, 0.0, 0.0,\
                          {}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
                          0.0, 0.0, 0.0, 'no', 'no', 'no', 0.0, 0.0, 'None', 0.0,\
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            RESULTS.append(result)
        RESULTS = np.array(RESULTS)
        

    else:
        ############################################################################
        # Creates the list of arguments to be used in the parallelization
        ############################################################################

        ARG_LIST = []
        for i_s, s in enumerate(STARNAME):
            ARG_LIST.append((s, SET_DEBUG, HOLD_PARAMS, ABUNDANCE,\
                            ERROR, FILE_PARAMS, FILE_COEFS_ERROR, PLOTS_COEFS_ERROR,\
                            CHECK_T, COLORS_FROM_FILE, NAME_FILE_WITH_COLORS,\
                            USE_COORDS, FILE_WITH_COORDS,\
                            SET_BOUNDARIES, FILE_WITH_BOUNDARIES, MASS_FROM_FILE,\
                            TRY_WITH_VT_HOLD, RECOMPUTE_WITH_LOGG_P, AGE_LIMITS,\
                            PATH_SPECTRA, IS_GIANT, PATH_ISO, SKIP_MASS, PATH_EW,\
                            SKIP_VSINI, ONLY_EW, LINELIST, EW_RECOMP, REFIT, PLOT_MASS,\
                            MINIMIZATION, READ_MODE, (i_s+1, ntotal), LINELIST_AB,\
                            EW_RECOMP_AB, PLOTEW))

        ############################################################################
        # Run the computation of the parameters, with the stars split into
        # different cores.
        ############################################################################

        POOL = mp.Pool(NCORES, maxtasksperchild=10)
        RESULTS = POOL.map(multi_run_wrapper, ARG_LIST, chunksize=1)
        RESULTS = np.array(RESULTS)
        POOL.close()

        I_NONE = np.where(RESULTS.T[0] is None)[0]
        RESULTS = np.delete(RESULTS, I_NONE, 0)

        del I_NONE, POOL

    ############################################################################
    # Writes the output in each variable
    ############################################################################

    STAR, T, LOGG, XMETAL, MICRO, EXCEPTION, NFEI, NFEII, AB_FEI, AB_FEII,\
        ABUND_DICT, ERR_T, ERR_T2, ERR_LOGG, ERR_MET, ERR_VT,\
        ERR_VT2, VS, ERR_VS, VM, ERR_VM, MASS, ERR_MASS1, ERR_MASS2, AGE,\
        ERR_AGE1, ERR_AGE2, LOGG_ISO, ERR_LOGG_ISO1, ERR_LOGG_ISO2, RADIUS,\
        ERR_RADIUS1, ERR_RADIUS2, LOGL, ERR_LOGL1, ERR_LOGL2, USE_TC, USE_VT, USE_LOGG_P,\
        T_C, ERR_T_C, RELATION_TEMP, A_V, EEP, ERR_EEP1, ERR_EEP2,\
        P_PREMS, P_MS, P_RGB, P_HB, P_POSTHB = RESULTS.transpose()

    ############################################################################
    # Options for output file
    ############################################################################

    if OUTPUT_FORMAT == 'FITS':
        SAVE_AS_ASCII = False
        SAVE_AS_FITS = True
    elif OUTPUT_FORMAT == 'ASCII':
        SAVE_AS_ASCII = True
        SAVE_AS_FITS = False
    else:
        SAVE_AS_ASCII = True
        SAVE_AS_FITS = True

    ############################################################################
    # Sets the instrument and rv columns
    ############################################################################

    NAME_STAR = STAR[:]
    INSTRUMENT = np.array(['nofound' for i in range(len(STAR))])
    RV = np.zeros(len(STAR))
    for i, s in enumerate(STAR):
        INDICES = [indx for indx, c in enumerate(s) if c == '_']
        if INDICES:
            NAME_INST = s[INDICES[-1]+1:]
            if NAME_INST == 'o':
                NAME_INST = s[INDICES[-2]+1:]
                NAME_STAR[i] = s[:INDICES[-2]]
            else:
                NAME_STAR[i] = s[:INDICES[-1]]
            if NAME_INST == 'harps':
                INSTRUMENT[i] = 'HARPS'
            elif NAME_INST == 'feros':
                INSTRUMENT[i] = 'FEROS'
            elif NAME_INST == 'feros_o':
                INSTRUMENT[i] = 'FEROS'
            elif NAME_INST == 'uves':
                INSTRUMENT[i] = 'UVES'
            elif NAME_INST in ('hires', 'HIRES'):
                INSTRUMENT[i] = 'HIRES'
            elif NAME_INST == 'coralie':
                INSTRUMENT[i] = 'CORALIE'
            elif NAME_INST == 'aat':
                INSTRUMENT[i] = 'AAT'
            elif NAME_INST == 'pfs':
                INSTRUMENT[i] = 'PFS'
            else:
                INSTRUMENT[i] = NAME_INST.upper()

        try:
            h = fits.getheader('Spectra/%s_%s_res.fits' % (NAME_STAR[i], NAME_INST), 0)
            RV[i] = h['RV']
            del h
        except (KeyError, IOError):
            pass

        if os.path.samefile(PATH_SPECTRA, './Spectra') is False:
            if os.path.isfile(os.path.join('./Spectra',\
                                           '%s_%s_res.fits' % (NAME_STAR[i], NAME_INST))):
                os.unlink(os.path.join('./Spectra', '%s_%s_res.fits' % (NAME_STAR[i], NAME_INST)))

    ############################################################################
    # Writes the results in the output files
    ############################################################################

    ORIGINAL_COLUMNS = np.array(['Starname', 'Instrument', 'RV', '[Fe/H]', 'err_[Fe/H]', \
    'Temperature', 'err_T', 'logg', 'err_logg', 'vt', 'err_vt', 'nFeI', 'nFeII',\
    'exception', 'vsini', 'err_vsini', 'vmac', 'err_vmac',\
    'Mass', 'err_mass1', 'err_mass2', 'Age', 'err_age1', 'err_age2', 'iso_logg', \
    'err_iso_logg1', 'err_iso_logg2', 'Radius', 'err_radius1', 'err_radius2', \
    'logL', 'err_logL1', 'err_logL2', 'eep', 'err_eep1', 'err_eep2',\
    'P_preMS', 'P_MS', 'P_RGB', 'P_HB', 'P_postHB', 'use_Tc', 'use_vt', 'use_iso_logg', \
    'T_photo', 'err_T_photo', 'T_photo_relation', 'Av'])
            

    DATA = Table([NAME_STAR, INSTRUMENT, RV, XMETAL, ERR_MET, T, ERR_T,\
            LOGG, ERR_LOGG, MICRO, ERR_VT, NFEI, NFEII, EXCEPTION, VS, ERR_VS,\
            VM, ERR_VM, \
            MASS, ERR_MASS1, ERR_MASS2, AGE, ERR_AGE1, ERR_AGE2, LOGG_ISO, ERR_LOGG_ISO1,\
            ERR_LOGG_ISO2, RADIUS, ERR_RADIUS1, ERR_RADIUS2, LOGL, ERR_LOGL1, ERR_LOGL2,\
            EEP, ERR_EEP1, ERR_EEP2, P_PREMS, P_MS, P_RGB, P_HB, P_POSTHB,\
            USE_TC, USE_VT, USE_LOGG_P, T_C, ERR_T_C, RELATION_TEMP, A_V],\
            names=['Starname', 'Instrument', 'RV', '[Fe/H]', 'err_[Fe/H]', \
            'Temperature', 'err_T', 'logg', 'err_logg', 'vt', 'err_vt', 'nFeI', 'nFeII',\
            'exception', 'vsini', 'err_vsini', 'vmac', 'err_vmac',\
            'Mass', 'err_mass1', 'err_mass2', 'Age', 'err_age1', 'err_age2', 'iso_logg', \
            'err_iso_logg1', 'err_iso_logg2', 'Radius', 'err_radius1', 'err_radius2', \
            'logL', 'err_logL1', 'err_logL2', 'eep', 'err_eep1', 'err_eep2',\
            'P_preMS', 'P_MS', 'P_RGB', 'P_HB', 'P_postHB', 'use_Tc', 'use_vt', 'use_iso_logg', \
            'T_photo', 'err_T_photo', 'T_photo_relation', 'Av'],\
            dtype=['str', 'str', 'float', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', 'float', 'int', \
            'int', 'int', 'float', 'float', 'float', 'float', \
            'float', 'float', 'float', \
            'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',\
            'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',\
            'float', 'float', 'str', 'str',\
            'str', 'float', 'float', 'str', 'float'])

    DATA['[Fe/H]'].format = '%8.2f'
    DATA['Temperature'].format = '.1f'
    DATA['Mass'].format = DATA['Age'].format = DATA['Radius'].format = DATA['logL'].format = '%8.2f'
    DATA['logg'].format = DATA['vt'].format = \
                DATA['vsini'].format = DATA['err_vsini'].format = \
                DATA['vmac'].format = DATA['err_vmac'].format = \
                DATA['err_T'].format = \
                DATA['err_logg'].format = \
                DATA['err_[Fe/H]'].format = \
                DATA['err_vt'].format = '%.3f'
    DATA['err_mass1'].format = DATA['err_mass2'].format = \
                DATA['err_age1'].format = DATA['err_age2'].format = DATA['iso_logg'].format = \
                DATA['err_iso_logg1'].format = DATA['err_iso_logg2'].format = \
                DATA['err_radius1'].format = DATA['err_radius2'].format = \
                DATA['err_logL1'].format = DATA['err_logL2'].format = '%.3f'
    DATA['RV'].format = '%.3f'
    DATA['T_photo'].format = '%.1f'
    DATA['err_T_photo'].format = '%.3f'
    DATA['Av'].format = '%.3f'
    DATA['eep'].format = DATA['err_eep1'].format = DATA['err_eep2'].format = '%.2f'
    DATA['P_preMS'].format = DATA['P_MS'].format = DATA['P_RGB'].format = \
                DATA['P_HB'].format = DATA['P_postHB'].format = '%.3f'
    
    if ABUNDANCE:
        nameions = np.unique(ascii.read('Spectra/%s' % LINELIST_AB, comment='-')['ele'])
        data_ab = Table()
        for n in nameions:
            if n in ['FeI', 'FeII']:
                data_ab.add_column(Column([-99 for _ in STAR],\
                                          name='%s_ab' % n, dtype='float', format='%8.2f'))
                data_ab.add_column(Column([-99 for _ in STAR],\
                                          name='e_%s_ab' % n, dtype='float', format='%.3f'))
                data_ab.add_column(Column([0 for _ in STAR], name='n%s_ab' % n, dtype='int'))
            else:
                data_ab.add_column(Column([-99 for _ in STAR], name=n,
                                          dtype='float', format='%8.2f'))
                data_ab.add_column(Column([-99 for _ in STAR],
                                          name='e_%s' % n, dtype='float', format='%.3f'))
                data_ab.add_column(Column([0 for _ in STAR],
                                          name='n%s' % n, dtype='int'))
        
        for i, s in enumerate(STAR):
            abdict = ABUND_DICT[i]
            for n in nameions:
                if n in abdict:
                    if n in ['FeI', 'FeII']:
                        data_ab['%s_ab' % n][i] = abdict[n][0]
                        data_ab['e_%s_ab' % n][i] = abdict[n][1]
                        data_ab['n%s_ab' % n][i] = abdict[n][2]
                    else:
                        data_ab['%s' % n][i] = abdict[n][0]
                        data_ab['e_%s' % n][i] = abdict[n][1]
                        data_ab['n%s' % n][i] = abdict[n][2]
            
        DATA = hstack([DATA, data_ab])

    DATA_SUMMARY = DATA['Starname', '[Fe/H]', 'err_[Fe/H]', 'Temperature', 'err_T', 'logg',\
                        'err_logg', 'logL', 'err_logL1', 'err_logL2', 'vsini', 'err_vsini',\
                        'Mass', 'err_mass1', 'err_mass2',\
                        'exception', 'use_Tc', 'use_vt', 'use_iso_logg', 'Instrument']

    if OUT_COLUMNS != 'All':
        COLUMNS = np.array(['Starname', 'Instrument', 'exception'] + OUT_COLUMNS.split(','))
        NOT_INCLUDED = ORIGINAL_COLUMNS[np.isin(ORIGINAL_COLUMNS, COLUMNS, invert=True)]
        DATA.remove_columns(NOT_INCLUDED)

    if SAVE_AS_ASCII:
        ascii.write(DATA, './output/%s_full.dat' % FILE_OUTPUT, \
                    format='fixed_width', delimiter=None, overwrite=True)
        ascii.write(DATA_SUMMARY, './output/%s.dat' % FILE_OUTPUT,\
                    format='fixed_width', delimiter=None, overwrite=True)

    if SAVE_AS_FITS:
        DATA.write('./output/%s.fits' % FILE_OUTPUT, format='fits', overwrite=True)


TIME_SECONDS = (time.time() - START)

print('\n\n\t Time in seconds: %.4f' % TIME_SECONDS)
print('\t Time in minutes: %.4f' % (TIME_SECONDS/60.))
print('\t Time in hours:   %.4f' % (TIME_SECONDS/3600.))

print('\n\t*****************************************\n')
print('\t Thank you for using SPECIES\n')
print('\t*****************************************\n')
