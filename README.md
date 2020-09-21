# SPECIES

Spectroscopic Parameters and atmosphEric ChemIstriEs of Stars (SPECIES) is a code meant to compute stellar parameters and abundances by using high resolution echelle spectra. The whole calculation is done automatically, with the stellar spectrum being the only mandatory input. It handles data from several spectrographs (HARPS, FEROS, UVES, HIRES, PFS, CORALIE so far) and more than one star at the same time.

Please cite Soto & Jenkins 2018, http://adsabs.harvard.edu/abs/2018A%26A...615A..76S if your use SPECIES for your work.

Authors: Maritza Soto and James Jenkins.

(**SPECIES went through some major changes. If you installed SPECIES before September 2020, please update all the scripts and make sure you include the MOOGPATH in your bash file (see the wiki). It is not necessary to reinstall ARES nor MOOG**).

# Computation of parameters

The atmospheric parameters (temperature, metallicity, surface gravity and microturbulence velocity) are computed by measuring the equivalent widths of several iron lines, done using ARES (Sousa et al. 2008). These are then given to MOOG (Sneden 1973), which solves the radiative transfer equation assuming local thermodynamic equilibrium (LTE) conditions. The atmospheric parameters are then derived through an iterative process that stops when no correlation is found between the line abundances with the excitation potential and the equivalent width. The atmospheric models are obtained from interpolation through a grid of ATLAS9 models (Castelli & Kurucz 2004).

Chemical abundances are obtained for 11 elements: Na, Mg, Al, Si, Ca, Ti, Cr, Mn, Ni, Cu and Zn. These elements are the default ones included in our `Spectra/lines_ab.dat` linelist file, but that file can be modified or another linelist file can be use (please refer to the wiki --> Running SPECIES). Rotational and macroturbulence velocity are found by temperature relations, and line fitting, measuring the profiles of five absorption lines.

Finally, physical parameters like mass and age are computed by interpolating throught a grid of MIST evolutionary models, using the metallicity, temperature and surface gravity found previously, as well as their uncertainties, as constrains to the likelihood function. It uses a bayesian approach to obtain the final values, which are taken as the mean and standard deviation of the Gaussian profile adjusted to the resulting chains.

More detail about the method and results from SPECIES can be found in Soto & Jenkins 2018 (http://adsabs.harvard.edu/abs/2018A%26A...615A..76S) and Soto et al. 2020.

# Installation and use

SPECIES is written mostly in Python, with the exception of MOOG, written in fortran. Installation instructions for MOOG, as well as required packages, are found in the Wiki page. Please refer to the Wiki for usage instructions as well. 

Please contact me if you have any questions or issues when running SPECIES!
