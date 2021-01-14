#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: config.
date: Jan 2021
language: PYTHON 3.7
DESCRIPTION: configuration file used by polaris_to_nautilus.py file. 
"""
import numpy as np

#------ Parameters ------
AU = 1.496e+11 # au in m
n_species = 16 # number of grain sizes

#------ wavelength list ------
wave_min = 9.120e-08 # smallest wave [m]
wave_max = 3.992e-07 # largest wave [m]
step = 0.1e-8 # step [m]
wavelengths = np.arange(wave_min, wave_max+step, step) # Define the wavelength list at which flux will be computed.


#------ Path to models ------
path_to_nautilus="../../results/nautilus/model1/" # Paths to nautilus disk model
path_to_polaris="../../results/polaris/model1/" # Paths to polaris disk model

#------ POLARIS ------
disk ="disk/" # subfolder in path_to_polaris to store flux, T and density with same coordinates as NAUTILUS model. Not used by NAUTILUS directly but these values must match with those inside NAUTILUS model. 
density_file="density.txt" # File name where densities are written
temperature_file="temperature.txt" # File name where temperatures are written 
stellar="stellar" # name of stellar flux
interstellar="interstellar" # name of interstellar flux


#------ NAUTILUS ------
static="1D_static.dat" # NAUTILUS 1D_static.dat. Only used to match NAUTILUS coordinates with POLARIS coordinates.
isrf_flux="flux_isrf_dust.txt" # File name of interstellar flux with format readable by NAUTILUS
srf_flux="flux_srf_dust.txt" # File name of stellar flux with format readable by NAUTILUS
total_flux="local_flux_dust.txt" # File name of interstellar+stellar flux with format readable by NAUTILUS
flux="flux/" #subfolder where the flux files are stored for NAUTILUS. (WARNING: always write the final slash / in the name)





