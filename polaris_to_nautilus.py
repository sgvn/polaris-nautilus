#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: polaris_to_nautilus.
authors: Sacha Gavino, adapted from Julia Kobus's routines
date: Jan 2021
language: PYTHON 3.7
SHORT DESCRIPTION:  Read POLARIS outputs to compute input flux and temperatures in the correct format for NAUTILUS. Stores the files both in NAUTILUS model and in POLARIS model   
"""
import os
import sys
import argparse
import struct
import shutil

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from astropy.constants import h, c

from config import * # import configuration file



#--- FUNCTIONS ---
#                #
#-----------------
def read_from_binary(input_file,n_byte):
    if n_byte == 2:
        out = struct.unpack('H',input_file.read(2))[0]
    elif n_byte == 8:
        out = struct.unpack('d',input_file.read(8))[0]
        
    return out
    
    
def read_polaris(grid_path):
	with open(grid_path, 'rb') as f:
		
		temp_ID = 2		# ID for dust temp: 2
		den_ID = 29		# ID for dust mass density: 29
		rf_ID = 33		# # ID for radiation field absolute value: 33
		
		grid_id = read_from_binary(f,2)
		n_quant = read_from_binary(f,2)
		
		quants = np.zeros(n_quant,dtype=int)
		for i_quant in range(n_quant):
			quants[i_quant] = read_from_binary(f,2)
		
		unique, counts = np.unique(quants, return_counts=True)
		quant_dict = dict(zip(unique, counts))
		
		n_gas_density = quant_dict[28]
		n_dust_density = quant_dict[den_ID]
		n_temp = quant_dict[temp_ID]
		if rf_ID in quant_dict:
			n_rad_lambdas = quant_dict[rf_ID]
		else:
			n_rad_lambdas = 0
		
		if n_gas_density != n_dust_density:
			print('Number of gas densities do not match number of dust densities')
			exit()
		
		rmin = read_from_binary(f,8)
		rmax = read_from_binary(f,8)
		
		if grid_id == 30:
			n_r = read_from_binary(f,2)
			n_phi = read_from_binary(f,2)
			n_theta = read_from_binary(f,2)
			
			temperature = np.zeros((n_r+1,n_phi,n_theta,n_temp))
			density = np.zeros((n_r+1,n_phi,n_theta,n_dust_density))
			rad_field_abs = np.zeros((n_r+1,n_phi,n_theta, n_rad_lambdas))
			
			sf_r = read_from_binary(f,8)
			sf_phi = read_from_binary(f,8)
			sf_theta = read_from_binary(f,8)
			
			#~ r_coords_bounds = np.zeros(n_r+1)
			#~ r_coords_mid = np.zeros(n_r)
			
			#~ dr1 = (rmax - rmin) * (sf_r-1.0)/ (sf_r**n_r - 1.0)
			#~ r_coords_bounds[0] = rmin
			#~ # set outer ring radii:
			#~ # r_coords_bounds( 1 ): r_in + dr1
			#~ # r_coords_bounds(n_r): r_ou
			#~ for i_r in range(1,n_r+1):
				#~ r_coords_bounds[i_r] = r_coords_bounds[0] + dr1 * (sf_r**i_r - 1.0) / (sf_r-1.0)
				#~ r_coords_mid[i_r-1] = (r_coords_bounds[i_r]+r_coords_bounds[i_r-1]) / 2
			
			for i_r in range(n_r):
				for i_phi in range(n_phi):
					for i_theta in range(n_theta):
						i_dens = 0
						i_temp = 0
						i_rad = 0
						for i_quant in range(n_quant):
							# ID for dust temp: 2
							if quants[i_quant]==temp_ID:
								temperature[i_r,i_phi,i_theta,i_temp] = read_from_binary(f,8)
								i_temp += 1
							
							#~ # ID for gas mass density: 28
							#~ elif quants[i_quant]==28:
								#~ density[i_r,i_phi,i_theta,i_dens] = read_from_binary(f,8)
								#~ i_dens += 1
								
							# ID for dust mass density: 29	
							elif quants[i_quant]==den_ID:
								density[i_r,i_phi,i_theta,i_dens] = read_from_binary(f,8)
								i_dens += 1
							
							# ID for radiation field absolute value: 33
							elif quants[i_quant]==rf_ID:
								rad_field_abs[i_r,i_phi,i_theta,i_rad] = read_from_binary(f,8)
								i_rad += 1
							
							else:
								tmp = read_from_binary(f,8)
			
			for i_quant in range(n_quant):
				# ID for dust temp: 2
				if quants[i_quant]==temp_ID:
					temperature[n_r,:,:] = read_from_binary(f,8)
				
				#~ # ID for gas mass density: 28
				#~ elif quants[i_quant]==28:
					#~ density[n_r,:,:] = read_from_binary(f,8)
				# ID for dust mass density: 29
				elif quants[i_quant]==den_ID:
					density[n_r,:,:] = read_from_binary(f,8)
				
				# ID for gas mass density: 33
				elif quants[i_quant]==rf_ID:
					rad_field_abs[n_r,:,:] = read_from_binary(f,8)
				
				else:
					tmp = read_from_binary(f,8)
    
    
		elif grid_id == 40:
			zmax = read_from_binary(f,8)
			
			n_rho = read_from_binary(f,2)
			n_phi = read_from_binary(f,2)
			n_z = read_from_binary(f,2)
		
			temperature = np.zeros((n_rho+n_z,n_phi,n_z))
			density = np.zeros((n_rho+n_z,n_phi,n_z))
			rad_field_abs = np.zeros((n_rho+n_z,n_phi,n_z))
			
			sf_r = read_from_binary(f,8)
			sf_phi = read_from_binary(f,8)
			sf_z = read_from_binary(f,8)
			
			
			for i_rho in range(n_rho):
				for i_phi in range(n_phi):
					for i_z in range(n_z):
						for i_quant in range(n_quant):
							# ID for dust temp: 2
							if quants[i_quant]==temp_ID:
								temperature[i_rho,i_phi,i_z] = read_from_binary(f,8)
							
							#~ # ID for gas mass density: 28
							#~ elif quants[i_quant]==28:
								#~ density[i_rho,i_phi,i_z] = read_from_binary(f,8)
							# ID for dust mass density: 29
							elif quants[i_quant]==den_ID:
								density[i_rho,i_phi,i_z] = read_from_binary(f,8)
							
							# ID for radiation field absolute value: 33
							elif quants[i_quant]==rf_ID:
								rad_field_abs[i_rho,i_phi,i_z] = read_from_binary(f,8)
							
							else:
								tmp = read_from_binary(f,8)
			
			for i_z in range(n_z):
				for i_quant in range(n_quant):
					# ID for dust temp: 2
					if quants[i_quant]==temp_ID:
						temperature[n_rho+i_z,:,:] = read_from_binary(f,8)
					
					#~ # ID for gas mass density: 28
					#~ elif quants[i_quant]==28:
						#~ density[n_rho+i_z,:,:] = read_from_binary(f,8)
					# ID for dust mass density: 29
					elif quants[i_quant]==den_ID:
						density[n_rho+i_z,:,:] = read_from_binary(f,8)
					
					# ID for gas mass density: 33
					elif quants[i_quant]==rf_ID:
						rad_field_abs[n_rho+i_z,:,:] = read_from_binary(f,8)
					
					else:
						tmp = read_from_binary(f,8)

	return temperature, density, rad_field_abs, rmin / AU, rmax / AU, n_r, sf_r, n_theta, sf_theta

def in_model_space(r, z, rmin, rmax):
	"""
	Test if the given position is inside model space
	"""
	
	radius = np.sqrt(r**2 + z**2)
	ims = (radius >= rmin) & (radius <= rmax)
	
	flag = False
	for num, is_in_model_space in enumerate(ims):
		if not is_in_model_space:
			print('coordinate (' + str(r[num]) + '/' + str(z[num]) + ') not in model space')
			
def calc_cell_number(r, z, rmin, rmax, n_r, n_theta, sf_r, sf_theta):		
	# r cell
	r_coords_bounds = np.zeros(n_r+1)
	r_coords_mid = np.zeros(n_r)
	dr1 = (rmax - rmin) * (sf_r-1.0)/ (sf_r**n_r - 1.0)
	r_coords_bounds[0] = rmin
	for i_r in range(1,n_r+1):
		r_coords_bounds[i_r] = r_coords_bounds[0] + dr1 * (sf_r**i_r - 1.0) / (sf_r-1.0)
		r_coords_mid[i_r-1] = (r_coords_bounds[i_r]+r_coords_bounds[i_r-1]) / 2
		
	#theta cells
	theta_coords_bounds = np.zeros(n_theta+1)
	for i_theta in range(0, int(n_theta/2)+1):
		theta_coords_bounds[i_theta] = np.pi/2 - (sf_theta**(n_theta / 2 - i_theta) -1)*np.pi/2 /(sf_theta**(n_theta / 2) -1)
	for i_theta in range(int(n_theta/2+1), n_theta + 1):
		theta_coords_bounds[i_theta] = np.pi/2 + (sf_theta**(i_theta - n_theta / 2) -1)*np.pi/2 /(sf_theta**(n_theta / 2) -1)
	theta_coords_mid = (theta_coords_bounds[1:] + theta_coords_bounds[:-1]) / 2	
	
	#Check if all coordinates are in model space:
	in_model_space(r, z, rmin, rmax)
	
	#find cell number:	
	n_cell = []
	midpoint_r = []
	midpoint_z =[]
	radius = np.sqrt(r**2 + z**2)
	theta = np.arctan2(z,r) + np.pi/2
	for rval, theta_val in zip(radius, theta):
		if rval>=r_coords_bounds[-1]:
			nr = n_r
		else:
			nr = np.amin(np.argwhere(np.array(r_coords_bounds)>rval))
			
		if theta_val>=theta_coords_bounds[-1]:
			ntheta = n_theta
		else:
			ntheta = np.amin(np.argwhere(np.array(theta_coords_bounds)>theta_val))
		midpoint_r.append(r_coords_mid[nr-1] * np.cos(theta_coords_mid[ntheta-1] - np.pi/2))
		midpoint_z.append(r_coords_mid[nr-1] * np.sin(theta_coords_mid[ntheta-1] - np.pi/2))
		n_cell.append((nr-1)*n_theta+ntheta - 1)
	n_cell = np.array(n_cell)
	midpoint_r = np.array(midpoint_r)
	midpoint_z = np.array(midpoint_z)
	return n_cell, midpoint_r, midpoint_z, nr, ntheta		
	
def savedensity(density_path, n_cell, r, z, midpoint_r, midpoint_z, density_data):
	if len(n_cell)==1:
		for i in range(n_species):
			print('grain size interval ' + str(i+1) + ' : den = ' + str(density_data[n_cell,i]) + ' kg m^-3')
		print("cell midpoint: r = " + str(midpoint_r[0]) + " au, z = " + str(midpoint_z[0]) + " au")
	else:
		den_array = np.zeros((len(n_cell), n_species + 4))
		den_array[:,0]=r
		den_array[:,1]=z
		den_array[:,2]=midpoint_r
		den_array[:,3]=midpoint_z
		den_array[:,4:] = density_data[n_cell,:]
		np.savetxt(density_path, den_array, header="r[au]    z[au]  nearest cell r[au] nearest cell z[au]  number density of grain size interval 1 [kg m^-3]    ...    number density of grain size interval 16 [m^-3]")

def savetemp(temp_path, n_cell, r, z, midpoint_r, midpoint_z, temp_data):
	if len(n_cell)==1:
		for i in range(n_species):
			print('grain size interval ' + str(i+1) + ' : T = ' + str(temp_data[n_cell,i]) + ' K')
		print("cell midpoint: r = " + str(midpoint_r[0]) + " au, z = " + str(midpoint_z[0]) + " au")
	else:
		temp_array = np.zeros((len(n_cell), n_species + 4))
		temp_array[:,0]=r
		temp_array[:,1]=z
		temp_array[:,2]=midpoint_r
		temp_array[:,3]=midpoint_z
		temp_array[:,4:]=temp_data[n_cell,:]
		np.savetxt(temp_path, temp_array, header="r[au]    z[au]	nearest cell r[au] nearest cell z[au]    temperature of grain size interval 1 [K]    ...    temperature of grain size interval 16 [K]")	
		
			
def saveradfield(path_to_wv, wavelengths, radfield, n_cell, r, z, midpoint_r, midpoint_z, rad_field, scaling = None):
	reshape_flux=pd.DataFrame()
	wave = np.array([
					7.692e-08,
					8.320e-08,
					9.000e-08,
					9.735e-08,
					1.053e-07,
					1.139e-07,
					1.232e-07,
					1.333e-07,
					1.442e-07,
					1.559e-07,
					1.687e-07,
					1.824e-07,
					1.973e-07,
					2.134e-07,
					2.309e-07,
					2.497e-07,
					2.701e-07,
					2.922e-07,
					3.161e-07,
					3.419e-07,
					3.698e-07,
					4.000e-07
					])
	for wavelength in wavelengths:
		if wavelength in wave:           
			waveindex = (np.abs(wave-wavelength)).argmin()		
			rad_field_s = rad_field[n_cell,waveindex]
		else:
			index_1 = np.where(wave < wavelength)[0][-1]
			index_2 = np.where(wave > wavelength)[0][0]
			rad_field_1 = rad_field[n_cell,index_1]
			rad_field_2 = rad_field[n_cell,index_2]
			rad_field_s = interp1d([wave[index_1],wave[index_2]], np.vstack([rad_field_1, rad_field_2]), axis=0)(wavelength)
		if scaling is not None:
			scaling_data = np.loadtxt(scaling)
			wavelengths_scaling = scaling_data[:,0] * 1e-9
			factor = scaling_data[:,1]
			scalingfactor = interp1d(wavelengths_scaling, factor)(wavelength)
			print("scaling", scalingfactor)
			rad_field_s = rad_field_s * scalingfactor 

		if len(n_cell)==1:
		    print(radfield +' RF = ' + str(rad_field_s) + 'W m^-1 m^-2', "for wavelength = ", wavelength, "m")
		    print("cell midpoint: r = " + str(midpoint_r[0]) + " au, z = " + str(midpoint_z[0]) + " au")      
		else:
		    rad_array = np.zeros((len(n_cell), 5))
		    rad_array[:,0]=r
		    rad_array[:,1]=z
		    rad_array[:,2]=midpoint_r
		    rad_array[:,3]=midpoint_z
		    rad_array[:,4]=rad_field_s
		    wavelength_nm = wavelength*1e+9
		    try:
		        os.mkdir(path_to_wv+'%.3f' % (wavelength_nm))
		    except FileExistsError:
		        pass
		    np.savetxt(path_to_wv+'%.3f/%s_rad_field_%.3em.txt' % (wavelength_nm ,radfield, wavelength), rad_array, header="r[au]    z[au]	nearest cell r[au] nearest cell z[au]    radiation field W m^-1 m^-2")
		
		#SAVE FLUX IN NAUTILUS MODEL
		name = ['r', 'z', 'near_r', 'near_z', 'flux']
		rad_array = pd.DataFrame(rad_array)
		rad_array.columns = name
		naut_flux = rad_array['flux']*1e-13*(float(wavelength)/(h*c))
		reshape_flux = pd.concat([reshape_flux, naut_flux], axis=1)


	return reshape_flux.transpose()

def read_coordinates_from_file(coordinates, radius):
    """
    read requested coordinates from file
    """
    coords = np.loadtxt(coordinates, comments="!")
	#r = coords[0]
    z = coords[:, 0]
    spt = len(z)
    r = np.full(spt, radius)
	
    return r, z


#_________________________________________________#
#                      MAIN                       #
#_________________________________________________#
if __name__ == "__main__":
	#--- PARSER ---
	parser = argparse.ArgumentParser(description='COUPLING ROUTINE FOR NAUTILUS AND POLARIS. Reads POLARIS binary files and extract flux values suited for NAUTILUS. WARNING!: give only the file names, not the paths')
	
	parser.add_argument('-t', '--temperaturefile', type=str,
						help='Enter the Polaris-file name (only the name), from which density and temperature values are to be extracted. ')
						
	parser.add_argument('-r', '--radfieldfile', type=str,
						help='Enter the Polaris-file name (only the name), from which the stellar radiation field is to be extracted. ')
						
	parser.add_argument('-isrf', '--isrffile', type=str,
						help='Enter the Polaris-file name (only the name), from which the interstellar radiation field is to be extracted. ')
													
	parser.add_argument('--scaling', dest='scaling',  type=str, default=None,
                                                help='Optional: enter file with scaling factors for scaling the stellar radition field, WARNING: wavelengths in this file (left column) must be in nm')
	args = parser.parse_args()
	#------------


	if args.temperaturefile is None and args.radfieldfile is None and args.isrffile is None:
		print('\n')
		print('---------------------------------------------------------------------------------------------------------------------')
		print("Please enter at least one polaris binary file.") 
		print("The binary files are output from polaris.")
		print("Optional: you can add a scaling factor for the stellar flux. See help.")
		#print("The wavelengths file contains a list of waves [m] and must be located in the same folder as this routine.\n")
		print("Here is an example to write in a command line:")
		print("$python polaris_to_nautilus.py -t grid_Temp.dat -r grid_SRF.dat -isrf grid_ISRF.dat \n")
		print("If you want to get extra help, write in the command line something like that:")
		print("$python polaris_to_nautilus.py -h")
		print('---------------------------------------------------------------------------------------------------------------------')
		print('\n')
	else:

		#--Create subfolder in the polaris results model--
		if not os.path.exists(path_to_polaris+disk):
			os.makedirs(path_to_polaris+disk)
		else:
			shutil.rmtree(path_to_polaris+disk)           # Removes all the subdirectories!
			os.makedirs(path_to_polaris+disk)

		#--Create subfolder in the polaris results model AND check if already exist.--
		#if os.path.exists(path_to_polaris+disk): #check if the folder exists, if yes, exit
		#	print('WARNING: the subfolder "%s" already exists in: %s.' %(disk, path_to_polaris))
		#	print('Please choose another name or remove "%s".' %disk)
		#	sys.exit(0)
		#else:
		#	os.makedirs(path_to_polaris+disk)

		radii = os.listdir ("%s" %path_to_nautilus)
		#radii = np.asarray(next(os.walk("%s" %path_to_nautilus))[1])
		radii.sort()
		for radius in radii:
			if radius.endswith('AU'): # loop through all the files and folders
				print("radius: ", radius)
				rad = radius.replace('AU','')
				rad = float(rad)
				radius = radius + "/"

				#----> LOOP STARTS HERE:<-------
				reshape_isrf=pd.DataFrame()
				reshape_srf=pd.DataFrame()
				reshape_tot=pd.DataFrame()

				os.makedirs(path_to_polaris+disk+radius)

				if not os.path.exists(path_to_nautilus+radius+flux):
					os.makedirs(path_to_nautilus+radius+flux)
				else:
					shutil.rmtree(path_to_nautilus+radius+flux)           # Removes all the subdirectories!
					os.makedirs(path_to_nautilus+radius+flux)

				r, z = read_coordinates_from_file(path_to_nautilus+radius+static, rad)


				print("Read Polaris-files")
				if args.temperaturefile is not None:
					temperature, density, rf , rmin, rmax, n_r, sf_r, n_theta, sf_theta = read_polaris(path_to_polaris+args.temperaturefile)
				if args.radfieldfile is not None:
					t, d, rad_field_abs , rmin, rmax, n_r, sf_r, n_theta, sf_theta = read_polaris(path_to_polaris+args.radfieldfile)
				if args.isrffile is not None:
					t, d, is_rad_field_abs , rmin, rmax, n_r, sf_r, n_theta, sf_theta = read_polaris(path_to_polaris+args.isrffile)


		
				print("Calc cell numbers")
				n_cell, midpoint_r, midpoint_z, nr, ntheta = calc_cell_number(r, z, rmin, rmax, n_r, n_theta, sf_r, sf_theta) 
	

				if args.temperaturefile is not None:
					temperature = temperature[:,:,:,:16].reshape((n_r+1)*n_theta, 16)
					savetemp(path_to_polaris+disk+radius+temperature_file, n_cell, r, z, midpoint_r, midpoint_z, temperature) 

		
				if args.radfieldfile is not None:
					rad_field = rad_field_abs.reshape((n_r+1)*n_theta, len(rad_field_abs[0,0,0,:]))
					flux_srf = saveradfield(path_to_polaris+disk+radius, wavelengths, stellar, n_cell, r, z, midpoint_r, midpoint_z, rad_field, args.scaling)
					flux_srf.to_csv(path_to_nautilus+radius+flux+srf_flux, sep=' ', float_format='%.5e', index=False, header=False)  

		
				if args.isrffile is not None:
					is_rad_field = is_rad_field_abs.reshape((n_r+1)*n_theta, len(is_rad_field_abs[0,0,0,:]))
					flux_isrf = saveradfield(path_to_polaris+disk+radius, wavelengths, interstellar, n_cell, r, z, midpoint_r, midpoint_z, is_rad_field) 
					flux_isrf.to_csv(path_to_nautilus+radius+flux+isrf_flux, sep=' ', float_format='%.5e', index=False, header=False)

				#total flux:
				if args.isrffile is not None and args.radfieldfile is not None:
					flux_tot = flux_srf + flux_isrf
					flux_tot.to_csv(path_to_nautilus+radius+flux+total_flux, sep=' ', float_format='%.5e', index=False, header=False)

				density = density.reshape((n_r+1)*n_theta, 16)
				savedensity(path_to_polaris+disk+radius+density_file, n_cell, r, z, midpoint_r, midpoint_z, density)
