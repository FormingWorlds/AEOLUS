#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:51:59 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import pandas as pd
from scipy import interpolate
import pathlib
import glob, re, os
import numpy as np

# String sorting not based on natsorted package
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def InterpolateStellarLuminosity(L_sun, AU, star_mass, time, mean_distance, albedo, Sfrac):

    
    # File name
    fname = str(pathlib.Path(__file__).parent.parent.absolute())+"/luminosity_tracks/Lum_m"+str(star_mass)+".txt"
    
    # If file exists, just interpolate needed time
    if os.path.isfile(fname):

        luminosity_df           = pd.read_csv(fname)
        star_age                = (time["star"])/1e+6   # Myr
        ages                    = luminosity_df["age"]*1e+3         # Myr
        luminosities            = luminosity_df["lum"]              # L_sol

        # Interpolate luminosity for current time
        interpolate_luminosity  = interpolate.interp1d(ages, luminosities)
        interpolated_luminosity = interpolate_luminosity([star_age])
        interpolated_luminosity = interpolated_luminosity[0]

    # Else: interpolate from 2D grid
    else:

        # Find all luminosity tracks
        lum_tracks = natural_sort(glob.glob(str(pathlib.Path(__file__).parent.parent.absolute())+"/luminosity_tracks/"+"Lum_m*.txt"))
        # Define data arrays for interpolation later on
        xy_age_mass = []
        z_lum       = []

        # Fill the arrays
        for lum_track in lum_tracks:

            # Read the specific data file
            luminosity_df    = pd.read_csv(lum_track)

            # Cut the string to get the mass of the star
            star_mass_lst    = float(lum_track[-7:-4])

            # Read out age and luminosity
            age_list         = list(luminosity_df["age"]*1e+3)
            luminosity_list  = list(luminosity_df["lum"])

            mass_list        = np.ones(len(age_list))*star_mass_lst

            # Fill the arrays
            zip_list = list(zip(age_list, mass_list))
            xy_age_mass.extend(zip_list)
            z_lum.extend(luminosity_list)

        # Bring arrays in numpy shape
        xy_age_mass = np.array(xy_age_mass)
        z_lum       = np.array(z_lum)

        # Define the interpolation grids
        grid_x, grid_y = np.mgrid[0.7:10000:100j, 0.1:1.4:100j]

        # Interpolate the luminosity from the 2D grid
        #print("dimension of xi (time['star']/1e+6, star_mass) = ", np.shape((time["star"]/1e+6, star_mass)))
        interpolated_luminosity = interpolate.griddata(xy_age_mass, z_lum, (time["star"]/1e+6, star_mass), method='linear', rescale=True)

    # Stellar constant, W m-2
    S_0    = interpolated_luminosity * L_sun / ( 4. * np.pi * (mean_distance*AU)**2. )

    # Scale instellation by fixed fraction
    S_0    = S_0 * Sfrac

    # Mean flux averaged over surface area, W m-2
    toa_heating             = ( 1. - albedo ) * S_0 / 4.
    
    return toa_heating