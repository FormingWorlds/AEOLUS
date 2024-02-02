#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:20:27 2023

@authors: 
Mark Hammond (MH)
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)

AEOLUS radiative-convective model, using SOCRATES for radiative-transfer.
"""

import matplotlib as mpl
mpl.use('Agg')

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.solve_pt import RadConvEqm
from modules.solve_pt import *
from modules.plot_flux_balance import plot_fluxes
from utils.socrates import CleanOutputDir

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start AEOLUS")

    start = t.time()
    ##### Settings

    # Planet 
    time = { "planet": 0., "star": 4e+9 } # yr,
    star_mass     = 0.1*1.0                 # M_sun, mass of star
    mean_distance = 0.0252*1.0                 # au, orbital distance
    pl_radius     = 1.1*6.371e6             # m, planet radius
    pl_mass       = 1.05*5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    T_surf        = 300.0                # K
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    # P_surf       = 100 * 1e5
    # vol_partial = {}
    # vol_mixing = { 
    #                 "CO2"  : 0.00417,
    #                 "H2O"  : 0.03,
    #                 "N2"   : 0.78084,
    #                 "H2"   : 0.03, 
    #                 "CH4"  : 0.000187, 
    #                 "O2"   : 0.20946, 
    #                 "O3"   : 0.0000006, 
    #                 "He"   : 0.00000524 , 
    #             }
    
    # OR:
    # Define volatiles by partial pressures
    P_surf = 0.0
    vol_mixing = {}
    vol_partial = {
        "H2O" : 0.0354 * 1.0e5, #ga.p_sat('H2O',T_surf), #1e3,
        "NH3" : 0.,
        "CO2" : 0.,
        "CH4" : 0.0,
        "CO" : 0.,
        "O2" : 0.0e4,
        "N2" : 1.0e5,
        "H2" : 0.0e3
        }

    # Stellar heating on/off
    stellar_heating = True

    # Rayleigh scattering on/off
    rscatter = True

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = False

    # Tropopause calculation
    trppD = True   # Calculate dynamically?
    trppT = 30.0     # Fixed tropopause value if not calculated dynamically

    # Water lookup tables enabled (e.g. for L vs T dependence)
    water_lookup = False
    
    # Surface temperature time-stepping
    surf_dt = False
    cp_dry = False
    # Options activated by surf_dt
    cp_surf = 1e5         # Heat capacity of the ground [J.kg^-1.K^-1]
    mix_coeff_atmos = 1e6 # mixing coefficient of the atmosphere [s]
    mix_coeff_surf  = 1e6 # mixing coefficient at the surface [s]

    # Cloud radiation
    do_cloud = False
    # Options activated by do_cloud
    re   = 1.0e-5 # Effective radius of the droplets [m] (drizzle forms above 20 microns)
    lwm  = 0.8    # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
    clfr = 0.8    # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Set up dirs
    dirs = {
            "aeolus": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+"/output/"
            }
    
    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Create atmosphere object
    atm = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, re, lwm, clfr,
                vol_mixing=vol_mixing, vol_partial=vol_partial, calc_cf=calc_cf, trppT=trppT, req_levels=100, water_lookup=water_lookup, do_cloud=do_cloud)

    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.toa_heating = 0.
    else:
        atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance)
        print("Instellation:", round(atm.toa_heating), "W/m^2")

    # Move/prepare spectral file
    print("Inserting stellar spectrum")

    StellarSpectrum.InsertStellarSpectrum(
        dirs["aeolus"]+"/spectral_files/Reach/Reach",
        dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]+"runtime_spectral_file"
    )

    # Set up atmosphere with general adiabat
    atm_dry, atm = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 

    OPR = atm.flux_up_total[0]
    ASR = atm.flux_down_total[0]
    imbalance = ASR - OPR

    print("T_surf, OLR, ASR, imbalance = ", T_surf, OPR, ASR, imbalance)
    
    # Plot abundances w/ TP structure
    if (cp_dry):
        ga.plot_adiabats(atm_dry,filename="output/dry_ga.pdf")
        atm_dry.write_PT(filename="output/dry_pt.tsv")
        plot_fluxes(atm_dry,filename="output/dry_fluxes.pdf")

    ga.plot_adiabats(atm,filename="output/moist_ga.pdf")
    atm.write_PT(filename="output/moist_pt.tsv")
    atm.write_ncdf("output/moist_atm.nc")
    plot_fluxes(atm,filename="output/moist_fluxes.pdf")

    # Copy this file to the output folder
    shutil.copy(__file__, os.path.join("output", os.path.basename(__file__)))

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    end = t.time()
    print("Runtime:", round(end - start,2), "s")

