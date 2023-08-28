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

import time as t
import os
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.radcoupler import RadConvEqm

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.SocRadModel as SocRadModel
from utils.atmosphere_column import atmos

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start AEOLUS")

    start = t.time()
    ##### Settings

    # Constants
    L_sun                   = 3.828e+26        # W, IAU definition
    AU                      = 1.495978707e+11  # m
    
    # Planet 
    time = { "planet": 0., "star": 4567e+6 } # yr,
    # time_current  = 0                 # yr, time after start of MO
    # time_offset   = 4567e+6           # yr, time relative to star formation
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 1.0                 # au, orbital distance
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass
    pl_albedo     = 0.3                 # [1], planet albedo

    # Boundary conditions for pressure & temperature
    T_surf        = 300.0                # K
    P_top         = 10.0                  # Pa

    # Define volatiles by mole fractions
    # P_surf       = 50 * 1e5
    # vol_mixing = { 
    #                 "CO2"  : 1.0 - 5e-2,
    #                 "H2O"  : 5e-2 - 1e-6,
    #                 "N2"   : 1e-6,  
    #                 "H2"   : 0.0, 
    #                 "NH3"  : 0.0,
    #                 "CH4"  : 0.0, 
    #                 "O2"   : 0.0, 
    #                 "CO"   : 0.0, 
    #                 # # No thermodynamic data, RT only
    #                 # "O3"   : 0.01, 
    #                 # "N2O"  : 0.01, 
    #                 # "NO"   : 0.01, 
    #                 # "SO2"  : 0.01, 
    #                 # "NO2"  : 0.01, 
    #                 # "HNO3" : 0.01, 
    #                 # "He"   : 0.01, 
    #                 # "OCS"  : 0.01,
    #             }
    # vol_partial = {}

    # OR:
    # Define volatiles by partial pressures
    P_surf = 1.0e5
    vol_mixing = {}
    vol_partial    = { 
                          "H2O" :  P_surf,
                          "NH3" :  0.,
                          "CO2" :  0.,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  0.,
                          "H2"  :  0.
                        }

    # Stellar heating on/off
    stellar_heating = True
    
    # False: interpolate luminosity from age and mass tables. True: define a custom instellation.
    custom_ISR = False

    # Rayleigh scattering on/off
    rscatter = True

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 70.0     # Fixed tropopause value if not calculated dynamically
    
    # Surface temperature time-stepping
    surf_dt = False
    cp_dry = False
    # Options activated by surf_dt
    cp_surf = 1e5         # Heat capacity of the ground [J.kg^-1.K^-1]
    mix_coeff_atmos = 1e6 # mixing coefficient of the atmosphere [s]
    mix_coeff_surf  = 1e6 # mixing coefficient at the surface [s]

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0
    
    # Database used
    database = 'ExoMol_010'                # 'ExoMol_resolution' or 'HITRAN_resolution'
    spec_bn  = 'sp_b318_'+database         # spectral basename

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, pl_albedo, vol_mixing=vol_mixing, vol_partial=vol_partial, calc_cf=calc_cf, trppT=trppT)

    # Compute stellar heating
    S_0, atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.toa_heating = 0.
    else:
        print("TOA heating:", round(atm.toa_heating), "W/m^2")
        
    # Compute heat flux
    dirs = {"output": os.getcwd()+"/output/"+database, "rad_conv": os.getcwd()}
    atm_dry, atm_moist = RadConvEqm(dirs, spec_bn, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 
    
    # Plot abundances w/ TP structure
    if (cp_dry):
        ga.plot_adiabats(atm_dry,filename=dirs["output"]+'/'+str(atm.ts)+'K_'+str(atm.ps/1e5)+'bar_dry_ga.pdf')
        atm_dry.write_PT(filename=dirs["output"]+'/'+str(atm.ts)+'K_'+str(atm.ps/1e5)+'bar_dry_pt.tsv')

    ga.plot_adiabats(atm_moist,filename=dirs["output"]+'/'+str(atm.ts)+'K_'+str(atm.ps/1e5)+'bar_moist_ga.pdf')
    atm_moist.write_PT(filename=dirs["output"]+'/'+str(atm.ts)+'K_'+str(atm.ps/1e5)+'bar_moist_pt.tsv')


    end = t.time()
    print("Runtime:", round(end - start,2), "s")

