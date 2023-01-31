#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:20:27 2023

@authors: 
Mark Hammond (MH)
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import sys
if (sys.path[-1]!='./modules'):
    sys.path.append('./modules')
    
import time as t
from atmosphere_column import atmos
from stellar_luminosity import InterpolateStellarLuminosity
from radcoupler import RadConvEqm
import os
import numpy as np

try:
    #import phys
    import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
    import SocRadModel
    from atmosphere_column import atmos
except:
    #import atm_rad_conv.phys as phys
    import atm_rad_conv.GeneralAdiabat as ga
    import atm_rad_conv.SocRadModel as SocRadModel
    from atm_rad_conv.atmosphere_column import atmos

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    start = t.time()
    ##### Settings

    # Constants
    L_sun                   = 3.828e+26        # W, IAU definition
    AU                      = 1.495978707e+11  # m
    
    # Planet age and orbit
    time = { "planet": 0., "star": 4567e+6 } # yr,
    # time_current  = 0                 # yr, time after start of MO
    # time_offset   = 4567e+6           # yr, time relative to star formation
    star_mass     = 0.08    #1.0                 # M_sun, mass of star
    mean_distance = 0.01154 #1.0                 # au, orbital distance

    # Surface pressure & temperature
    
    T_surf        = 1225.0                # K

    # # Volatile molar concentrations: must sum to ~1 !
    # P_surf        = 210e+5              # Pa
    # vol_list = { 
    #               "H2O"  : 100e5/P_surf,
    #               "CO2"  : 100e5/P_surf,
    #               "H2"   : 0., 
    #               "NH3"  : 100e5/P_surf,
    #               "N2"   : 10e5/P_surf,  
    #               "CH4"  : 0., 
    #               "O2"   : 0., 
    #               "CO"   : 0., 
    #               # # No thermodynamic data, RT only
    #               # "O3"   : 0.01, 
    #               # "N2O"  : 0.01, 
    #               # "NO"   : 0.01, 
    #               # "SO2"  : 0.01, 
    #               # "NO2"  : 0.01, 
    #               # "HNO3" : 0.01, 
    #               # "He"   : 0.01, 
    #               # "OCS"  : 0.01,
    #             }
    
    # Partial pressure guesses
    #P_surf      = "calc" # does not do anything  
     # Volatiles considered
    vol_list    = { 
                          "H2O" :  1.00e+6,
                          "NH3" :  0.,
                          "CO2" :  0.,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  0.,
                          "H2"  :  0.
                        }

    P_surf      = sum(vol_list.values())

    # Stellar heating on/off
    stellar_heating = True
    
    # False: interpolate luminosity from age andm mass tables. True: define a custom instellation.
    custom_ISR = True

    # Rayleigh scattering on/off
    rscatter = True

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = True

    # Surface temperature time-stepping
    surf_dt = True
    # Options activated by surf_dt
    cp_surf = 1e5         # Heat capacity of the ground [J.kg^-1.K^-1]
    mix_coeff_atmos = 1e6 # mixing coefficient of the atmosphere [s]
    mix_coeff_surf  = 1e6 # mixing coefficient at the surface [s]

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, vol_list, calc_cf=calc_cf)

    # Compute stellar heating
    if custom_ISR: # If the age is not known with enough accuracy, define the ISR manually
        L_star          = 0.000553*L_sun                                           # Trappist-1 luminosity
        atm.toa_heating = ( 1. - 0. ) * (L_star/(4.*np.pi*(mean_distance*AU)**2.)) # Tidally-locked, zero albedo
    else:
        atm.toa_heating = InterpolateStellarLuminosity(L_sun, AU, star_mass, time, mean_distance, atm.albedo_pl, Sfrac)
    
    
    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.toa_heating = 0.
    else:
        print("TOA heating:", round(atm.toa_heating), "W/m^2")
        
    # Compute heat flux
    atm_dry, atm_moist = RadConvEqm({"output": os.getcwd()+"/output", "rad_conv": os.getcwd()}, time, atm, [], [], standalone=True, cp_dry=True, trpp=False, calc_cf=calc_cf, rscatter=rscatter, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 
    
    # Plot abundances w/ TP structure
    #ga.plot_adiabats(atm_moist)

    end = t.time()
    print(end - start)

    atm_print = atm_dry # atm_moist
    print("TOA heating - atm_print:", round(atm_print.toa_heating), "W/m^2")
    print("atm_print.flux_down_total[0] = ", atm_print.flux_down_total[0], "W/m^2")
    ClausiusClapeyron = [ga.Tdew('H2O',p) for p in atm_print.p]

    print("OPR = ", atm_print.flux_up_total[0], "ASR = ", atm_print.flux_down_total[0], "ASR-OPR = ", atm_print.flux_down_total[0]-atm_print.flux_up_total[0])

    np.savetxt(f'data_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.pl,atm_print.tmpl,atm_print.SW_flux_up,atm_print.LW_flux_up,atm_print.SW_flux_down,atm_print.LW_flux_down,atm_print.SW_flux_net,atm_print.LW_flux_net,atm_print.flux_up_total,atm_print.flux_down_total,atm_print.net_flux)))

    np.savetxt(f'data_heating_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.SW_heating,atm_print.LW_heating,atm_print.net_heating)))

    np.savetxt(f'data_spectral_TOA_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.LW_spectral_flux_up[:,0],atm_print.SW_spectral_flux_up[:,0],atm_print.net_spectral_flux[:,0])))

    np.savetxt(f'data_spectral_Surface_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.LW_spectral_flux_up[:,-1],atm_print.SW_spectral_flux_up[:,-1],atm_print.net_spectral_flux[:,-1])))
