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
RB: Adapted for looping over the surface temperature. 
"""

import matplotlib as mpl
mpl.use('Agg')

import time as t
import os, shutil
import numpy as np
from scipy.optimize import brentq

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.solve_pt import RadConvEqm
from modules.plot_flux_balance import plot_fluxes
from modules.relative_humidity import compute_Rh

from utils.socrates import CleanOutputDir
import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum

import utils.phys as phys

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start AEOLUS")

    start = t.time()
    ##### Settings

    # Planet 
    time = { "planet": 0., "star": 4e+9 } # yr,
    star_mass     = 1.0 #0.1*1.0                 # M_sun, mass of star
    mean_distance = 1.0 #0.0252*1.0              # au, orbital distance
    pl_radius     = 6.371e6 #1.1*6.371e6             # m, planet radius
    pl_mass       = 5.972e24 #1.05*5.972e24           # kg, planet mass

    # Boundary conditions for pressure & temperature
    T_surf        = 500.0                # K
    P_top         = 1.0                  # Pa

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
    clfr = 0.0    # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    def calculate_imbalance(index, T_surf):
    
        pH2O                    = 0.0354 * 1.0e5 #min(260e5,ga.p_sat('H2O',T_surf))       # Pa
        pCO2                    = 0.                            # Pa
        pH2                     = 0.                            # Pa
        pN2                     = 1e+5                          # Pa
        pCH4                    = 0.                            # Pa
        pO2                     = 0.                            # Pa
        pCO                     = 0.                            # Pa
        pHe                     = 0.                            # Pa
        pNH3                    = 0.                            # Pa
        P_surf                  = pH2O + pCO2 + pH2 + pN2 + pCH4 + pO2 + pCO + pHe + pNH3  # Pa

        # Define volatiles by mole fractions
        vol_partial = {}
        vol_mixing = { 
                        "H2O" : pH2O / P_surf,   
                        "CO2" : pCO2 / P_surf,
                        "H2"  : pH2  / P_surf,
                        "N2"  : pN2  / P_surf,
                        "CH4" : pCH4 / P_surf,
                        "O2"  : pO2  / P_surf,
                        "CO"  : pCO  / P_surf,
                        "He"  : pHe  / P_surf,
                        "NH3" : pNH3 / P_surf,
                    }
        
        # OR:
        # Define volatiles by partial pressures
        # P_surf = 0.0
        # vol_mixing = {}
        # vol_partial = {
        #     "H2O" : ga.p_sat('H2O',T_surf), #4e3,
        #     "NH3" : 0.0,
        #     "CO2" : 0.0,#4.0e2,
        #     "CH4" : 0.0,
        #     "CO" : 0.0,#3.0e2,
        #     "O2" : 0.0,#0.0e4,
        #     "N2" : 1.0e5,
        #     "H2" : 0.0
        #     }
    
        # Create atmosphere object
        atm = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, re, lwm, clfr,
                    vol_mixing=vol_mixing, vol_partial=vol_partial, calc_cf=calc_cf, trppT=trppT, req_levels=100, water_lookup=water_lookup, do_cloud=do_cloud)


        # Set stellar heating on or off
        if stellar_heating == False: 
            atm.toa_heating = 0.
        else:
            _, atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)
            print("TOA heating:", round(atm.toa_heating), "W/m^2")

        # Move/prepare spectral file
        print("Inserting stellar spectrum")

        StellarSpectrum.InsertStellarSpectrum(
            dirs["aeolus"]+"/spectral_files/Reach/Reach",
            dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
            dirs["output"]+"runtime_spectral_file"
        )

        # Set up atmosphere with general adiabat
        atm_dry, atm = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 

        outgoing = atm.flux_up_total[0] # LW UP + SW UP from scattering
        incoming = atm.flux_down_total[0] # LW DOWN + SW DOWN
        imbalance = incoming - outgoing

        print("T_surf, outgoing, incoming, imbalance = ", T_surf, outgoing, incoming, imbalance)

        # Plot abundances w/ TP structure
        if (cp_dry):
            ga.plot_adiabats(atm_dry,filename=dirs["output"]+"dry_ga.pdf")
            atm_dry.write_PT(filename=dirs["output"]+"dry_pt.tsv")
            plot_fluxes(atm_dry,filename=dirs["output"]+"dry_fluxes.pdf")

        ga.plot_adiabats(atm,filename=dirs["output"]+"moist_ga.pdf")
        atm.write_PT(filename=dirs["output"]+"moist_pt.tsv")
        atm.write_ncdf(dirs["output"]+"moist_atm.nc")
        plot_fluxes(atm,filename=dirs["output"]+"moist_fluxes.pdf")

        return atm #,imbalance
    
    def minimize_imbalance(bracket, tolerance=1.0):
        result = brentq(calculate_imbalance, a=bracket[0], b=bracket[1], xtol=tolerance)
        
        if result is None:
            print("Brent's method did not converge.")
        else:
            print(f"Minimized imbalance achieved! T_surf: {result}, Imbalance: {calculate_imbalance(result)}")


    loop = False
    if loop == False:

        # ------------------- Day side -------------------
        # Set up dirs
        dirs = {
                "aeolus": os.getenv('AEOLUS_DIR')+"/",
                "output": os.getenv('AEOLUS_DIR')+f"/output_day/"
                }
        
        # Tidy directory
        if os.path.exists(dirs["output"]):
            shutil.rmtree(dirs["output"])
        os.mkdir(dirs["output"])

        # Set up the dayside cloud
        re   = 1.0e-5 
        lwm  = 0.8   
        clfr = 0.0  

        Ts_day = 300.0 #ga.Tdew('H2O',sum(vol_partial.values()))
        stellar_heating = True
        atm = calculate_imbalance(0, Ts_day)
        print("Ts_day = ", Ts_day)
        OLR_day = atm.LW_flux_up[0]
        OSR_day = atm.SW_flux_up[0]
        OPR_day = atm.flux_up_total[0]
        NET_day = atm.net_flux[0]
        print("OLR_day = ", OLR_day)
        print("OSR_day = ", OSR_day)
        print("OPR_day = ", OPR_day)
        print("NET_day = ", NET_day)
        ISR = atm.SW_flux_down[0]
        print("ISR = ", ISR)

        flux_surf_down = atm.flux_down_total[-1]
        print("flux_surf_down = ", flux_surf_down)
        Ts_night = ( flux_surf_down / ((1.-atm.albedo_s)*4.*phys.sigma) )**(1./4.) 
        print("Ts_night = ", Ts_night)

        # ------------------- Night side -------------------
        dirs = {
        "aeolus": os.getenv('AEOLUS_DIR')+"/",
        "output": os.getenv('AEOLUS_DIR')+f"/output_night/"
        }
        
        if os.path.exists(dirs["output"]):
            shutil.rmtree(dirs["output"])
        os.mkdir(dirs["output"])

        # Set up the nightside cloud
        re   = 1.0e-5 
        lwm  = 0.8   
        clfr = 0.0  

        stellar_heating = False
        atm = calculate_imbalance(0, Ts_night)

        ISR_night = atm.SW_flux_down[0]
        print("ISR_night = ", ISR_night) # should be zero

        OLR_night = atm.LW_flux_up[0]
        OSR_night = atm.SW_flux_up[0]
        OPR_night = atm.flux_up_total[0]
        NET_night = atm.net_flux[0]
        print("OLR_night = ", OLR_night)
        print("OSR_night = ", OSR_night)
        print("OPR_night = ", OPR_night)
        print("NET_night = ", NET_night)

        # ------------------- Albedo from balance -------------------
        # OPR_day + OPR_night - S(1-alpha_p) = 0 ; ASR = S(1-alpha_p). S = ISR.
        ASR = OLR_day + OLR_night
        print("ASR = ", ASR)
        planetary_albedo = 1.0-(ASR/ISR)
        print("planetary_albedo = ", planetary_albedo)

    else:
        T_surf_array = np.arange(200.0, 3000.0+50.0, 50)
        #T_surf_array = np.arange(200.0, 700.0+50.0, 50)
        asr_albedo_file = []
        for i, temperature in enumerate(T_surf_array):

            # ------------------- Day side -------------------
            dirs = {
                    "aeolus": os.getenv('AEOLUS_DIR')+"/",
                    "output": os.getenv('AEOLUS_DIR')+f"/200_to_3000_K_full_condensate_retention_clear_sky/output_day_{int(temperature)}K/"
                    }
            
            # Tidy directory
            if os.path.exists(dirs["output"]):
                shutil.rmtree(dirs["output"])
            os.mkdir(dirs["output"])

            # Set up the dayside cloud
            re   = 1.0e-5 
            lwm  = 0.8   
            clfr = 0.0  

            print("index, temperature = ", i, temperature)
            stellar_heating = True
            atm = calculate_imbalance(i, temperature)
        
            OLR_day = atm.LW_flux_up[0]
            OSR_day = atm.SW_flux_up[0]
            OPR_day = atm.flux_up_total[0]
            NET_day = atm.net_flux[0]
            ISR     = atm.SW_flux_down[0]
            print("ISR = ", ISR)
            flux_surf_down = atm.flux_down_total[-1]
            Ts_night = ( flux_surf_down / ((1.-atm.albedo_s)*4.*phys.sigma) )**(1./4.) 
            print("Ts_day, Ts_night = ", temperature, Ts_night)

            # ------------------- Night side -------------------
            dirs = {
            "aeolus": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+f"/200_to_3000_K_full_condensate_retention_clear_sky/output_night_{int(temperature)}K/"
            }
            
            if os.path.exists(dirs["output"]):
                shutil.rmtree(dirs["output"])
            os.mkdir(dirs["output"])

            # Set up the nightside cloud
            re   = 1.0e-5 
            lwm  = 0.8   
            clfr = 0.0  

            stellar_heating = False
            atm = calculate_imbalance(0, Ts_night)

            OLR_night = atm.LW_flux_up[0]
            OSR_night = atm.SW_flux_up[0]
            OPR_night = atm.flux_up_total[0]
            NET_night = atm.net_flux[0]

            # ------------------- Analysis -------------------
            ASR = NET_day + NET_night
            planetary_albedo = 1.0-(ASR/ISR)

            asr_albedo_file.append([temperature, ASR, planetary_albedo])

        filename = f"asr_albedo.txt"
        with open(filename, "w") as file:
            for temperature, ASR, planetary_albedo in asr_albedo_file:
                file.write(f"{temperature}\t{ASR}\t{planetary_albedo}\n")

    # Copy this file to the output folder
    shutil.copy(__file__, os.path.join(dirs["output"], os.path.basename(__file__)))

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    end = t.time()
    print("Runtime:", round(end - start,2), "s")

