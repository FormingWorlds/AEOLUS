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
RB: Adapted for looping over the surface temperature and simple day-night scheme.
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
import utils.phys as phys

####################################
##### Global constants
####################################

# AU definition [m]
AU = 149597870700.
# Universal gravitational constant [m^3.kg^-1.s^-2]
G_universal = 6.67430e-11
# Universal gas constant [J.mol^-1.K^-1]
R = 8.31446261815324 
# Number of seconds in a day
day_s = 24.0*3600.0

####################################
##### Sun/Earth parameters
####################################

L_Sun        = 3.828e26    # Solar luminosity [W]
M_Sun        = 1.988500e30 # Solar mass [kg]
radius_Earth = 6.371e6     # Earth radius [m]
mass_Earth   = 5.9722e24   # Earth mass [kg]
grav_Earth   = 9.807       # Earth surface gravity [m.s^-2]

def surface_gravity(mass,radius):
	g = G_universal*np.array(mass)*mass_Earth/(np.array(radius)*radius_Earth)**2
	if g.ndim > 0:
		g = np.sort(g)
	return g

def angular_velocity(orbital_period):
	return 2.0*np.pi/(orbital_period*day_s)

# Sources: TRAPPIST-1: doi:10.3847/psj/abd022 ; Proxima-b: doi:10.3847/2041-8205/831/2/l16 (Figure 2) ; Teegarden: doi:10.1051/0004-6361/201935460 (Figure 12). 
planets = { 'Earth': {'Stellar Mass': M_Sun, 'Stellar Luminosity': L_Sun, 'Radius': radius_Earth, 'Gravity': grav_Earth, 'Semi-major axis': AU, 'Angular velocity': angular_velocity(1.0)},
            'TRAPPIST-1b': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 1.116*radius_Earth, 'Gravity': 1.102*grav_Earth, 'Semi-major axis': 0.01154*AU, 'Angular velocity': angular_velocity(1.510826)},
			'TRAPPIST-1d': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 0.788*radius_Earth, 'Gravity': 0.624*grav_Earth, 'Semi-major axis': 0.02227*AU, 'Angular velocity': angular_velocity(4.049219)},
			'TRAPPIST-1e': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 0.920*radius_Earth, 'Gravity': 0.817*grav_Earth, 'Semi-major axis': 0.02925*AU, 'Angular velocity': angular_velocity(6.101013)},
			'Proxima-b':   {'Stellar Mass': 0.1221*M_Sun, 'Stellar Luminosity': 0.001567*L_Sun, 'Radius': np.array([0.94,1.4])*radius_Earth, 'Gravity': surface_gravity(1.07,[0.94,1.4]), 'Semi-major axis': 0.04856*AU, 'Angular velocity': angular_velocity(11.1868)},
			'Teegarden-b':   {'Stellar Mass': 0.09455120*M_Sun, 'Stellar Luminosity': 0.00073*L_Sun, 'Radius': 1.02*radius_Earth, 'Mass': 1.05*mass_Earth, 'Gravity': surface_gravity(1.05,1.02), 'Semi-major axis': 0.0252*AU, 'Angular velocity': angular_velocity(4.9100)},
			'Teegarden-c':   {'Stellar Mass': 0.09455120*M_Sun, 'Stellar Luminosity': 0.00073*L_Sun, 'Radius': 1.04*radius_Earth, 'Mass': 1.11*mass_Earth, 'Gravity': surface_gravity(1.11,1.04), 'Semi-major axis': 0.0443*AU, 'Angular velocity': angular_velocity(11.409)}}

atmospheres = { 'Earth': {'H2O': 1e-3, 'CO2': 3.50e-4, 'O3': 0.07e-6, 'N2O': 0.31e-6, 'CO': 0.10e-6, 'CH4': 1.70e-6, 'O2': 0.20947, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.02e-6, 'NH3': 1.0e-7, 'HNO3': 0.0, 'N2': 0.78084, 'H2': 0.53e-6, 'He': 5.24e-6, 'OCS': 0.0},
                'Pure Steam': {'H2O': 1.0, 'CO2': 0.0, 'O3': 0.0, 'N2O': 0.0, 'CO': 0.0, 'CH4': 0.0, 'O2': 0.0, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.0, 'NH3': 0.0, 'HNO3': 0.0, 'N2': 0.0, 'H2': 0.0, 'He': 0.0, 'OCS': 0.0}}

molecules = {'Molecular Weight': {'H2O': 0.018, 'CO2': 0.044, 'O3': 0.048, 'N2O': 0.044, 'CO': 0.028, 'CH4': 0.016, 'O2': 0.032, 'NO': 0.030, 'SO2': 0.064, 'NO2': 0.046, 'NH3': 0.017, 'HNO3': 0.063, 'N2': 0.028, 'H2': 0.002, 'He': 0.004, 'OCS': 0.060},
			 'Isobaric Heat Capacity': {'H2O': 1864.0, 'CO2': 849.0, 'O3': 819.375, 'N2O': 877.364, 'CO': 1040.0, 'CH4': 2232.0, 'O2': 918.0, 'NO': 995.0, 'SO2': 624.0625, 'NO2': 805.0, 'NH3': 2175.0, 'HNO3': 849.365, 'N2': 1040.0, 'H2': 14310.0, 'He': 5197.5, 'OCS': 41.592}}

molecular_weights = list(molecules['Molecular Weight'].values())
heat_capacities   = list(molecules['Isobaric Heat Capacity'].values())

planet     = 'Teegarden-b'
atmosphere = 'Earth'

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start AEOLUS")

    start = t.time()
    ##### Settings

    # Planet 
    time = { "planet": 0., "star": 4e+9 } # yr,
    star_mass     = max(0.1,planets[planet]['Stellar Mass']/M_Sun)*1.0 # M_sun, mass of star. Minimum supported is 0.1.
    mean_distance = planets[planet]['Semi-major axis']/AU * 1.0        # au, orbital distance
    pl_radius     = planets[planet]['Radius']                          # m, planet radius
    pl_mass       = planets[planet]['Mass']                            # kg, planet mass

    # Boundary conditions for pressure
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
    do_cloud = True
    # Options activated by do_cloud
    re   = 1.0e-5 # Effective radius of the droplets [m] (drizzle forms above 20 microns)
    lwm  = 0.8    # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
    clfr = 0.8    # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Set up dirs
    if os.environ.get('AEOLUS_DIR') == None:
        raise Exception("Environment variables not set! Have you sourced AEOLUS.env?")

    def calculate_imbalance(T_surf, index=0):

        pH2O                    = min(260e5,ga.p_sat('H2O',T_surf))       # Pa
        pCO2                    =      150.0 # Pa
        pH2                     = 0.                            # Pa
        pN2                     = 1e+5                          # Pa
        pCH4                    =      0.96673 # Pa
        pO2                     =      0.25e5 # Pa
        pO3                     = 0.0014                             # Pa   
        pCO                     = 0. #10.                            # Pa
        pHe                     = 0. #10.                            # Pa
        pNH3                    = 0. #10.                            # Pa
        P_surf                  = pH2O + pCO2 + pH2 + pN2 + pCH4 + pO2 + pO3 + pCO + pHe + pNH3  # Pa

        # Define volatiles by mole fractions
        vol_partial = {}
        vol_mixing = { 
                        "H2O" : pH2O / P_surf,   
                        "CO2" : pCO2 / P_surf,
                        "H2"  : pH2  / P_surf,
                        "N2"  : pN2  / P_surf,
                        "CH4" : pCH4 / P_surf,
                        "O2"  : pO2  / P_surf,
                        "O3"  : pO3  / P_surf,
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
            atm.instellation = 0.
            atm.alpha_cloud  = 0.0 # Dry nightside
        else:
            atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)
            print("Instellation:", round(atm.instellation), "W/m^2")
            atm.alpha_cloud  = 0.9 # Moist dayside

        # Move/prepare spectral file
        print("Inserting stellar spectrum")

        StellarSpectrum.InsertStellarSpectrum(
            dirs["aeolus"]+"/spectral_files/Reach/Reach_cloud/Reach",
            #dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
            dirs["aeolus"]+"/rad_trans/socrates_code/data/solar/trappist1",
            dirs["output"]+"runtime_spectral_file"
        )

        # Set up atmosphere with general adiabat
        atm_dry, atm = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 

        outgoing = atm.flux_up_total[0]   # LW UP + SW UP from scattering
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

        return atm
    
    loop = True
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
        re   = 1.0e-5 # Effective radius of the droplets [m] (drizzle forms above 20 microns)
        lwm  = 0.8    # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
        clfr = 0.4    # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

        Ts_day = 300.0 #ga.Tdew('H2O',sum(vol_partial.values()))
        stellar_heating = True
        atm = calculate_imbalance(Ts_day)
        print("atm.alpha_cloud, day = ", atm.alpha_cloud)
        print("atm.tmp day = ", atm.tmp)
        print("atm.p day = ", atm.p)
        print("Dayside tropopause pressure = ", atm.trppP)

        print("Ts_day = ", Ts_day)
        OLR_day = atm.LW_flux_up[0]
        OSR_day = atm.SW_flux_up[0]
        OPR_day = atm.flux_up_total[0]
        NET_day = atm.net_flux[0]
        ASR_day = atm.SW_flux_down[0] - atm.SW_flux_up[0]
        print("OLR_day = ", OLR_day)
        print("OSR_day = ", OSR_day)
        print("OPR_day = ", OPR_day)
        print("NET_day = ", NET_day)
        print("ASR_day = ", ASR_day)
        ISR = atm.SW_flux_down[0]
        print("ISR = ", ISR)

        flux_surf_down = atm.LW_flux_down[-1]
        print("LW_flux_down[-1] = ", flux_surf_down)
        print("flux_down_total[-1] = ", atm.flux_down_total[-1])
        print("SW_flux_down = ", atm.SW_flux_down)
        Ts_night = ( flux_surf_down / ((1.-atm.albedo_s)*phys.sigma) )**(1./4.) 
        if Ts_night >= Ts_day: # if the nightside surface is warmer than the dayside surface, assume that efficient nightside IR escape to space will bring it down
            Ts_night = Ts_day
        np.savetxt('ts_night.txt', [Ts_night])
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
        atm = calculate_imbalance(Ts_night)
        print("atm.alpha_cloud, night = ", atm.alpha_cloud)
        print("atm.tmp night = ", atm.tmp)
        print("atm.p night = ", atm.p)
        print("Nightside tropopause pressure = ", atm.trppP)

        ISR_night = atm.SW_flux_down[0]
        print("ISR_night = ", ISR_night) # should be zero

        OLR_night = atm.LW_flux_up[0]
        OSR_night = atm.SW_flux_up[0]
        OPR_night = atm.flux_up_total[0]
        NET_night = atm.net_flux[0]
        ASR_night = atm.SW_flux_down[0] - atm.SW_flux_up[0]
        print("OLR_night = ", OLR_night)
        print("OSR_night = ", OSR_night)
        print("OPR_night = ", OPR_night)
        print("NET_night = ", NET_night)
        print("ASR_night = ", ASR_night)

    else:
        T_surf_array = np.concatenate((np.arange(200.0,290.0+10,10),np.arange(300.0,3000.0+50.0,50))) #np.arange(200.0, 3000.0+50.0, 50)

        for i, temperature in enumerate(T_surf_array):
             
            # ------------------- Day side -------------------
            # Set up dirs
            dirs = {
                    "aeolus": os.getenv('AEOLUS_DIR')+"/",
                    "output": os.getenv('AEOLUS_DIR')+f"/200_to_3000_K_full_condensate_retention_cloudy_40/output_day_{int(temperature)}K/"
                    }
            
            # Tidy directory
            if os.path.exists(dirs["output"]):
                shutil.rmtree(dirs["output"])
            os.makedirs(dirs["output"])

            # Set up the dayside cloud
            re   = 1.0e-5 # Effective radius of the droplets [m] (drizzle forms above 20 microns)
            lwm  = 0.8    # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
            clfr = 0.4    # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

            stellar_heating = True
            atm = calculate_imbalance(temperature)

            print("temperature = ", temperature)
            OLR_day = atm.LW_flux_up[0]
            OSR_day = atm.SW_flux_up[0]
            OPR_day = atm.flux_up_total[0]
            NET_day = atm.net_flux[0]
            ASR_day = atm.SW_flux_down[0] - atm.SW_flux_up[0]
            print("OLR_day = ", OLR_day)
            print("OSR_day = ", OSR_day)
            print("OPR_day = ", OPR_day)
            print("NET_day = ", NET_day)
            print("ASR_day = ", ASR_day)
            ISR = atm.SW_flux_down[0]
            print("ISR = ", ISR)

            flux_surf_down = atm.LW_flux_down[-1]
            print("flux_surf_down = ", flux_surf_down)
            Ts_night = ( flux_surf_down / ((1.-atm.albedo_s)*phys.sigma) )**(1./4.) 
            if Ts_night >= temperature: # if the nightside surface is warmer than the dayside surface, assume that efficient nightside IR escape to space will bring it down
                Ts_night = temperature
            np.savetxt('ts_night.txt', [Ts_night])
            print("Ts_night = ", Ts_night)

            # ------------------- Night side -------------------
            dirs = {
            "aeolus": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+f"/200_to_3000_K_full_condensate_retention_cloudy_40/output_night_{int(temperature)}K/"
            }
            
            if os.path.exists(dirs["output"]):
                shutil.rmtree(dirs["output"])
            os.makedirs(dirs["output"])

            # Set up the nightside cloud
            re   = 1.0e-5 
            lwm  = 0.8   
            clfr = 0.0  

            stellar_heating = False
            atm = calculate_imbalance(Ts_night)

            ISR_night = atm.SW_flux_down[0]
            print("ISR_night = ", ISR_night) # should be zero

            OLR_night = atm.LW_flux_up[0]
            OSR_night = atm.SW_flux_up[0]
            OPR_night = atm.flux_up_total[0]
            NET_night = atm.net_flux[0]
            ASR_night = atm.SW_flux_down[0] - atm.SW_flux_up[0]
            print("OLR_night = ", OLR_night)
            print("OSR_night = ", OSR_night)
            print("OPR_night = ", OPR_night)
            print("NET_night = ", NET_night)
            print("ASR_night = ", ASR_night)



    # Copy this file to the output folder
    shutil.copy(__file__, os.path.join("output", os.path.basename(__file__)))

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    end = t.time()
    print("Runtime:", round(end - start,2), "s")