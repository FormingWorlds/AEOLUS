#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:01:27 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
import numpy as np
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere
from modules.water_cloud import simple_cloud
from modules.relative_humidity import compute_Rh, update_for_constant_RH
from modules.make_inversion import inversion

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates

#rh_computed = False # for the relative humidity

def compute_moist_adiabat(atm, dirs, standalone, trppD, calc_cf=False, rscatter=False, do_cloud=False):
    """Compute moist adiabat case

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py
        dirs : dict
            Named directories
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        calc_cf : bool
            Calculate contribution function?
        rscatter : bool
            Include Rayleigh scattering?
            
    """

    atm_moist = copy.deepcopy(atm)
    
    atm_moist = ga.general_adiabat(atm_moist) # Build initial general adiabat structure

    if atm_moist.instellation == 0.: # create an inversion on the nightside
        # Read dayside arrays
        data = np.loadtxt('dayside.txt')
        data_l = np.loadtxt('dayside_l.txt')
        tmp_day  = data[:, 0] 
        tmpl_day = data_l[:, 0] 
        p_day    = data[:, 1]  
        pl_day   = data_l[:, 1]
        ts_night = np.loadtxt('ts_night.txt')
        atm_moist.tmp = inversion(tmp_day, atm_moist.ts, ts_night, p_day, atm_moist.p, True)
        atm_moist.tmpl = inversion(tmpl_day, atm_moist.ts, ts_night, pl_day, atm_moist.pl, True)
    else:
        np.savetxt('dayside.txt', np.column_stack((atm_moist.tmp, atm_moist.p)))  # save the dayside TP to use on the nightside.
        np.savetxt('dayside_l.txt', np.column_stack((atm_moist.tmpl, atm_moist.pl)))
    
    atm_moist.rh = compute_Rh(atm_moist)

    if do_cloud:
        atm_moist = simple_cloud(atm_moist) # Before radiation, set up the cloud for Socrates using the current PT profile

    # Run SOCRATES
    atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    # Calculate tropopause
    if (trppD == True) or (atm_moist.trppT > atm_moist.minT):
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist,trppD, verbose=standalone)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        if do_cloud:
            atm_moist = simple_cloud(atm_moist) # Update cloud location after previous PT changes
        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist
