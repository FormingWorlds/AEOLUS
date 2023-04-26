#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:01:27 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.SocRadModel as SocRadModel

def compute_moist_adiabat(atm, dirs, standalone, trpp, calc_cf=False, rscatter=False):

    atm_moist = copy.deepcopy(atm)

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(atm_moist)

    # Run SOCRATES
    atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    if trpp == True:
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=calc_cf, rscatter=rscatter)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist
