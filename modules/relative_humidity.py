#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:18:48 2024

@authors:    
Ryan Boukrouche (RB)
"""

import numpy as np
from utils.find_lcl import find_intersection
import utils.GeneralAdiabat as ga
    
def compute_Rh(atm):
    """Compute the relative humidity.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

    """

    #partial_p = atm.p_vol['H2O']-atm.x_cond['H2O']*atm.p # Remove the condensates so that only the vapor is accounted for.
    partial_p = atm.x_gas['H2O']*atm.p
    p_sat = [ga.p_sat('H2O', atm.tmp[i]) for i in range(len(atm.p))]
    R_h = partial_p/p_sat
    return R_h

def update_for_constant_RH(atm, R_h):
    """Compute the partial pressure and mixing ratio of the condensible that keeps the relative humidity constant below the LCL. Also update the total pressure.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

        R_h : float
            Relative humidity

    """

    p_sat = [ga.p_sat('H2O', atm.tmp[i]) for i in range(len(atm.p))]
    #print("update_for_constant_RH: p_sat = ", p_sat)

    i_lcl, lcl = find_intersection(p_sat, atm.p_vol['H2O'], tolerance=1e0)
    #print("update_for_constant_RH: atm.p[i_lcl:] = ", atm.p[i_lcl:])
    #print("update_for_constant_RH: atm.p_vol['H2O'][i_lcl:] = ", atm.p_vol['H2O'][i_lcl:])
    #print("update_for_constant_RH: R_h[i_lcl:] = ", R_h[i_lcl:])
    #print("update_for_constant_RH: atm.ps = ", atm.ps)
    #print("update_for_constant_RH: atm.x_gas['H2O'][i_lcl:] = ", atm.x_gas['H2O'][i_lcl:])

    #atm.p_vol['H2O'][i_lcl:] = R_h[i_lcl:] * p_sat[i_lcl:]  # Updating the partial pressure of H2O from the LCL down
    atm.p_vol['H2O'] = R_h * p_sat
    #print("update_for_constant_RH: updated atm.p_vol['H2O'][i_lcl:] = ", atm.p_vol['H2O'][i_lcl:])
    atm.p = sum(atm.p_vol[molecule] for molecule in atm.vol_list.keys()) # Updating the total pressure
    #print("update_for_constant_RH: updated atm.p[i_lcl:] = ", atm.p[i_lcl:])
    atm.ps = float(sum(atm.p_vol[molecule][-1] for molecule in atm.vol_list.keys())) # Not forgetting the surface pressure
    #print("update_for_constant_RH: updated atm.ps = ", atm.ps)
    #atm.x_gas['H2O'][i_lcl:] = R_h[i_lcl:] * p_sat[i_lcl:] / atm.p[i_lcl:] # Updating the mixing ratio of H2O
    #print("atm.p_vol['H2O'][-1], atm.x_gas['H2O'], p_sat[-1] = ", atm.p_vol['H2O'][-1], atm.x_gas['H2O'][-1], p_sat[-1])

    return atm