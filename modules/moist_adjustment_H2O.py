#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:10:01 2023

@authors:    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)
"""

import numpy as np
import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles

def moist_adj(atm, conv_timescale, nb_convsteps = 10):
    """Moist adjustment routine, as formulated by RB.

    Assumes constant latent heat and does not conserve moist enthalpy. 

    Parameters
    ----------
        atm : atmos
            Atmosphere object
        conv_timescale : float
            Convective timescale

        nb_convsteps : int
            Number of passes

    Returns
    ----------
        dT_conv : np.ndarray
            Array of temperature changes introduced by adjustment

    """

    nlev = len(atm.p)
    
    dT_conv = np.zeros(nlev)
    Tmid_cc = atm.tmp.copy()

    if not ("H2O" in atm.vol_list.keys()):
        return dT_conv
        
    for _ in range(nb_convsteps):
        did_adj = False
        #------------------- L = cst -------------------
        for i in range(nlev-1): #Downward pass
            pp_h2o = atm.p_vol["H2O"][i]
            if (pp_h2o > 1e-10):
                continue
            if (Tmid_cc[i] < ga.Tdew('H2O',pp_h2o)):
                Tmid_cc[i]=ga.Tdew('H2O',pp_h2o)
                did_adj = True

        for i in range(-2,1,-1): #Upward pass
            pp_h2o = atm.p_vol["H2O"][i]
            if (pp_h2o > 1e-10):
                continue
            if (Tmid_cc[i] < ga.Tdew('H2O',pp_h2o)):
                Tmid_cc[i]=ga.Tdew('H2O',pp_h2o)
                did_adj = True

        # If no adjustment required, exit the loop
        if (did_adj == False):
            break

    # Change in temperature is Tmid_cc - Tmid
    dT_conv[:] = (Tmid_cc[:] - atm.tmp[:])/conv_timescale
    return dT_conv
