#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 6 11:22:05 2023

@authors:    
Ryan Boukrouche (RB)
"""

import numpy as np
import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles

def simple_cloud(atm):
    """Simple water cloud routine.

    Define the three SOCRATES cloud input arrays with user-defined values where cloud condensation occurs. These values can vary with height. 
    If the general adiabat is used, the liquid water mass mixing ratio is provided by x_cond.
    The cloud top is defined by the troposphere.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

    """

    nlev = len(atm.p)
    Tmid = atm.tmp.copy()

    if np.sum(atm.x_cond['H2O']) == 0.0: # If x_cond = 0, no condensation scheme to provide lwm

        for i in range(nlev-1): 
            pp_h2o = atm.p_vol["H2O"][i]

            if (Tmid[i] <= ga.Tdew('H2O',pp_h2o) and atm.p[i] > atm.trppP):
                atm.re[i]   = atm.effective_radius
                atm.lwm[i]  = atm.liquid_water_fraction
                atm.clfr[i] = atm.cloud_fraction
            else:
                atm.re[i]   = 0.0
                atm.lwm[i]  = 0.0
                atm.clfr[i] = 0.0

    else: # a condensation scheme is providing x_cond
            
        for i in range(nlev-1): 
            pp_h2o = atm.p_vol["H2O"][i]

            if (atm.x_cond['H2O'][i] > 0.0 and atm.p[i] > atm.trppP):
                atm.re[i]   = atm.effective_radius
                atm.lwm[i]  = atm.x_cond['H2O'][i]
                atm.clfr[i] = atm.cloud_fraction
            else:
                atm.re[i]   = 0.0
                atm.lwm[i]  = 0.0
                atm.clfr[i] = 0.0
            
    return atm