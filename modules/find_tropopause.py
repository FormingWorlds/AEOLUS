#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:59:40 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)
"""

import numpy as np

def find_tropopause(atm_moist, dynamic: bool, verbose=True):
    """Computes tropopause location via two methods: dynamically based on heating rates, or by a set temperature value.

    Parameters
    ----------
        atm_moist : atmos
            Atmosphere (moist) object from atmosphere_column.py
        dynamic : bool
            Calculate dynamically with heating rates? If not, use fixed value of trppT.
            
    """

    # Heating criterion
    if dynamic:

        if verbose:
            print("TROPOPAUSE SET BY HEATING")

        # Find tropopause index
        trpp_idx   = 0 
        signchange = ((np.roll(np.sign(atm_moist.net_heating), 1) - np.sign(atm_moist.net_heating)) != 0).astype(int)[1:]
        signchange_indices = np.nonzero(signchange)[0]

        # Criteria for "significant " heating
        DeltaT_max_sign  = 50.
        DeltaT_at_trpp   = 30.
        DeltaT_mean_sign = 10.
        dZ_strato        = round(len(atm_moist.net_heating)*0.02)
        
        # If heating sign change below TOA -> tropopause
        if np.size(signchange_indices) > 0 and np.max(atm_moist.net_heating) > DeltaT_max_sign:

            # First guess: maximum heating index
            # print(np.argmax(atm_moist.net_heating))
            max_heat_idx = np.argmax(atm_moist.net_heating)

            # Lower height while heating still significant
            while atm_moist.net_heating[max_heat_idx] > DeltaT_at_trpp and max_heat_idx < len(atm_moist.p)-1:
                max_heat_idx += 1

            trpp_idx     = max_heat_idx


            # # First guess: uppermost sign change (below TOA)
            # if atm_moist.net_heating[signchange_indices[0]-1] > 0 and np.mean(atm_moist.net_heating[:signchange_indices[0]-1]) > DeltaT_mean_sign:
            #     trpp_idx = signchange_indices[0]

            # # Decrease trpp height (== increase idx) while heating in trpp layer is significant
            # for idx, sgn_idx_top in enumerate(signchange_indices):

            #     if idx < np.size(signchange_indices)-1:
            #         sgn_idx_down = signchange_indices[idx+1]
                
            #         # Check mean and max heating and extent of layer above trpp idx
            #         if np.mean(atm_moist.net_heating[sgn_idx_top:sgn_idx_down]) > DeltaT_mean_sign and np.max(atm_moist.net_heating[sgn_idx_top:sgn_idx_down]) > DeltaT_max_sign and abs(sgn_idx_down-sgn_idx_top) > dZ_strato:
            #             trpp_idx = sgn_idx_down

            # Only consider tropopause if deeper than X% below TOA
            if trpp_idx < dZ_strato: 
                trpp_idx = 0

            # Only consider if mean heating above tropopause significant
            if np.mean(atm_moist.net_heating[0:trpp_idx]) < DeltaT_mean_sign:
                trpp_idx = 0


        # # If heating everywhere (close to star) & heating is significant
        # if np.size(signchange_indices) <= 1 and np.mean(atm_moist.net_heating) > DeltaT_mean_sign:
        #     trpp_idx = np.size(atm_moist.tmp)-1

        # If significant tropopause found or isothermal atmosphere from stellar heating
        if trpp_idx != 0:

            # # Print tropopause index for debugging
            # print("Tropopause @ (index, P/Pa, T/K):", trpp_idx, round(atm_moist.pl[trpp_idx],3), round(atm_moist.tmpl[trpp_idx],3))
        
            atm_moist.trppidx   = trpp_idx                  # index
            atm_moist.trppP     = atm_moist.pl[trpp_idx]    # pressure 
            atm_moist.trppT     = atm_moist.tmpl[trpp_idx]  # temperature

    # Temperature criterion
    else:

        if verbose:
            print("TROPOPAUSE SET BY CONTANT VALUE OF", atm_moist.trppT, "K")

        # Tropopause is triggered
        if np.any(atm_moist.tmpl <= atm_moist.trppT) or np.any(atm_moist.tmp <= atm_moist.trppT):
            # Loop from bottom upwards to find first level where temperature drops below required value
            for i in range(len(atm_moist.tmpl)-1,-1,-1):
                if atm_moist.tmpl[i] < atm_moist.trppT:
                    # Found tropopause, so break loop and set result
                    atm_moist.trppidx = i
                    atm_moist.trppP   = atm_moist.pl[i] 
                    break

            # Make-safe the edge case where trppT > T_surf
            atm_moist.trppidx = int(min(atm_moist.trppidx, len(atm_moist.tmpl)-1))
            atm_moist.trppidx = int(min(atm_moist.trppidx, len(atm_moist.tmp)-1))
            atm_moist.trppP   = atm_moist.pl[atm_moist.trppidx] 

        # Tropopause isn't triggered
        else:
            atm_moist.trppidx   = -1
            atm_moist.trppP     = 1.0e-10

    return atm_moist
