#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 11:57:03 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
    
def set_stratosphere(atm):

    trpp_idx = int(atm.trppidx)
    trpp_prs = atm.trppP
    trpp_tmp = atm.trppT

    # Standard nodes
    for prs_idx, prs in enumerate(atm.p):
        if prs < trpp_prs:
            atm.tmp[prs_idx] = trpp_tmp

    # Staggered nodes
    for prsl_idx, prls in enumerate(atm.pl):
        if prls < trpp_prs:
            atm.tmpl[prsl_idx] = trpp_tmp

    # Set mixing ratios to same as tropopause
    for idx in reversed(range(0, trpp_idx)):
    
        atm.cp[idx] = 0.

        # Volatile abundances
        for vol in atm.vol_list.keys():
        
            # Saturation vapor pressure
            p_vol_sat     = ga.p_sat(vol, atm.tmp[idx])

            # If still condensible
            if atm.p[idx] > p_vol_sat:

                cond_diff            = (atm.x_cond[vol][idx] - atm.x_cond[vol][trpp_idx])

                atm.xc[idx]          -= cond_diff
                atm.xv[idx]          += cond_diff

                atm.x_cond[vol][idx] = atm.x_cond[vol][trpp_idx]
                atm.x_gas[vol][idx]  = atm.x_gas[vol][trpp_idx]
                atm.p_vol[vol][idx]  = atm.x_gas[vol][idx] * atm.p[idx]

            # If not anymore
            else:
                atm.xc[idx]          -= atm.x_cond[vol][idx]
                atm.x_gas[vol][idx]  = atm.x_gas[vol][trpp_idx]
                atm.xd[idx]          += atm.x_gas[vol][idx] + atm.x_cond[vol][idx]
                atm.xv[idx]          -= atm.x_gas[vol][idx]

                atm.x_cond[vol][idx] = 0.
                
                atm.p_vol[vol][idx]  = atm.x_gas[vol][trpp_idx] * atm.p[idx]
                

            # Renormalize cp w/ molar concentration
            atm.cp[idx]   += (atm.x_gas[vol][idx] + atm.x_cond[vol][idx]) * ga.cpv(vol, atm.tmp[idx]) / (atm.xd[idx] + atm.xv[idx] + atm.xc[idx]) # w/ cond

    return atm
