#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:05:39 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
import numpy as np

from modules.dry_adiabat_setup import dry_adiabat_atm
from modules.moist_adjustment_H2O import moist_adj
from modules.dry_adjustment import DryAdj
from modules.simple_boundary import simple_boundary_tend
from modules.water_cloud import simple_cloud

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates
import utils.phys as phys

# Time integration for n steps
def compute_dry_adiabat(atm, dirs, standalone, calc_cf=False, rscatter=False, pure_steam_adj=False, surf_dt=False, cp_surf=1e5, mix_coeff_atmos=1e6, mix_coeff_surf=1e6, do_cloud=False):

    # Dry adiabat settings 
    rad_steps   = 100  # Maximum number of radiation steps
    conv_steps  = 30   # Number of convective adjustment steps (per radiation step)
    dT_max      = 20.  # K, Maximum temperature change per radiation step
    T_floor     = 10.  # K, Temperature floor to prevent SOCRATES crash

    # Build general adiabat structure
    atm                 = ga.general_adiabat(copy.deepcopy(atm))

    # Copy moist pressure arrays for dry adiabat
    atm_dry             = dry_adiabat_atm(atm)

    # Initialise previous OLR and TOA heating to zero
    PrevOLR_dry         = 0.
    PrevMaxHeat_dry     = 0.
    PrevTemp_dry        = atm.tmp * 0.
    # Initialize the surface temperature tendency to 0 
    dT_surf             = 0

    # Time stepping
    for i in range(0, rad_steps):

        # atm.toa_heating = 0

        # Compute radiation, midpoint method time stepping
        try:
            if do_cloud:
                atm_dry         = simple_cloud(atm_dry) # Before radiation, set up the cloud for Socrates using the current PT profile
            atm_dry         = socrates.radCompSoc(atm_dry, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud)
            dT_dry          = atm_dry.net_heating * atm_dry.dt
    
            # Limit the temperature change per step
            dT_dry          = np.where(dT_dry > dT_max, dT_max, dT_dry)
            dT_dry          = np.where(dT_dry < -dT_max, -dT_max, dT_dry)
    
            # Do the surface balance
            if surf_dt:
                # Net surface flux (positive downward, negative upward)
                net_Fs = atm_dry.SW_flux_down[len(atm_dry.tmp)] + atm_dry.LW_flux_down[len(atm_dry.tmp)] - atm_dry.LW_flux_up[len(atm_dry.tmp)] - atm_dry.SW_flux_up[len(atm_dry.tmp)]
                dT_surf += net_Fs/cp_surf 
                dT_dry, dT_surf = simple_boundary_tend(len(dT_dry)-1,atm_dry.tmp,atm_dry.ts,dT_dry,dT_surf,mix_coeff_atmos,mix_coeff_surf)
                # Apply surface heating
                atm_dry.ts      += dT_surf * atm_dry.dt

                # Update temperature of lowest cell
                k_turb          = .1
                dT_dry[-1] += (atm_dry.ts - atm_dry.tmp[-1])/mix_coeff_atmos
                dT_surf    -= (atm_dry.ts - atm_dry.tmp[-1])/mix_coeff_surf
                atm_dry.tmp[-1] += dT_dry[-1] * k_turb * (atm_dry.tmp[-1] - atm_dry.ts) 
                
            # Apply heating
            atm_dry.tmp     += dT_dry 
            
            # Pure steam convective adjustment
            if pure_steam_adj:
                dT_moist = moist_adj(atm_dry,1000.)
                atm_dry.tmp     += dT_moist
            
            # Dry convective adjustment
            for iadj in range(conv_steps):
                atm_dry     = DryAdj(atm_dry)
    
            # Temperature floor to prevent SOCRATES crash
            if np.min(atm_dry.tmp) < T_floor:
                atm_dry.tmp = np.where(atm_dry.tmp < T_floor, T_floor, atm_dry.tmp)
    
            # Convergence criteria
            dTglobal_dry    = abs(round(np.max(atm_dry.tmp-PrevTemp_dry[:]), 4))
            dTtop_dry       = abs(round(atm_dry.tmp[0]-atm_dry.tmp[1], 4))
    
            # Break criteria
            dOLR_dry        = abs(round(atm_dry.LW_flux_up[0]-PrevOLR_dry, 6))
            dbreak_dry      = (0.01*(phys.sigma*atm_dry.ts**4)**0.5)
    
            # Inform during runtime
            if i % 2 == 1 and standalone == True:
                print("Dry adjustment step %d:" % (i+1))
                print("\tOLR = " + str(atm_dry.LW_flux_up[0]) + " W/m^2", "dOLR = " + str(dOLR_dry) + " W/m^2")
                print("\tT_top = " + str(atm_dry.tmp[0]) + " K,",  "dT_top = " + str(dTtop_dry) + " K")
                print("\tT_bot = " + str(atm_dry.tmp[-1]) + " K,",  "dT_max = " + str(dTglobal_dry) + " K")
                print("\tT_surf = " +str(atm_dry.ts) + " K" )
                
    
            # Reduce timestep if heating is not converging
            if dTglobal_dry < 0.05 or dTtop_dry > dT_max:
                atm_dry.dt  = atm_dry.dt*0.99
                # if standalone == True:
                #     print("Dry adiabat not converging -> dt_new =", round(atm_dry.dt,5), "days")
    
            # Sensitivity break condition
            if (dOLR_dry < dbreak_dry) and i > 5:
                if standalone == True:
                    print("Timestepping break ->", end=" ")
                    print("dOLR/step =", dOLR_dry, "W/m^2, dTglobal_dry =", dTglobal_dry)
                break    # break here
        except:
            if standalone == True:
                print("Socrates cannot be executed properly, T profile:", atm_dry.tmp)
            break    # break here

        PrevOLR_dry       = atm_dry.LW_flux_up[0]
        PrevMaxHeat_dry   = abs(np.max(atm_dry.net_heating))
        PrevTemp_dry[:]   = atm_dry.tmp[:]

    return atm_dry
