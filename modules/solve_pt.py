#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:05:00 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)
"""

import pickle as pkl
import scipy.optimize as optimise

from modules.compute_moist_adiabat import compute_moist_adiabat
from modules.dry_adiabat_timestep import compute_dry_adiabat
from utils.atmosphere_column import atmos
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere
import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates

def RadConvEqm(dirs, time, atm, standalone:bool, cp_dry:bool, trppD:bool, calc_cf:bool, rscatter:bool, do_cloud:bool, 
               pure_steam_adj=False, surf_dt=False, cp_surf=1e5, mix_coeff_atmos=1e6, mix_coeff_surf=1e6):
    """Sets the atmosphere to a temperature profile using the general adiabat. 
    
    Optionally does radiative time-stepping, but this is deprecated.

    Parameters
    ----------
        dirs : dict
            Named directories
        time : dict
            Dictionary of time values, including stellar age and evolution of planet
        atm : atmos
            Atmosphere object from atmosphere_column.py
        standalone : bool
            Running AEOLUS as standalone code?
        cp_dry : bool
            Compute dry adiabat case
        trppD : bool 
            Calculate tropopause dynamically?
        calc_cf : bool
            Calculate contribution function?
        rscatter : bool
            Include Rayleigh scattering
        do_cloud : bool
            Include cloud radiation
        pure_steam_adj : bool
            Use pure steam adjustment?
        surf_dt : float
            Timestep to use for T_surf timestepping cases
        cp_surf : float
            Surface heat capacity in T_surf timestepping cases
        mix_coeff_atmos : float
            Mixing coefficient (atmosphere) for T_surf timestepping cases?
        mix_coeff_surf : float
            Mixing coefficient (surface) for T_surf timestepping cases?
            
    """

    ### Moist/general adiabat

    atm_moist = compute_moist_adiabat(atm, dirs, standalone, trppD, calc_cf, rscatter, do_cloud)

    ### Dry adiabat
    if cp_dry == True:

        # Compute dry adiabat  w/ timestepping
        atm_dry   = compute_dry_adiabat(atm, dirs, standalone, calc_cf, rscatter, pure_steam_adj, surf_dt, cp_surf, mix_coeff_atmos, mix_coeff_surf, do_cloud)

        if standalone == True:
            print("Net, OLR => moist:", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print("| dry:", str(round(atm_dry.net_flux[0], 3)), str(round(atm_dry.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print()
    else: 
        atm_dry = {}
    
    # Plot
    if standalone == True:
        #plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs)
        # Save to disk
        with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
            pkl.dump(atm_moist, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    return atm_dry, atm_moist


def MCPA(dirs, atm, standalone:bool, trppD:bool, rscatter:bool):
    """Calculates the temperature profile using the multiple-condensible pseudoadiabat.

    Prescribes a stratosphere, and also calculates fluxes.

    Parameters
    ----------
        dirs : dict
            Named directories
        atm : atmos
            Atmosphere object from atmosphere_column.py
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include rayleigh scattering?
            
    """

    ### Moist/general adiabat
    return compute_moist_adiabat(atm, dirs, standalone, trppD, False, rscatter)

def MCPA_CBL(dirs, atm_inp, trppD:bool, rscatter:bool, atm_bc:int=0, T_surf_guess:float=-1, T_surf_max:float=-1, method:int=0):
    """Calculates the temperature profile using the multiple-condensible pseudoadiabat and steps T_surf to conserve energy.

    Prescribes a stratosphere, and also calculates fluxes. Finds T_surf and fluxes such that the conductive BL conduction
    equation is satisfied (Fourier's law).

    Parameters
    ----------
        dirs : dict
            Named directories
        atm_inp : atmos
            Atmosphere object from atmosphere_column.py
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include rayleigh scattering?
        
        atm_bc : int
            Where to measure boundary condition flux (0: TOA, 1: Surface).
        T_surf_guess : float
            Surface temperature guess (-1 to disable)
        T_surf_max : float
            Surface temperature ceiling (-1 to disable)
        method : int
            Root finding method (0: secant, 1: brentq)
    """

    # Store constants
    #    Passed into atmos() constructor ...
    trppT = atm_inp.trppT
    minT = atm_inp.minT; maxT = atm_inp.maxT
    Psurf = atm_inp.ps; ptop = atm_inp.ptop
    pl_r = atm_inp.planet_radius; pl_m = atm_inp.planet_mass
    vol_list = atm_inp.vol_list
    nlev_save = atm_inp.nlev_save
    #    Passed later ...
    attrs = {}
    for a in ["alpha_cloud", "instellation", "zenith_angle", "albedo_pl", 
                "inst_sf", "skin_k", "skin_d", "tmp_magma"]:
        attrs[a] = getattr(atm_inp,a)

    # Calculate conductive flux for a given atmos object 'a'
    def skin(a):
        return a.skin_k / a.skin_d * (a.tmp_magma - a.ts)
    
    # Initialise a new atmos object
    def ini_atm(Ts):
        _atm = atmos(Ts, Psurf, ptop, 
                   pl_r, pl_m, 
                   vol_mixing=vol_list, trppT=trppT, 
                   minT=minT, maxT=maxT, 
                   req_levels=nlev_save)
        for a in attrs.keys():
            setattr(_atm,a,attrs[a])
        return _atm
    
    # We want to optimise this function (returns residual of F_atm and F_skn, given T_surf)
    def func(x):

        print("Evaluating at T_surf = %.1f K" % x)
        atm_tmp = ga.general_adiabat(ini_atm(x))

        # Calculate tropopause
        if (trppD == True) or (atm_tmp.trppT > atm_tmp.minT):

            if trppD:
                atm_tmp = socrates.radCompSoc(atm_tmp, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
        
            # Find tropopause index
            atm_tmp = find_tropopause(atm_tmp,trppD, verbose=False)

            # Reset stratosphere temperature and abundance levels
            atm_tmp = set_stratosphere(atm_tmp)

        # Calculate final fluxes from T(p)
        atm_tmp = socrates.radCompSoc(atm_tmp, dirs, recalc=True, calc_cf=False, rscatter=rscatter)

        if atm_bc == 0:
            F_atm = atm_tmp.net_flux[0]  
        else:
            F_atm = atm_tmp.net_flux[-1]  

        F_skn = skin(atm_tmp)
        print("    F_atm = %+.2e W m-2      F_skn = %+.2e W m-2" % (F_atm,F_skn))

        del atm_tmp
        return float(F_skn - F_atm)
    
    print("Solving for global energy balance with conductive lid (T_magma = %.1f K)" % attrs["tmp_magma"])

    # Use an 'initial guess' method
    if method == 0:
        x0 = atm_inp.ts          # guess 1
        if T_surf_guess < minT:  # guess 2 (if not disabled)
            x1 = attrs["tmp_magma"] * 0.8
        else:
            x1 = T_surf_guess
        r = optimise.root_scalar(func, method='secant', x0=x0, x1=x1, xtol=1e-5, maxiter=20)

    # Use a 'bracketing' method
    elif method == 1:
        bracket = [800.0, maxT]
        if T_surf_max > minT:
            bracket = [bracket[0], min(T_surf_max, maxT)]
        r = optimise.root_scalar(func, method='brentq', bracket=bracket, xtol=1e-2, maxiter=20)

    else:
        raise Exception("Invalid solution method chosen (%d)" % method)

    # Extract solution
    T_surf_sol = float(r.root)
    succ  = bool(r.converged)

    if not succ:
        print("WARNING: Did not find solution for surface skin balance")
    else:
        print("Found surface solution")
        
    # Check bounds on T_surf
    T_surf = max(T_surf_sol, minT)
    T_surf = min(T_surf,     maxT)
    if T_surf_max > minT:
        T_surf = min(T_surf, T_surf_max)

    if T_surf != T_surf_sol:
        print("T_surf limits activated")
        print("    Found T_surf = %g K" % T_surf_sol)
        print("    Using T_surf = %g K" % T_surf)

    # Get atmosphere state from solution value
    atm = compute_moist_adiabat(ini_atm(T_surf), dirs, False, trppD, False, rscatter)
    atm.ts = T_surf

    if atm_bc == 0:
        F_atm = atm.net_flux[0]  
    else:
        F_atm = atm.net_flux[-1]  

    F_olr = atm.LW_flux_up[0]

    print("    T_surf = %g K"       % T_surf)
    print("    F_atm  = %.4e W m-2" % F_atm)
    print("    F_skn  = %.4e W m-2" % skin(atm))
    print("    F_olr  = %.4e W m-2" % F_olr)

    return atm

