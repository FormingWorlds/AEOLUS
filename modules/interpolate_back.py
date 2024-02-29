#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 17:14:08 2024

@authors:    
Ryan Boukrouche (RB)
"""

import numpy as np
from scipy.interpolate import CubicSpline

def interpolate_back_to_10000(atm):
    """Interpolate the General adiabat arrays of length 100 back to 10000 to input them in the general adiabat again.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

    """

    # Interpolate the arrays to length 10000. A cublic spline is noticably more accurate than a linear interpolation here.
    pressure_interp  = np.logspace(np.log10(np.min(atm.p)), np.log10(np.max(atm.p)), 10000) # array of indices corresponding to the array of length 10000

    x = np.log(atm.p) #if np.all(atm.p[1:] > atm.p[:-1]) else np.flip(np.log(atm.p))

    y = atm.rh 
    cs = CubicSpline(x, y) # CubicSpline requires x to be strictly increasing
    atm.rh = np.flip(cs(np.log(pressure_interp)))

    y = atm.xd 
    cs = CubicSpline(x, y) 
    atm.xd = np.flip(cs(np.log(pressure_interp)))

    y = atm.tmp 
    cs = CubicSpline(x, y) 
    atm.tmp = np.flip(cs(np.log(pressure_interp)))

    y = atm.z 
    cs = CubicSpline(x, y) 
    atm.z = np.flip(cs(np.log(pressure_interp)))

    y = atm.grav_z
    cs = CubicSpline(x, y) 
    atm.grav_z = np.flip(cs(np.log(pressure_interp)))

    y = atm.mu 
    cs = CubicSpline(x, y) 
    atm.mu = np.flip(cs(np.log(pressure_interp)))

    y = atm.xc 
    cs = CubicSpline(x, y) 
    atm.xc = np.flip(cs(np.log(pressure_interp)))

    y = atm.xv 
    cs = CubicSpline(x, y) 
    atm.xv = np.flip(cs(np.log(pressure_interp)))

    y = atm.cp 
    cs = CubicSpline(x, y) 
    atm.cp = np.flip(cs(np.log(pressure_interp)))

    y = atm.rho 
    cs = CubicSpline(x, y) 
    atm.rho = np.flip(cs(np.log(pressure_interp)))

    atm.pl = np.logspace(np.log10(np.min(atm.p)), np.log10(np.max(atm.p)), len(pressure_interp)+1)

    for vol in atm.vol_list.keys():

        y = atm.p_vol[vol] 
        cs = CubicSpline(x, y) 
        atm.p_vol[vol] = np.flip(cs(np.log(pressure_interp)))

        y = atm.x_gas[vol] 
        cs = CubicSpline(x, y) 
        atm.x_gas[vol] = np.flip(cs(np.log(pressure_interp)))

        y = atm.x_cond[vol] 
        cs = CubicSpline(x, y) 
        atm.x_cond[vol] = np.flip(cs(np.log(pressure_interp)))

    atm.p = np.flip(pressure_interp)

    return atm