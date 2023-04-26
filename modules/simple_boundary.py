#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:03:32 2023

@authors: 
Mark Hammond (MH)  
Ryan Boukrouche (RB)
"""

def simple_boundary_tend(npz,pt,ts,t_dt,ts_dt,mix_coeff_atmos=1e6,mix_coeff_surf=1e6):

    t_dt[npz] += (ts + (ts - pt[npz]))/mix_coeff_atmos
    ts_dt -= (ts - pt[npz])/mix_coeff_surf
    
    return t_dt, ts_dt