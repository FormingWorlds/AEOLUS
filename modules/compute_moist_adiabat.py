#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:01:27 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
import numpy as np
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere
from modules.water_cloud import simple_cloud
from modules.relative_humidity import compute_Rh, update_for_constant_RH
from modules.interpolate_back import interpolate_back_to_10000

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates

#rh_computed = False # for the relative humidity

def compute_moist_adiabat(atm, dirs, standalone, trppD, calc_cf=False, rscatter=False, do_cloud=False):
    """Compute moist adiabat case

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py
        dirs : dict
            Named directories
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        calc_cf : bool
            Calculate contribution function?
        rscatter : bool
            Include Rayleigh scattering?
            
    """

    #global rh_computed

    atm_moist = copy.deepcopy(atm)
    
    atm_moist = ga.general_adiabat(atm_moist) # Build initial general adiabat structure

    #constant_relative_humidity = True
    #if constant_relative_humidity:

    atm_moist.rh = compute_Rh(atm_moist)

#     # Manual PT arrays for debugging:
#     atm_moist.tmp = np.array([212.812465,   214.1396898,  215.4835828,  216.84445996, 218.22264527,
#  219.618471,   221.03227801, 222.46441596, 223.91521321, 225.38502182,
#  226.87424909, 228.38329802, 229.91256698, 231.46246507, 233.03341252,
#  234.62584102, 236.24017727, 237.87684253, 239.53632097, 241.21912765,
#  242.92575792, 244.65672126, 246.41254177, 248.19375873, 250.00092709,
#  251.83456086, 253.69524662, 255.5836457,  257.50038184, 259.44609763,
#  261.42145528, 263.42713727, 265.46384725, 267.53226789, 269.63312456,
#  271.76723027, 274.06767509, 276.56162405, 279.12319219, 281.73267226,
#  284.39142105, 287.1008195,  289.86227323, 292.67732633, 295.54761194,
#  298.47477153, 301.46051231, 304.50661046, 307.61491471, 310.78735008,
#  314.02582391, 317.33242214, 320.7094174,  324.15908181, 327.68378628,
#  331.28600604, 334.96832634, 338.73344868, 342.58412062, 346.52325396,
#  350.55401223, 354.67966987, 358.90361793, 363.22941092, 367.66077685,
#  372.20162787, 376.856026,   381.62819125, 386.52270023, 391.54441511,
#  396.69835965, 401.98982567, 407.42439113, 413.00793976, 418.74668227,
#  424.64700821, 430.71584452, 436.96069113, 443.38931676, 450.00995412,
#  456.83133502, 463.8627288,  471.11398428, 478.59543404, 486.31814734,
#  494.2941347,  502.53613741, 511.05768718, 519.87324915, 528.99830375,
#  538.44943741, 548.24435278, 558.40196878, 568.94291081, 579.88948609,
#  591.26556288, 603.09691901, 615.41143665, 628.23932149, 642.39276694])

#     atm_moist.p = np.array([1.08944150e+00, 1.30273406e+00, 1.55778537e+00, 1.86277102e+00,
#  2.22746725e+00, 2.66356427e+00, 3.18504105e+00, 3.80861337e+00,
#  4.55426965e+00, 5.44591168e+00, 6.51212079e+00, 7.78707398e+00,
#  9.31163951e+00, 1.11346869e+01, 1.33146533e+01, 1.59214169e+01,
#  1.90385368e+01, 2.27659313e+01, 2.72230809e+01, 3.25528583e+01,
#  3.89261078e+01, 4.65471221e+01, 5.56601904e+01, 6.65574296e+01,
#  7.95881472e+01, 9.51700391e+01, 1.13802578e+02, 1.36083025e+02,
#  1.62725572e+02, 1.94584240e+02, 2.32680246e+02, 2.78234749e+02,
#  3.32707983e+02, 3.97846072e+02, 4.75736999e+02, 5.68877532e+02,
#  6.80253265e+02, 8.13434313e+02, 9.72689755e+02, 1.16312448e+03,
#  1.39084281e+03, 1.66314419e+03, 1.98875716e+03, 2.37811914e+03,
#  2.84371103e+03, 3.40045722e+03, 4.06620405e+03, 4.86229184e+03,
#  5.81423894e+03, 6.95255973e+03, 8.31374275e+03, 9.94142031e+03,
#  1.18877671e+04, 1.42151728e+04, 1.69982415e+04, 2.03261838e+04,
#  2.43056759e+04, 2.90642791e+04, 3.47545291e+04, 4.15588252e+04,
#  4.96952770e+04, 5.94246961e+04, 7.10589560e+04, 8.49709895e+04,
#  1.01606743e+05, 1.21499470e+05, 1.45286826e+05, 1.73731308e+05,
#  2.07744695e+05, 2.48417276e+05, 2.97052799e+05, 3.55210261e+05,
#  4.24753881e+05, 5.07912861e+05, 6.07352835e+05, 7.26261323e+05,
#  8.68449900e+05, 1.03847638e+06, 1.24179091e+06, 1.48491068e+06,
#  1.77562881e+06, 2.12326420e+06, 2.53896018e+06, 3.03604178e+06,
#  3.63044278e+06, 4.34121655e+06, 5.19114671e+06, 6.20747754e+06,
#  7.42278720e+06, 8.87603208e+06, 1.06137955e+07, 1.26917810e+07,
#  1.51765977e+07, 1.81478958e+07, 2.17009192e+07, 2.59495591e+07,
#  3.10300044e+07, 3.71051072e+07, 4.43696031e+07, 5.30563534e+07])
    

    #print("atm_moist.rh after computing = ", atm_moist.rh)
    #p_sat = [ga.p_sat('H2O', atm_moist.tmp[i]) for i in range(len(atm_moist.p))]
    #print("p_sat after computing = ", p_sat)
    #print("atm_moist.p after computing = ", atm_moist.p)
    #print("atm_moist.p_vol['H2O'] after computing = ", atm_moist.p_vol['H2O'])
    #print("atm_moist.p_vol['N2'] after computing = ", atm_moist.p_vol['N2'])
    #print("atm_moist.x_gas['H2O'] after computing = ", atm_moist.x_gas['H2O'])
    #print("atm_moist.x_gas['N2'] after computing = ", atm_moist.x_gas['N2'])
    #print("atm_moist.x_cond['H2O'] after computing = ", atm_moist.x_cond['H2O'])

        # R_h = np.ones_like(atm_moist.rh) # Evaporating water to saturation throughout the column: R_h = 1.
        # if (np.all(atm_moist.rh != R_h) and atm_moist.tmp[-1] >= ga.Tdew('H2O', atm_moist.p[-1])): 
        #     atm_moist = update_for_constant_RH(atm_moist, R_h)

        # atm_moist.rh = compute_Rh(atm_moist)

        # atm_moist = interpolate_back_to_10000(atm_moist) # interpolate the arrays back to length 10000 to build the updated general adiabat.

        # print("--------------------------------------------------------------------------")
        # print("--------------------SECOND CALL TO GENERAL ADIABAT------------------------")
        # print("--------------------------------------------------------------------------")

        # # Build general adiabat structure with updated arrays for constant relative humidity.
        # atm_moist = ga.general_adiabat(atm_moist)
        # atm_moist.rh = compute_Rh(atm_moist) 
        # print("atm_moist.rh after adjusting = ", atm_moist.rh)
        # print("p_sat after adjusting = ", p_sat)
        # print("atm_moist.p after adjusting = ", atm_moist.p)
        # print("atm_moist.p_vol['H2O'] after adjusting = ", atm_moist.p_vol['H2O'])
        # print("atm_moist.x_gas['H2O'] after adjusting = ", atm_moist.x_gas['H2O'])
        # print("atm_moist.x_cond['H2O'] after adjusting = ", atm_moist.x_cond['H2O'])

    if do_cloud:
        atm_moist = simple_cloud(atm_moist) # Before radiation, set up the cloud for Socrates using the current PT profile

    # Run SOCRATES
    atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    # Calculate tropopause
    if (trppD == True) or (atm_moist.trppT > atm_moist.minT):
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist,trppD, verbose=standalone)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        if do_cloud:
            atm_moist = simple_cloud(atm_moist) # Update cloud location after previous PT changes
        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist
