
'''
Created 01/12/19

@authors: 
Ryan Boukrouche (RB)
Tim Lichtenberg (TL)
RJ Graham (RJ)
This file builds a self-consistent temperature profile from either Li, Ingersoll & Oyafuso 2018 or 
Graham, Boukrouche, Lichtenberg, & Pierrehumbert 2020, from the ground up using the Runge-Kutta 4 scheme 
'''

import time
import numpy as np
import scipy.interpolate as spint
import math
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import utils.water_tables as wt

from utils.cp_funcs import *
from utils.ClimateUtilities import *
from utils.atmosphere_column import atmos
import utils.phys as phys

# Color definitions: 
# https://www.codecademy.com/articles/seaborn-design-ii
# https://python-graph-gallery.com/python-colors/
# https://matplotlib.org/tutorials/colors/colormaps.html
# https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
no_colors   = 7
vol_colors = {
    "H2O"            : [mpl.colormaps["PuBu"](i) for i in np.linspace(0,1.0,no_colors)],
    "CO2"            : [mpl.colormaps["Reds"](i) for i in np.linspace(0,1.0,no_colors)],
    "H2"             : [mpl.colormaps["Greens"](i) for i in np.linspace(0,1.0,no_colors)],
    "N2"             : [mpl.colormaps["Purples"](i) for i in np.linspace(0,1.0,no_colors)],
    "O2"             : [mpl.colormaps["Wistia"](i) for i in np.linspace(0,1.0,no_colors+2)],
    "CH4"            : [mpl.colormaps["RdPu"](i) for i in np.linspace(0,1.0,no_colors)],
    "CO"             : [mpl.colormaps["pink_r"](i) for i in np.linspace(0,1.0,no_colors)],
    "S"              : [mpl.colormaps["YlOrBr"](i) for i in np.linspace(0,1.0,no_colors)],
    "He"             : [mpl.colormaps["Greys"](i) for i in np.linspace(0,1.0,no_colors)],
    "NH3"            : [mpl.colormaps["summer"](i) for i in np.linspace(0,1.0,no_colors)],
    "mixtures"       : [mpl.colormaps["Set3"](i) for i in np.linspace(0,1.0,no_colors)],
    "H2O-CO2"        : mpl.colormaps["Set3"](1.0/no_colors),
    "CO2-H2O"        : mpl.colormaps["Set3"](1.0/no_colors),
    "H2O-H2"         : mpl.colormaps["Set3"](2.0/no_colors),
    "H2-H2O"         : mpl.colormaps["Set3"](2.0/no_colors),
    "CO-H2"          : mpl.colormaps["Set3"](3.0/no_colors),
    "H2-CO"          : mpl.colormaps["Set3"](3.0/no_colors),
    "H2-CO2"         : mpl.colormaps["Set3"](4.0/no_colors),
    "CO2-H2"         : mpl.colormaps["Set3"](4.0/no_colors),
    "CH4-H2"         : mpl.colormaps["Set3"](5.0/no_colors),
    "H2-CH4"         : mpl.colormaps["Set3"](5.0/no_colors),
    "H2-CH4"         : mpl.colormaps["Set3"](5.0/no_colors),
    "H2-N2"          : mpl.colormaps["Set2"](0.0/no_colors),
    "N2-H2"          : mpl.colormaps["Set2"](0.0/no_colors),
    "CO2-N2"         : mpl.colormaps["Set2"](1.0/no_colors),
    "N2-CO2"         : mpl.colormaps["Set2"](1.0/no_colors),
    "black_1"        : "#000000",
    "black_2"        : "#323232",
    "black_3"        : "#7f7f7f",
    "qgray"          : "#768E95",
    "qgray2"         : "#888888",
    "qblue"          : "#4283A9", # http://www.color-hex.com/color/4283a9
    "qgreen"         : "#62B4A9", # http://www.color-hex.com/color/62b4a9
    "qred"           : "#E6767A",
    "qturq"          : "#2EC0D1",
    "qorange"        : "#ff7f0e",
    "qmagenta"       : "#9A607F",
    "qyellow"        : "#EBB434",
    "qgray_dark"     : "#465559",
    "qblue_dark"     : "#274e65",
    "qgreen_dark"    : "#3a6c65",
    "qred_dark"      : "#b85e61",
    "qturq_dark"     : "#2499a7",
    "qmagenta_dark"  : "#4d303f",
    "qyellow_dark"   : "#a47d24",
    "qgray_light"    : "#acbbbf",
    "qblue_light"    : "#8db4cb",
    "qgreen_light"   : "#a0d2cb",
    "qred_light"     : "#eb9194",
    "qturq_light"    : "#57ccda",
    "qmagenta_light" : "#c29fb2",
    "qyellow_light"  : "#f1ca70",
}

# Volatile Latex names
vol_latex = {
    "H2O"     : r"H$_2$O",
    "CO2"     : r"CO$_2$",
    "H2"      : r"H$_2$" ,
    "CH4"     : r"CH$_4$",
    "CO"      : r"CO",
    "N2"      : r"N$_2$",
    "S"       : r"S",
    "O2"      : r"O$_2$",
    "He"      : r"He",
    "NH3"     : r"NH$_3$",
    "H2O-CO2" : r"H$_2$O–CO$_2$",
    "H2O-H2"  : r"H$_2$O–H$_2$",
    "H2O-CO"  : r"H$_2$O–CO",
    "H2O-CH4" : r"H$_2$O–CH$_4$",
    "H2O-N2"  : r"H$_2$O–N$_2$",
    "H2O-O2"  : r"H$_2$O–O$_2$",
    "H2-H2O"  : r"H$_2$–H$_2$O",
    "H2-CO"   : r"H$_2$–CO",
    "H2-CH4"  : r"H$_2$–CH$_4$",
    "H2-CO2"  : r"H$_2$–CO$_2$",
    "H2-N2"   : r"H$_2$–N$_2$",
    "H2-O2"   : r"H$_2$-O$_2$",
    "CO2-N2"  : r"CO$_2$–N$_2$",
    "CO2-H2O" : r"CO$_2$–H$_2$O",
    "CO2-CO"  : r"CO$_2$–CO",
    "CO2-CH4"  : r"CO$_2$–CH$_4$",
    "CO2-O2"  : r"CO$_2$–O$_2$",
    "CO2-H2"  : r"CO$_2$–H$_2$",
    "CO-H2O" : r"CO–H$_2$O",
    "CO-CO2" : r"CO–CO$_2$",
    "CO-H2"  : r"CO–H$_2$",
    "CO-CH4" : r"CO–CH$_4$",
    "CO-N2"  : r"CO–N$_2$",
    "CO-O2"  : r"CO–O$_2$",
    "CH4-H2O" : r"CH$_4$–H$_2$O",
    "CH4-CO2" : r"CH$_4$–CO$_2$",
    "CH4-H2"  : r"CH$_4$–H$_2$",
    "CH4-CO"  : r"CH$_4$–CO",
    "CH4-CH4" : r"CH$_4$–CH$_4$",
    "CH4-N2"  : r"CH$_4$–N$_2$",
    "CH4-O2"  : r"CH$_4$–O$_2$",
    "N2-H2O" : r"N$_2$–H$_2$O",
    "N2-CO2" : r"N$_2$–CO$_2$",
    "N2-H2"  : r"N$_2$–H$_2$",
    "N2-CO"  : r"N$_2$–CO",
    "N2-CH4" : r"N$_2$–CH$_4$",
    "N2-N2"  : r"N$_2$–N$_2$",
    "N2-O2"  : r"N$_2$–O$_2$",
    "O2-H2O" : r"O$_2$–H$_2$O",
    "O2-CO2" : r"O$_2$–CO$_2$",
    "O2-H2"  : r"O$_2$–H$_2$",
    "O2-CO"  : r"O$_2$–CO",
    "O2-CH4" : r"O$_2$–CH$_4$",
    "O2-N2"  : r"O$_2$–N$_2$",
    "O2-O2"  : r"O$_2$–O$_2$",
}

def atm_z(atm, idx):

    # Mean temperature below
    #T_mean_down = np.mean( atm.tmp[:idx+1] )

    # # Average temperature weighted by amount of substance?
    # T_mean_down = np.average( atm.tmp[:idx+1], weights=)

    # # Use gas phase mean molar mass
    # atm.mu[idx] = atm.mu[idx]

    # print(atm.grav_z[idx], atm.mu[idx], atm.p[idx], atm.p[idx+1])
    

    mu_c = 0.
    if atm.xc[idx] > 0:
        for vol in atm.vol_list.keys():
            mu_c += atm.x_cond[vol][idx]*phys.molar_mass[vol] / atm.xc[idx]
    
    # Integration
    atm.rho[idx] =  atm.p[idx]/phys.R_gas/atm.tmp[idx] * ( atm.mu[idx] + mu_c*atm.alpha_cloud*atm.xc[idx]/(atm.xv[idx]+atm.xd[idx]) )
    #     cell centres
    dp = atm.p[idx+1]-atm.p[idx]
    dz =  -dp / atm.rho[idx]
    atm.z[idx+1] = atm.z[idx] + dz
    #     cell edges 
    dp = atm.pl[idx+1]-atm.pl[idx]
    dz =  -dp / atm.rho[idx]
    atm.zl[idx+1] = atm.zl[idx] + dz
    if dz<0:
        print("WARNING: dz  = %g < 0 " % dz)
        print("         (dp  = %g)" % dp)
        print("         (rho = %g)" % atm.rho[idx])
    
    # Next gravity
    atm.grav_z[idx+1] = atm.grav_s * ((atm.planet_radius)**2) / ((atm.planet_radius+atm.z[idx+1])**2)

    # Calculate scale height in km
    # H = 1e-3 * phys.R_gas * atm.tmp[idx] / (atm.mu[idx] * atm.grav_z[idx])
    # print(H)

    # print(idx, T_mean_down, dz, atm.z[idx], atm.z[idx+1], atm.grav_z[idx], atm.grav_z[idx+1])

    return atm

## Saturation vapor pressure [Pa] for given temperature T [K]. 
## Assuming the ideal gas law and a constant latent heat of vaporization. 
## Select the molecule of interest with the switch argument (a string).
def p_sat(switch,T, water_lookup=False): 

    # Define volatile
    match switch:
        case 'H2O':
            e = phys.satvps_function(phys.water)
        case 'CH4':
            e = phys.satvps_function(phys.methane)
        case 'CO2':
            e = phys.satvps_function(phys.co2)
        case 'CO':
            e = phys.satvps_function(phys.co)
        case 'N2':
            e = phys.satvps_function(phys.n2)
        case 'O2':
            e = phys.satvps_function(phys.o2)
        case 'H2':
            e = phys.satvps_function(phys.h2)
        case 'He':
            e = phys.satvps_function(phys.he)
        case 'NH3':
            e = phys.satvps_function(phys.nh3)   
        case _:
            e = lambda T,water_lookup: 0.0
        
    # Return saturation vapor pressure
    svp = e(T,water_lookup=water_lookup)
    return float(svp)


## Dew point temperature [K] given a pressure p [Pa]. Select the molecule of interest with the switch argument (a string).
def Tdew(switch, p): 
    
    # Avoid math error for p = 0
    p = np.max([p, 1e-100])

    if switch == 'H2O':
        p_triple = 611.657 # Pa
        if p > p_triple:
            Tref = 373.15 # K, boiling point of H2O at 1 atm 
            pref = 1e5 # esat('H2O',Tref) returns 121806.3 Pa, should return 1e5    
            L_H2O=phys.water.L_vaporization*phys.water.MolecularWeight*1e-3
        else:
            Tref = 273.15 #K triple point temperature of H2O
            pref = p_triple #triple point pressure of H2O
            L_H2O=phys.water.L_sublimation*phys.water.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_H2O)*math.log(p/pref))
        #return (B(p)/(A(p)-numpy.log10(p/pref)))-C(p)
    if switch == 'CH4':
        Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
        pref = p_sat('CH4',Tref)
        L_CH4=phys.CH4.L_vaporization*phys.methane.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CH4)*math.log(p/pref))
    if switch == 'CO2':
        p_triple = 5.11e5
        if p > p_triple:
            Tref = 253. # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
            pref = p_sat('CO2',Tref)
            L_CO2=phys.CO2.L_vaporization*phys.co2.MolecularWeight*1e-3
        else:
            Tref = 216.58 # K; triple point
            pref = p_triple
            L_CO2 = phys.CO2.L_sublimation*phys.co2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CO2)*math.log(p/pref))
    if switch == 'CO':
        Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
        pref = p_sat('CO',Tref)
        L_CO=phys.CO.L_vaporization*phys.co.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_CO)*math.log(p/pref))
    if switch == 'N2':
        Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
        pref = p_sat('N2',Tref)
        L_N2=phys.N2.L_vaporization*phys.n2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_N2)*math.log(p/pref))
    if switch == 'O2':
        Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
        pref = p_sat('O2',Tref)
        L_O2=phys.O2.L_vaporization*phys.o2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_O2)*math.log(p/pref))
    if switch == 'H2':
        Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
        pref = p_sat('H2',Tref)
        L_H2=phys.H2.L_vaporization*phys.h2.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_H2)*math.log(p/pref))
    if switch == 'He':
        Tref = 4.22 # K, boiling point of He at 1 atm 
        pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5
        L_He=phys.He.L_vaporization*phys.he.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_He)*math.log(p/pref))
    if switch == 'NH3':
        Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
        pref = p_sat('NH3',Tref)
        L_NH3=phys.NH3.L_vaporization*phys.nh3.MolecularWeight*1e-3
        return Tref/(1.-(Tref*phys.R_gas/L_NH3)*math.log(p/pref))
    
    return 0.0


def get_beta(switch, T, water_lookup=False):
    if (switch == 'H2O') and water_lookup \
        and (T > phys.H2O.TriplePointT) \
        and (T< phys.H2O.CriticalPointT):
        beta = wt.lookup('phase_grad', T)
    else:
        L = L_heat(switch, T)
        beta = L/phys.R_gas/T

    return beta
## Molar latent heat [J mol-1] for gas phase considered given a temperature T [K]. 
## Select the molecule of interest with the switch argument (a string).

def get_T_crit(switch):
    match switch:
        case 'H2O':
            return phys.H2O.CriticalPointT
        case 'CH4':
            return phys.CH4.CriticalPointT
        case 'CO2':
            return phys.CO2.CriticalPointT
        case 'CO':
            return phys.CO.CriticalPointT
        case 'N2':
            return phys.N2.CriticalPointT
        case 'O2':
            return phys.O2.CriticalPointT
        case 'H2':
            return phys.H2.CriticalPointT
        case 'He':
            return phys.He.CriticalPointT
        case 'NH3':
            return phys.NH3.CriticalPointT
        case _:
            return 0.0
        
def L_heat(switch, T, water_lookup=False):

    match switch:
        case 'H2O':
            L_sublimation   = phys.H2O.L_sublimation
            
            if water_lookup:
                L_vaporization  = wt.lookup('L_vap', T)
            else:
                L_vaporization  = phys.H2O.L_vaporization

            MolecularWeight = phys.H2O.MolecularWeight
            T_triple        = phys.H2O.TriplePointT
            T_crit          = phys.H2O.CriticalPointT
        
        case 'CH4':
            L_sublimation   = phys.CH4.L_sublimation
            L_vaporization  = phys.CH4.L_vaporization
            MolecularWeight = phys.CH4.MolecularWeight
            T_triple        = phys.CH4.TriplePointT
            T_crit          = phys.CH4.CriticalPointT
            
        case 'CO2':
            L_sublimation   = phys.CO2.L_sublimation
            L_vaporization  = phys.CO2.L_vaporization
            MolecularWeight = phys.CO2.MolecularWeight
            T_triple        = phys.CO2.TriplePointT
            T_crit          = phys.CO2.CriticalPointT
            
        case 'CO':
            L_sublimation   = phys.CO.L_sublimation
            L_vaporization  = phys.CO.L_vaporization
            MolecularWeight = phys.CO.MolecularWeight
            T_triple        = phys.CO.TriplePointT
            T_crit          = phys.CO.CriticalPointT
                
        case 'N2':
            L_sublimation   = phys.N2.L_sublimation
            L_vaporization  = phys.N2.L_vaporization
            MolecularWeight = phys.N2.MolecularWeight
            T_triple        = phys.N2.TriplePointT
            T_crit          = phys.N2.CriticalPointT
            
        case 'O2':
            L_sublimation   = phys.O2.L_sublimation
            L_vaporization  = phys.O2.L_vaporization
            MolecularWeight = phys.O2.MolecularWeight
            T_triple        = phys.O2.TriplePointT
            T_crit          = phys.O2.CriticalPointT
        
        case 'H2':
            L_sublimation   = phys.H2.L_vaporization # No H2 sublimation
            L_vaporization  = phys.H2.L_vaporization
            MolecularWeight = phys.H2.MolecularWeight
            T_triple        = phys.H2.TriplePointT
            T_crit          = phys.H2.CriticalPointT
            
        case 'He':
            L_sublimation   = phys.He.L_vaporization  # No He sublimation
            L_vaporization  = phys.He.L_vaporization
            MolecularWeight = phys.He.MolecularWeight
            T_triple        = phys.He.TriplePointT
            T_crit          = phys.He.CriticalPointT
            
        case 'NH3':
            L_sublimation   = phys.NH3.L_sublimation
            L_vaporization  = phys.NH3.L_vaporization
            MolecularWeight = phys.NH3.MolecularWeight
            T_triple        = phys.NH3.TriplePointT
            T_crit          = phys.NH3.CriticalPointT
        case _:
            L_sublimation   = 0.0
            L_vaporization  = 0.0
            MolecularWeight = 0.0
            T_triple        = 0.0
            T_crit          = 0.0

    # Gas-solid transition
    if T <= T_triple:
        # Conversion from J.kg-1 to J.mol-1 (molecular weight is in g/mol in phys.f90)
        L_heat = L_sublimation*MolecularWeight*1e-3 
    # Gas-liquid transition
    elif T <= T_crit:
        L_heat = L_vaporization*MolecularWeight*1e-3 
    # Super-critical state
    else:
        # Numerical issus with 0.?
        L_heat = 0.         #1e-50
        # L_heat = L_vaporization*MolecularWeight*1e-3 

    # # No latent heat contribution in moist adiabat if below p_sat
    # if P < p_sat(switch, T):
    # # if T > Tdew(switch, psat):
    #     L_heat = 0.

    return L_heat  


satvph2o = phys.satvps_function(phys.H2O)
def slopeRay( logpa, logT ):
    eps     = phys.H2O.MolecularWeight/phys.air.MolecularWeight
    L       = phys.H2O.L_vaporization
    Ra      = phys.air.R
    Rc      = phys.H2O.R
    cpa     = phys.air.cp
    cpc     = phys.H2O.cp
    pa      = math.exp(logpa)
    T       = math.exp(logT)
    qsat    = eps*(satvph2o(T)/pa)
    print('qsat=%.3f'%qsat)
    num     = (1. + (L/(Ra*T))*qsat)*Ra
    print('numerator=%.3f'%num)
    den     = cpa + (cpc + (L/(Rc*T) - 1.)*(L/T))*qsat
    print('denominator=%.3f'%den)
    print('dlnT/dlnPa=%.3f'%(num/den))
    return num/den


# Dry adiabat function from p. 91 of RTP book as a comparison
def dry_adiabat( T_surf, P_array, cp_array ):

    T_dry   = np.zeros(len(P_array))
    P_surf  = np.amax(P_array)

    # Check if scalar or cp array
    if not isinstance(cp_array, list):
        cp_array = np.ones(len(P_array))*cp_array

    for idx, prs in enumerate(P_array):

        if cp_array[idx] > 0.:

            T_dry[idx] = T_surf * ( prs / P_surf ) ** ( phys.R_gas / cp_array[idx] )

    return T_dry

# Function for determining the dry adiabat pressure at a given temperature 
def dry_adiabat_pressure( P_surf, T_array, cp_array ):

    P_dry   = np.zeros(len(T_array))
    T_surf  = np.amax(T_array)

    # Check if scalar or cp array
    if not isinstance(cp_array, list):
        cp_array = np.ones(len(T_array))*cp_array

    for idx, tmp in enumerate(T_array):

        if cp_array[idx] > 0.:

            P_dry[idx] = P_surf * ( tmp / T_surf ) ** ( cp_array[idx] / phys.R_gas  )

    return P_dry

def dlnT_dlnP_d(lnP, lnT, atm):
    # T instead lnT
    tmp = math.exp(lnT)


    idx = int(np.amax(atm.ifatm))
    
    
    # Sum terms in equation
    num_sum     = 0.
    denom_sum  = 0. 

    
    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 
        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        eta_cond    = atm.x_cond[vol][idx] / atm.xd[idx]
        #print(eta_vol)
        # sums for saturated comps
        # HII added, ensure below the critical point !
        if ((np.isclose(atm.p_vol[vol][idx] ,p_sat(vol,tmp, water_lookup=atm.water_lookup)) or
             atm.p_vol[vol][idx]  > p_sat(vol,tmp, water_lookup=atm.water_lookup)) and
            tmp < get_T_crit(vol)):

            beta = get_beta(vol, tmp, water_lookup=atm.water_lookup)
            L = L_heat(vol,tmp, water_lookup=atm.water_lookup)
            L_RT = L/phys.R_gas/tmp
            #beta_vol    = L_heat(vol, tmp) / (R_gas * tmp) #RTP
            # Sum in numerator
            num_sum     += eta_vol * L_RT*phys.R_gas
            # Sum in denominator
            denom_sum += eta_vol*(cpv(vol,tmp)  + phys.R_gas*beta*(L_RT - 1) +
                                  atm.alpha_cloud * eta_cond * cp_cond(vol,tmp))
            
        # sums for subsaturated comps
        else: 
            #print(eta_vol)
            denom_sum += cpv(vol,tmp) * eta_vol #The eta_vol (x_vol/x_d) is there so that if there are multiple dry components the average dry specific heat is used
        
    
        
    
    # Collect terms
    numerator   = phys.R_gas + num_sum
    denominator = denom_sum

    # dlnT/dlnPd
    dlnTdlnPd = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnPd


# In[70]:


def moist_slope(lnP, lnT, atm):
    # T instead lnT
    tmp = math.exp(lnT)
    
    idx = int(np.amax(atm.ifatm))
    # Sum terms in equation
    num_sum     = 0.
    denom_sum  = 0. 

    
    
    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 
        
        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        
        #print(eta_vol)
        # sums for saturated comps
        # Ensure below the critical point for volatile before we
        # assume condensation occurs
        if ((np.isclose(atm.p_vol[vol][idx] ,p_sat(vol,tmp,water_lookup=atm.water_lookup)) or
             atm.p_vol[vol][idx]  > p_sat(vol,tmp,water_lookup=atm.water_lookup)) and
            tmp < get_T_crit(vol)):
            
            # Sum in numerator
            num_sum     += eta_vol
            
            # Sum in denominator
            beta = get_beta(vol, tmp, water_lookup=atm.water_lookup)
            denom_sum  += beta * eta_vol * dlnT_dlnP_d(lnP,lnT,atm)
             
    
    # Collect terms
    numerator   = 1 + num_sum
    denominator = 1 + denom_sum
    dlnPd_dlnP = numerator / denominator
    # dlnT/dlnPd
    dlnTdlnP = dlnPd_dlnP * dlnT_dlnP_d(lnP,lnT,atm)

    # Moist adiabat slope
    return dlnTdlnP




def condensation( atm, idx, wet_list, dry_list, prs_reset):

    # Temperature floor
    tmp = np.amax([atm.tmp[idx], atm.minT])
    
    
    if idx==0:
        
        # Renormalize mixing ratios to ensure sum == 1
        if sum(atm.vol_list.values()) != 1.:
            vol_list_new = {}
            for vol in atm.vol_list.keys():
                vol_list_new[vol]   = atm.vol_list[vol] / sum(atm.vol_list.values())
            for vol in atm.vol_list.keys():
                atm.vol_list[vol] = vol_list_new[vol]
        # Partial pressures scaled 
        for vol in atm.vol_list:
            atm.p_vol[vol][idx] = atm.vol_list[vol] * atm.p[idx]
            #atm.p_vol[vol][idx] = atm.vol_list[vol] * p_tot_pre_condensation
            # Mean gas phase molar mass
            atm.mu[idx]       += atm.vol_list[vol] * phys.molar_mass[vol]
    
    else:  
        
        # Total pressure & partial pressures of condensing phases
        p_cond_sum = 0
        for vol in atm.vol_list:
            if vol in wet_list:
                atm.p_vol[vol][idx] = p_sat(vol, atm.tmp[idx],water_lookup=atm.water_lookup)
                p_cond_sum += atm.p_vol[vol][idx]
                
        # Calculate the total partial pressure of dry species
        p_dry_tot = atm.p[idx] - p_cond_sum
        dry_frac_sum = 0
        # calculate sum of fractions of dry species for partial pressure calculation
        for vol in atm.vol_list:
            if vol in dry_list:
                dry_frac_sum += atm.x_gas[vol][idx-1]
        # Calculate the individual partial pressures of dry species
        for vol in atm.vol_list:
            if vol in dry_list:
                #print(atm.x_gas[vol][idx-1])
                atm.p_vol[vol][idx] = p_dry_tot * (atm.x_gas[vol][idx-1]/dry_frac_sum)
                
       
    # Calculate mean molar mass
    atm.mu[idx]   = 0.
    
    for vol in atm.vol_list.keys():
        atm.mu[idx]   += phys.molar_mass[vol] * atm.p_vol[vol][idx]
    atm.mu[idx] /= atm.p[idx]
        
    

    
    # Update condensate mixing ratios
    for vol in atm.vol_list.keys():
       
        # Condensate phase
        if atm.p_vol[vol][idx] < p_sat(vol, tmp,water_lookup=atm.water_lookup):
            atm.x_cond[vol][idx] = 0.
        else:
            #atm.x_cond[vol][idx] = atm.vol_list[vol] - ( atm.p_vol[vol][idx] / atm.p[idx] )
            if idx == 0:
                #atm.x_cond[vol][idx] = atm.vol_list[vol] - ( mu_old / atm.mu[idx] * atm.p_vol[vol][idx] / p_tot_pre_condensation )
                atm.x_cond[vol][idx] = 0.
            else:
                #mu_old = atm.mu[idx-1]
                atm.x_cond[vol][idx] = atm.x_cond[vol][idx-1] + atm.x_gas[vol][idx-1] - (1-atm.xc[idx-1])*(atm.p_vol[vol][idx]/atm.p[idx])

        #atm.x_cond[vol][idx] *= atm.alpha_cloud

        # Add to molar concentration of total condensed phase
        atm.xc[idx]     += atm.x_cond[vol][idx]
    for vol in atm.vol_list.keys():    
        # Gas phase molar concentration
        #atm.x_gas[vol][idx] = atm.p_vol[vol][idx] / atm.p[idx]
        if idx == 0:
            #atm.x_gas[vol][idx] =  mu_old / atm.mu[idx] * atm.p_vol[vol][idx] /  p_tot_pre_condensation
            atm.x_gas[vol][idx] = atm.vol_list[vol]
            
        else:
            atm.x_gas[vol][idx] =  ( 1-atm.xc[idx] ) * atm.p_vol[vol][idx] /  atm.p[idx]
            #atm.x_gas[vol][idx] =  mu_old / atm.mu[idx] * atm.p_vol[vol][idx] /  p_tot_pre_condensation
        
        # Add to molar concentration of total gas (dry or moist) phase
        # ! REVISIT ! keeping xd == 0 leads to a bug, why?
        if vol in dry_list:
            atm.xd[idx]          += atm.x_gas[vol][idx]
        if vol in wet_list:
            atm.xv[idx]          += atm.x_gas[vol][idx]
        
        # Mean cp of both gas phase and retained condensates
        atm.cp[idx]    += atm.x_gas[vol][idx] * cpv(vol, atm.tmp[idx]) + atm.x_cond[vol][idx] * cp_cond(vol, atm.tmp[idx]) * atm.alpha_cloud

    # Renormalize cp w/ molar concentration (= cp_hat, Eq. 1, Li+2018)
    atm.cp[idx]  = atm.cp[idx] / ( atm.xd[idx] + atm.xv[idx] )

    # Dry concentration floor
    atm.xd[idx]  = np.amax([atm.xd[idx], 1e-10])
    
    # Loop through species to determine wet_list and dry_list for next level
    
    for vol in atm.vol_list.keys():
        
        
        # Condensation if p_i > p_sat
        if atm.p_vol[vol][idx] >= p_sat(vol, tmp,water_lookup=atm.water_lookup):
            
            # Add condensing species to wet_list
            if vol not in wet_list:
                wet_list.append(vol)
                
            
            # Set species partial pressure to p_sat
            #atm.p_vol[vol][idx]  = p_sat(vol, tmp)
            
            # Add the species partial pressure to the condensing species sum
            #p_cond_sum += atm.p_vol[vol][idx]
            
            # If the species is in dry_list from previous iteration of while loop, remove it
            if vol in dry_list:
                
                dry_list.remove(vol)
            
        else:
            
            # Add a non-condensing species to the list of dry species
            # if it's present (>0) and if it's not in a list already
            if atm.vol_list[vol] > 0. and vol not in dry_list and vol not in wet_list:
                dry_list.append(vol)
    
    return atm, wet_list, dry_list

# Builds the generalized moist adiabat from slope(lnP, lnT, atm object)
def general_adiabat( atm ):
    # First, we check to see if the initial vol_list/Tsurf/Psurf  
    # combo gives a super-saturated surface state. If so, we remove material
    # until the super-saturated phase is just saturated; none of this condensate
    # is retained because the initial supersaturated state is unrealistic.
    # Adjust the partial pressures of the supersaturated species to their
    # saturation vapor pressures, then adjust Vol_list and Psurf.
    new_psurf = 0
    new_p_vol = {}
    wet_list = []
    dry_list = []
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] * atm.ps > p_sat(vol, atm.ts,water_lookup=atm.water_lookup):
            new_psurf += p_sat(vol,atm.ts,water_lookup=atm.water_lookup)
            new_p_vol[vol] = p_sat(vol,atm.ts,water_lookup=atm.water_lookup)
        else:
            new_psurf += atm.vol_list[vol] * atm.ps
            new_p_vol[vol] = atm.vol_list[vol] * atm.ps
            
    if new_psurf != atm.ps:
        Tsurf = atm.ts
        alpha = atm.alpha_cloud
        toa_heating = atm.toa_heating
        tmp_magma = atm.tmp_magma
        minT = atm.minT
        nlev_save = atm.nlev_save
        skin_d = atm.skin_d
        skin_k = atm.skin_k
        effective_radius = atm.effective_radius
        liquid_water_fraction = atm.liquid_water_fraction
        cloud_fraction = atm.cloud_fraction

        for vol in atm.vol_list.keys():
            atm.vol_list[vol] = new_p_vol[vol] / new_psurf

        atm = atmos(Tsurf, new_psurf, atm.ptop, atm.planet_radius, atm.planet_mass, effective_radius, liquid_water_fraction, cloud_fraction, vol_mixing=atm.vol_list, trppT=atm.trppT, minT=minT, req_levels=nlev_save)
        
        atm.alpha_cloud = alpha
        atm.toa_heating = toa_heating
        atm.tmp_magma = tmp_magma
        atm.skin_d = skin_d
        atm.skin_k = skin_k
        atm.effective_radius = effective_radius
        atm.liquid_water_fraction = liquid_water_fraction
        atm.cloud_fraction = cloud_fraction
        
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] * atm.ps == p_sat(vol,atm.ts,water_lookup=atm.water_lookup):
            wet_list.append(vol)
        elif atm.vol_list[vol] * atm.ps != p_sat(vol,atm.ts,water_lookup=atm.water_lookup) and atm.vol_list[vol] > 0:
            dry_list.append(vol)
    
    ### Initialization
    
    # Initialize the tuple solution
    moist_tuple     = [] #[tuple([np.log(atm.ps), atm.ts])] 
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.01

    # Integration counter
    idx             = 0  
    
    # Surface condensation
    atm, wet_list, dry_list  = condensation(atm, idx, wet_list, dry_list, prs_reset=False)

    # Create the integrator instance                                              
    int_slope       = integrator(moist_slope, np.log(atm.ps), np.log(atm.ts), step)

    # Update parameters used in the slope function dlntT/dlnP
    int_slope.setParams(atm)
    
    ### Integration of full general adiabat
    while atm.p[idx] > atm.ptop:
        
        # Execute the Runge-Kutta integrator, fill array of tuples
        moist_tuple.append(int_slope.next())

        # Fill next T,P values
        atm.p[idx+1]    = np.exp(int_slope.x)
        atm.tmp[idx+1]  = np.exp(int_slope.y)

        # Calculate next local gravity and height
        atm = atm_z(atm, idx)

        # print("RK4 step, idx:", idx, round(atm.p[idx+1],5), round(atm.tmp[idx+1],5))

        # Set next level to calculate
        idx             += 1
        atm.ifatm[idx]  = idx

        # Calculate condensation at next level
        atm, wet_list, dry_list    = condensation(atm, idx, wet_list, dry_list, prs_reset=False)

    # Interpolate
    atm = interpolate_atm(atm)
    
    # Rescale mixing ratios for plotting purposes
    atm.xd[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    atm.xv[:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] > 1e-10:
            atm.x_cond[vol][:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
            atm.x_gas[vol][:] *= 1 / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    atm.xc[:] *= atm.alpha_cloud / ( 1 - (1-atm.alpha_cloud) * atm.xc[:])
    
    
    return atm

# Interpolate and flip pressure, temperature and volatile grids to fixed size
def interpolate_atm(atm):

    # Trim array zeros
    atm_len     = int(np.max(atm.ifatm)+1)
    rest_len    = int(len(atm.p)-atm_len)
    atm.p       = np.flip(np.split(atm.p, [atm_len, rest_len])[0], axis=0)
    atm.tmp     = np.flip(np.split(atm.tmp, [atm_len, rest_len])[0], axis=0)

    # Interpolate staggered nodes
    atm.pl      = np.logspace(np.log10(np.min(atm.p)), np.log10(np.max(atm.p)), atm.nlev_save+1)
    atm.tmpl    = np.interp(atm.pl, atm.p, atm.tmp)

    # Interpolate atmosphere nodes
    prs_itp     = (atm.pl[1:] + atm.pl[:-1]) / 2.
    tmp_itp     = (atm.tmpl[1:] + atm.tmpl[:-1]) / 2.

    # Trim level-dependent quantities
    atm.grav_z  = np.flip(np.split(atm.grav_z, [atm_len,rest_len])[0], axis=0)
    atm.z       = np.flip(np.split(atm.z, [atm_len,rest_len])[0], axis=0)
    atm.mu      = np.flip(np.split(atm.mu, [atm_len,rest_len])[0], axis=0)
    atm.xd      = np.flip(np.split(atm.xd, [atm_len, rest_len])[0], axis=0)
    atm.xv      = np.flip(np.split(atm.xv, [atm_len, rest_len])[0], axis=0)
    atm.xc      = np.flip(np.split(atm.xc, [atm_len, rest_len])[0], axis=0)
    atm.cp      = np.flip(np.split(atm.cp, [atm_len, rest_len])[0], axis=0)
    atm.rho     = np.flip(np.split(atm.rho, [atm_len, rest_len])[0], axis=0)
    # Interpolate level-dependent quantities
    atm.grav_z  = np.interp(prs_itp, atm.p, atm.grav_z)
    atm.z       = np.interp(prs_itp, atm.p, atm.z)
    atm.mu      = np.interp(prs_itp, atm.p, atm.mu)
    atm.xd      = np.interp(prs_itp, atm.p, atm.xd)
    atm.xv      = np.interp(prs_itp, atm.p, atm.xv)
    atm.xc      = np.interp(prs_itp, atm.p, atm.xc)
    atm.cp      = np.interp(prs_itp, atm.p, atm.cp)
    atm.rho     = np.interp(prs_itp, atm.p, atm.rho)
    # Trim & interpolate species-dependent quantities
    for vol in atm.vol_list.keys():
        
        # Trim
        atm.p_vol[vol]   = np.flip(np.split(atm.p_vol[vol], [atm_len, rest_len])[0], axis=0)
        atm.x_gas[vol]   = np.flip(np.split(atm.x_gas[vol], [atm_len, rest_len])[0], axis=0)
        atm.x_cond[vol]  = np.flip(np.split(atm.x_cond[vol], [atm_len, rest_len])[0], axis=0)
        
        # Interpolate
        atm.pl_vol[vol]  = np.interp(atm.pl, atm.p, atm.p_vol[vol]) # staggered!
        atm.p_vol[vol]   = np.interp(prs_itp, atm.p, atm.p_vol[vol])
        atm.x_gas[vol]   = np.interp(prs_itp, atm.p, atm.x_gas[vol]) 
        atm.x_cond[vol]  = np.interp(prs_itp, atm.p, atm.x_cond[vol])

    # Rewrite atmosphere nodes
    atm.p       = prs_itp
    atm.tmp     = tmp_itp

    return atm


def plot_adiabats(atm,filename='output/general_adiabat.pdf'):

    # sns.set_style("ticks")
    # sns.despine()

    ls_moist    = 2.5
    ls_dry      = 2.0
    ls_ind      = 1.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))#, sharey=True)
    # sns.set_style("ticks")
    # sns.despine()
    
    
    # For reference p_sat lines
    T_sat_array    = np.linspace(20,3000,1000) 
    p_partial_sum  = np.zeros(len(atm.tmp))
    
    vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}

    # Individual species
    for vol in vol_list_sorted.keys():

        # Only if volatile is present
        if atm.vol_list[vol] > 1e-10:
    
            # Saturation vapor pressure for given temperature
            Psat_array = [ p_sat(vol, T,water_lookup=atm.water_lookup) for T in T_sat_array ]
            ax1.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+vol_latex[vol], lw=ls_ind, ls=":", color=vol_colors[vol][4])

            # Plot partial pressures
            ax1.semilogy(atm.tmp, atm.p_vol[vol], color=vol_colors[vol][4], lw=ls_ind, ls="-", label=r'$p$'+vol_latex[vol],alpha=0.99)

            # Sum up partial pressures
            p_partial_sum += atm.p_vol[vol]

            # Plot individual molar concentrations
            ax2.semilogy(atm.x_cond[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="--", label=vol_latex[vol]+" cond.")
            ax2.semilogy(atm.x_gas[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="-", label=vol_latex[vol]+" gas")
            
    # # Plot sum of partial pressures as check
    ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="--", label=r'$\sum p^\mathrm{i}$',alpha=0.99)
    
    # # # Dry adiabat function from RTB book
    # # ! careful: non-T dependent CPs used, can lead to differences with general adiabat !
    ax1.semilogy( dry_adiabat( atm.ts, atm.pl, atm.cp[-1]), atm.pl , color=vol_colors["black_3"], ls="-.", lw=ls_dry, label=r'Dry adiabat') # Functional form

    # General moist adiabat
    ax1.semilogy(atm.tmpl, atm.pl, color=vol_colors["black_1"], lw=ls_moist,label="General adiabat",alpha=0.99)

    # Phase molar concentrations
    ax2.semilogy(atm.xd+atm.xv,atm.p, color=vol_colors["black_2"], lw=ls_ind, ls=":", label=r"Gas phase")

    fs_l = 16
    fs_m = 14
    fs_s = 12

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature, $T$ (K)', fontsize=fs_l)
    ax1.set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    # ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes', fontsize=fs_l)
    ax1.legend(loc=1, ncol=np.min([len(atm.vol_list)+1,2]), fontsize=fs_s)
    ax1.set_xlim([0,np.max(atm.ts)])

    ax2.invert_yaxis()
    # ax2.set_title('Phase & species abundances', fontsize=fs_l)
    ax2.set_xlabel(r'Molar concentration, $X^{\mathrm{i}}_{\mathrm{phase}}$', fontsize=fs_l)
    ax2.set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    ax2.legend(loc=2, ncol=2, fontsize=fs_s)

    ax1.set_ylim(top=atm.ptop)
    ax1.set_ylim(bottom=atm.ps)
    ax2.set_ylim(top=atm.ptop)
    ax2.set_ylim(bottom=atm.ps)

    ax2.set_xscale("log")
    ax2.set_xlim([1e-4, 1.05])
    ax2.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax2.set_xticklabels(["$10^{-4}$", "0.001", "0.01", "0.1", "1"])
    # ax2.set_xlim(right=1.1)

    ax1.tick_params(axis='both', which='major', labelsize=fs_m)
    ax1.tick_params(axis='both', which='minor', labelsize=fs_m)
    ax2.tick_params(axis='both', which='major', labelsize=fs_m)
    ax2.tick_params(axis='both', which='minor', labelsize=fs_m)
    
    ax1.text(0.02, 0.015, 'A', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=ax1.transAxes)
    ax2.text(0.02, 0.015, 'B', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=ax2.transAxes)
    #fig.suptitle(r'$\alpha$=%.1f'%atm.alpha_cloud)
    fig.suptitle('$T_{surf}$=%.1f K, $T_{bot}$=%.1f K'% (atm.ts,atm.tmp[-1]))
    #plt.show()

    fig.savefig(filename, bbox_inches='tight')
    plt.close()

    return


####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    # Surface temperature & partial pressures
    T_surf                  = 400                           # K
    pH2O                    = p_sat('H2O',T_surf)           # Pa
    pCO2                    = 0.                            # Pa
    pH2                     = 0.                            # Pa
    pN2                     = 3e+5                          # Pa
    pCH4                    = 0.                            # Pa
    pO2                     = 0.                            # Pa
    pHe                     = 0.                            # Pa
    pNH3                    = 0.                            # Pa
    P_surf                  = pH2O + pCO2 + pH2 + pN2 + pCH4 + pO2 + pHe + pNH3  # Pa

    # Set fraction of condensate retained in column (0 = full rainout)
    alpha_cloud             = 0.0

    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass
    P_top         = 1.0                 # Pa
    
    # Volatile molar concentrations in the dictionary below are defined as fractions that must sum to one
    # The vanilla setting defines a water-saturated atmosphere with a 3 bar N2 background
    vol_list = { 
                  "H2O" : pH2O / P_surf,   
                  "CO2" : pCO2 / P_surf,
                  "H2"  : pH2  / P_surf,
                  "N2"  : pN2  / P_surf,
                  "CH4" : pCH4 / P_surf,
                  "O2"  : pO2  / P_surf,
                  "CO"  : pN2  / P_surf,
                  "He"  : pHe  / P_surf,
                  "NH3" : pNH3 / P_surf,
                }
    # Create atmosphere object
    atm                     = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, vol_mixing=vol_list)

    # Set fraction of condensate retained in column (0 = full rainout)
    atm.alpha_cloud         = alpha_cloud
    
    # Calculate moist adiabat + condensation
    atm                     = general_adiabat(atm)

    # Plot adiabat
    plot_adiabats(atm, filename="../output/general_adiabat.pdf")
    
