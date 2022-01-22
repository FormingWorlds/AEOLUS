
'''
Created 01/12/19

@authors: 
Ryan Boukrouche (RB)
Tim Lichtenberg (TL)
RJ Graham (RJ)
This file builds a self-consistent temperature profile from either Li, Ingersoll & Oyafuso 2018 or 
Graham, Boukrouche, Lichtenberg, & Pierrehumbert 2020, from the ground up using the Runge-Kutta 4 scheme 
from ClimateUtilities.py.
'''

import time
import numpy as np
import scipy.interpolate as spint
import math
# from scipy.integrate import odeint
# from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# import seaborn as sns
import copy
from matplotlib import cm
from cp_funcs import *
# Coupler-specific modules
try:
    from ClimateUtilities import *
    from atmosphere_column import atmos
    import phys
except:
    from atm_rad_conv.ClimateUtilities import *
    from atm_rad_conv.atmosphere_column import atmos
    import atm_rad_conv.phys as phys

# Color definitions: 
# https://www.codecademy.com/articles/seaborn-design-ii
# https://python-graph-gallery.com/python-colors/
# https://matplotlib.org/tutorials/colors/colormaps.html
# https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
no_colors   = 7
vol_colors = {
    "H2O"            : cm.get_cmap('PuBu', no_colors)(range(no_colors)),
    "CO2"            : cm.get_cmap("Reds", no_colors)(range(no_colors)),
    "H2"             : cm.get_cmap("Greens", no_colors)(range(no_colors)),
    "N2"             : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "O2"             : cm.get_cmap("Wistia", no_colors+2)(range(no_colors+2)),
    "CH4"            : cm.get_cmap("RdPu", no_colors)(range(no_colors)),
    "CO"             : cm.get_cmap("pink_r", no_colors)(range(no_colors)),
    "S"              : cm.get_cmap("YlOrBr", no_colors)(range(no_colors)),
    "He"             : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "NH3"            : cm.get_cmap("summer", no_colors)(range(no_colors)),
    "mixtures"       : cm.get_cmap("Set3", 9)(range(no_colors)),
    "H2O-CO2"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "CO2-H2O"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "H2O-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-H2O"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-CO"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "CO-H2"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "H2-CO2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "CO2-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "H2-CH4"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "CH4-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "H2-N2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "N2-H2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "CO2-N2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    "N2-CO2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    # "H2O"            : sns.color_palette("PuBu", no_colors),
    # "CO2"            : sns.color_palette("Reds", no_colors),
    # "H2"             : sns.color_palette("Greens", no_colors),
    # "N2"             : sns.color_palette("Purples", no_colors), # sns.cubehelix_palette(7)
    # "O2"             : sns.light_palette("darkturquoise", no_colors),
    # "CH4"            : sns.color_palette("RdPu", no_colors),
    # "CO"             : sns.light_palette("#731d1d", no_colors),
    # "S"              : sns.light_palette("#EBB434", no_colors),
    # "He"             : sns.color_palette("Greys", no_colors),
    # "NH3"            : sns.light_palette("teal", no_colors),
    # "mixtures"       : sns.color_palette("Set2", 11),
    # "H2O-CO2"        : sns.color_palette("Set2", 11)[9],
    # "CO2-H2O"        : sns.color_palette("Set2", 11)[9],
    # "H2O-H2"         : sns.color_palette("Set2", 11)[3],
    # "H2-H2O"         : sns.color_palette("Set2", 11)[3],
    # "H2-CO"          : sns.color_palette("Set2", 11)[7],
    # "CO-H2"          : sns.color_palette("Set2", 11)[7],
    # "H2-CO2"         : sns.color_palette("Set2", 11)[8],
    # "CO2-H2"         : sns.color_palette("Set2", 11)[8],
    # "H2-CH4"         : sns.color_palette("Set2", 11)[2],
    # "CH4-H2"         : sns.color_palette("Set2", 11)[2],
    # "H2-N2"          : sns.color_palette("Set2", 11)[4],
    # "N2-H2"          : sns.color_palette("Set2", 11)[4],
    # "CO2-N2"         : sns.color_palette("Set2", 11)[5],
    # "N2-CO2"         : sns.color_palette("Set2", 11)[5],
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

molar_mass      = {
          "H2O" : 0.01801528,           # kg mol−1
          "CO2" : 0.04401,              # kg mol−1
          "H2"  : 0.00201588,           # kg mol−1
          "CH4" : 0.01604,              # kg mol−1
          "CO"  : 0.02801,              # kg mol−1
          "N2"  : 0.028014,             # kg mol−1
          "O2"  : 0.031999,             # kg mol−1
          "SO2" : 0.064066,             # kg mol−1
          "H2S" : 0.0341,               # kg mol−1 
          "H"   : 0.001008,             # kg mol−1 
          "C"   : 0.012011,             # kg mol−1 
          "O"   : 0.015999,             # kg mol−1 
          "N"   : 0.014007,             # kg mol−1 
          "S"   : 0.03206,              # kg mol−1 
          "He"  : 0.0040026,            # kg mol−1 
          "NH3" : 0.017031,             # kg mol−1 
        }

def atm_z(atm, idx):

    # Mean temperature below
    #T_mean_down = np.mean( atm.tmp[:idx+1] )

    # # Average temperature weighted by amount of substance?
    # T_mean_down = np.average( atm.tmp[:idx+1], weights=)

    # # Use gas phase mean molar mass
    # atm.mu[idx] = atm.mu[idx]

    # print(atm.grav_z[idx], atm.mu[idx], atm.p[idx], atm.p[idx+1])
    dp = atm.p[idx+1]-atm.p[idx]
    
    mu_c = 0.
    if atm.xc[idx] > 0:
        for vol in atm.vol_list.keys():
            mu_c += atm.x_cond[vol][idx]*molar_mass[vol] / atm.xc[idx]
    
    # Integration
    #dz = - phys.R_gas * T_mean_down * np.log(atm.p[idx+1]/atm.p[idx]) / ( atm.mu[idx] * atm.grav_z[idx] )
    atm.rho[idx] =  atm.p[idx]/phys.R_gas/atm.tmp[idx] * ( atm.mu[idx] + mu_c*atm.alpha_cloud*atm.xc[idx]/(atm.xv[idx]+atm.xd[idx]) )
    dz =  -dp / atm.rho[idx]
    if dz<0:
        print(dz)
    # Next height
    atm.z[idx+1] = atm.z[idx] + dz

    # Next gravity
    atm.grav_z[idx+1] = atm.grav_s * ((atm.planet_radius)**2) / ((atm.planet_radius+atm.z[idx+1])**2)

    # print(idx, T_mean_down, dz, atm.z[idx], atm.z[idx+1], atm.grav_z[idx], atm.grav_z[idx+1])

    return atm

## Saturation vapor pressure [Pa] for given temperature T [K]. 
## Assuming the ideal gas law and a constant latent heat of vaporization. 
## Select the molecule of interest with the switch argument (a string).
def p_sat(switch,T): 
    
    # Define volatile
    if switch == 'H2O':
        e = phys.satvps_function(phys.water)
    if switch == 'CH4':
        e = phys.satvps_function(phys.methane)
    if switch == 'CO2':
        e = phys.satvps_function(phys.co2)
    if switch == 'CO':
        e = phys.satvps_function(phys.co)
    if switch == 'N2':
        e = phys.satvps_function(phys.n2)
    if switch == 'O2':
        e = phys.satvps_function(phys.o2)
    if switch == 'H2':
        e = phys.satvps_function(phys.h2)
    if switch == 'He':
        e = phys.satvps_function(phys.he)
    if switch == 'NH3':
        e = phys.satvps_function(phys.nh3)
    
    # Return saturation vapor pressure
    return e(T)

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
    
## Molar latent heat [J mol-1] for gas phase considered given a temperature T [K]. 
## Select the molecule of interest with the switch argument (a string).
def L_heat(switch, T, P):

    if switch == 'H2O':
        L_sublimation   = phys.H2O.L_sublimation
        L_vaporization  = phys.H2O.L_vaporization
        MolecularWeight = phys.H2O.MolecularWeight
        T_triple        = phys.H2O.TriplePointT
        T_crit          = phys.H2O.CriticalPointT
    
    if switch == 'CH4':
        L_sublimation   = phys.CH4.L_sublimation
        L_vaporization  = phys.CH4.L_vaporization
        MolecularWeight = phys.CH4.MolecularWeight
        T_triple        = phys.CH4.TriplePointT
        T_crit          = phys.CH4.CriticalPointT
        
    if switch == 'CO2':
        L_sublimation   = phys.CO2.L_sublimation
        L_vaporization  = phys.CO2.L_vaporization
        MolecularWeight = phys.CO2.MolecularWeight
        T_triple        = phys.CO2.TriplePointT
        T_crit          = phys.CO2.CriticalPointT
        
    if switch == 'CO':
        L_sublimation   = phys.CO.L_sublimation
        L_vaporization  = phys.CO.L_vaporization
        MolecularWeight = phys.CO.MolecularWeight
        T_triple        = phys.CO.TriplePointT
        T_crit          = phys.CO.CriticalPointT
            
    if switch == 'N2':
        L_sublimation   = phys.N2.L_sublimation
        L_vaporization  = phys.N2.L_vaporization
        MolecularWeight = phys.N2.MolecularWeight
        T_triple        = phys.N2.TriplePointT
        T_crit          = phys.N2.CriticalPointT
        
    if switch == 'O2':
        L_sublimation   = phys.O2.L_sublimation
        L_vaporization  = phys.O2.L_vaporization
        MolecularWeight = phys.O2.MolecularWeight
        T_triple        = phys.O2.TriplePointT
        T_crit          = phys.O2.CriticalPointT
    
    if switch == 'H2':
        L_sublimation   = phys.H2.L_vaporization # No H2 sublimation
        L_vaporization  = phys.H2.L_vaporization
        MolecularWeight = phys.H2.MolecularWeight
        T_triple        = phys.H2.TriplePointT
        T_crit          = phys.H2.CriticalPointT
        
    if switch == 'He':
        L_sublimation   = phys.He.L_vaporization  # No He sublimation
        L_vaporization  = phys.He.L_vaporization
        MolecularWeight = phys.He.MolecularWeight
        T_triple        = phys.He.TriplePointT
        T_crit          = phys.He.CriticalPointT
        
    if switch == 'NH3':
        L_sublimation   = phys.NH3.L_sublimation
        L_vaporization  = phys.NH3.L_vaporization
        MolecularWeight = phys.NH3.MolecularWeight
        T_triple        = phys.NH3.TriplePointT
        T_crit          = phys.NH3.CriticalPointT

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
    num     = (1. + (L/(Ra*T))*qsat)*Ra
    den     = cpa + (cpc + (L/(Rc*T) - 1.)*(L/T))*qsat
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


# dlnT/dlnP slope function 
def moist_slope(lnP, lnT, atm):
    
    # T instead lnT
    tmp = math.exp(lnT)

    # Find current atm index
    idx = int(np.amax(atm.ifatm))

    # Sum terms in equation
    num_sum     = 0.
    denom_sum1  = 0. 
    denom_sum2  = 0. 
    denom_sum3  = 0.

    # Calculate sums over volatiles
    for vol in atm.vol_list.keys(): 
            
        # Coefficients
        eta_vol     = atm.x_gas[vol][idx] / atm.xd[idx]
        beta_vol    = L_heat(vol, tmp, atm.p_vol[vol][idx]) / (phys.R_gas * tmp) 

        # Beta terms zero if below saturation vapor pressure
        if atm.p_vol[vol][idx] < p_sat(vol, tmp): beta_vol = 0.

        # Sum in numerator
        num_sum     += eta_vol * beta_vol

        # Sums in denominator
        denom_sum1  += eta_vol * (beta_vol**2.)
        denom_sum3  += eta_vol
                              
    # Sum 2 in denominator  
    denom_sum2  = num_sum ** 2.

    # Collect terms
    numerator   = 1. + num_sum
    denominator = (atm.cp[idx] / phys.R_gas) + (denom_sum1 + denom_sum2) / (1. + denom_sum3)

    # dlnT/dlnP
    dlnTdlnP = numerator / denominator

    # Moist adiabat slope
    return dlnTdlnP

def condensation( atm, idx, prs_reset):

    #print("BADGER1")
    # Temperature floor
    tmp = np.amax([atm.tmp[idx], 20.])
    #print("1. tmp = ", tmp)
    #print("BADGER2")
    
    
    if idx==0:
        #print("BADGER3")
        
        p_tot_pre_condensation = atm.p[idx]
        # Renormalize mixing ratios to ensure sum == 1
        if sum(atm.vol_list.values()) != 1.:
            vol_list_new = {}
            for vol in atm.vol_list.keys():
                vol_list_new[vol]   = atm.vol_list[vol] / sum(atm.vol_list.values())
            for vol in atm.vol_list.keys():
                atm.vol_list[vol] = vol_list_new[vol]
                #print("2. atm.vol_list[vol] = ", atm.vol_list[vol])
        #print("BADGER4")
    # Calculating the pressure an air parcel would have at this temperature IF 
    # the parcel moved adiabatically upward from the below level w/o condensation
    # e.g. following dry adiabat from previous level
    else:
        #print("BADGER5")
        # Calculating mean cp in layer below; multiplied  by xd+xv because this cp is renormalized by that sum during the calculation
        cp_local = atm.cp[idx-1] * ( atm.xd[idx-1] + atm.xv[idx-1] ) 
        #print("3. cp_local = ", cp_local)
        
        # Dry adiabatic pressure calculated from dP/dT = P/T * cp/R
        p_tot_pre_condensation = atm.p[idx-1] + ( tmp - atm.tmp[idx-1] ) * atm.p[idx-1] / atm.tmp[idx-1] * cp_local / phys.R_gas
        #print("4. p_tot_pre_condensation = ", p_tot_pre_condensation)
        #print("BADGER6")
    #print("BADGER7")
    
    
    # Partial pressures and molar mass: pre-condensation values
    for vol in atm.vol_list.keys():
        if idx == 0:
    
            #print("BADGER8")
            # Partial pressures scaled 
            # atm.p_vol[vol][idx] = atm.vol_list[vol] * atm.p[idx]
            atm.p_vol[vol][idx] = atm.vol_list[vol] * p_tot_pre_condensation
            # Mean gas phase molar mass
            atm.mu[idx]       += atm.vol_list[vol] * molar_mass[vol]
            #print("5. atm.p_vol[vol][idx] = ", atm.p_vol[vol][idx])
            #print("6. atm.mu[idx] = ", atm.mu[idx])
            #print("BADGER9")
        else:
            #print("BADGER10")
            # Sum up gaseous mixing ratios to normalize
            x_sum = 0
            
            for vol1 in atm.vol_list.keys():
                x_sum += atm.x_gas[vol1][idx-1]
                
            # Divide this volatile's mixing ratio by the sum of all mixing ratiosto acquire volume mixing ratio  
            x_vol = atm.x_gas[vol][idx-1] / x_sum
            
            # Scale the pre-condensation dry partial pressure by this volume mixing ratio
            atm.p_vol[vol][idx] = x_vol * p_tot_pre_condensation
            #print("7. atm.p_vol[vol][idx] = ", atm.p_vol[vol][idx])
            
            # Mean gas phase molar mass, pre-condensation
            atm.mu[idx] += x_vol * molar_mass[vol]
            #print("8. atm.mu[idx] = ", atm.mu[idx])
            #print("BADGER11")
    # Account for condensation
    
    # Total pressure of condensing phases
    p_cond_sum = 0
    
    # Molar mass before condensation
    mu_old = atm.mu[idx].copy()
    
    # Flag telling code to keep making new loops as long as new species condense
    condensation_flag = True
    
    # List of dry species
    dry_list = []
    
    # List of condensing species
    wet_list = []
    
    #print("BADGER12")
      
    while condensation_flag == True:
        
        # If new condensation does not occur, this flag will remain false and the loop
        # will end after the next iteration
        condensation_flag = False
        
        # Loop through species
        #print("BADGER13")
        for vol in atm.vol_list.keys():
            
            #print("BADGER14")
            # Condensation if p_i > p_sat
            if atm.p_vol[vol][idx] > p_sat(vol, tmp):
                #print("BADGER15")
                # Add condensing species to wet_list
                wet_list.append(vol)
                
                # Set condensation_flag to True -- this will force the loop
                # to iterate again to check if any new species are made to
                # condense by the condensation of this species
                condensation_flag = True
                
                # Set species partial pressure to p_sat
                atm.p_vol[vol][idx]  = p_sat(vol, tmp)
                #print("9. p_sat(vol, tmp) = ", p_sat(vol, tmp))
                #print("10. atm.p_vol[vol][idx] = ", atm.p_vol[vol][idx])
                
                # Add the species partial pressure to the condensing species sum
                p_cond_sum += atm.p_vol[vol][idx]
                
                # If the species is in dry_list from previous iteration of while loop, remove it
                if vol in dry_list:
                    
                    dry_list.remove(vol)
                #print("BADGER16")
            else:
                #print("BADGER17")
                # Add a non-condensing species to the list of dry species
                # if it's present (>0) and if it's not in a list already
                if atm.vol_list[vol] > 0. and vol not in dry_list and vol not in wet_list:
                    dry_list.append(vol)
                #print("BADGER18")
    
    #print(dry_list)
    # Calculate the total partial pressure of dry species
    p_dry_tot = atm.p[idx] - p_cond_sum
    
    # Sum of dry volatile fractions from vol_list
    dry_frac_sum = 0
    #print("BADGER19")
    if idx == 0:
        # Calculate dry_frac_sum
        for vol in dry_list:
            dry_frac_sum += atm.vol_list[vol]
    else:
        for vol in dry_list:
            dry_frac_sum += atm.x_gas[vol][idx-1]
        #print("BADGER20")
    # Calculate individual dry partial pressures by multiplying the total dry partial
    # pressure by the mixing ratio of each dry species w.r.t. the rest of the dry species
    if idx == 0:
        for vol in dry_list:
            #print(vol)
            atm.p_vol[vol][idx] = p_dry_tot * ( atm.vol_list[vol] / dry_frac_sum )
            #print("BADGER21")
            #print("11. atm.p_vol[vol][idx] = ", atm.p_vol[vol][idx])
    else:
        for vol in dry_list:
            atm.p_vol[vol][idx] = p_dry_tot * ( atm.x_gas[vol][idx-1] / dry_frac_sum )
    # Reset mean molar mass
    atm.mu[idx]   = 0.
    p_vol_sum       = 0.
    #print("BADGER22")
    # Calculate new mean molar mass
    for vol in atm.vol_list.keys():
        atm.mu[idx]   += molar_mass[vol] * atm.p_vol[vol][idx]
    atm.mu[idx] /= atm.p[idx]
    #print("BADGER23")
    #print("12. atm.mu[idx] = ", atm.mu[idx])
    
    # Update mixing ratios
    for vol in atm.vol_list.keys():

        # Mean molar mass
        p_vol_sum       += atm.p_vol[vol][idx]
        #print("13. p_sat(vol, tmp) = ", p_sat(vol, tmp))
        #print("14. atm.p_vol[vol][idx] = ", atm.p_vol[vol][idx])

        # Condensate phase
        if atm.p_vol[vol][idx] < p_sat(vol, tmp):
            atm.x_cond[vol][idx] = 0.
        else:
            #atm.x_cond[vol][idx] = atm.vol_list[vol] - ( atm.p_vol[vol][idx] / atm.p[idx] )
            if idx == 0:
                atm.x_cond[vol][idx] = atm.vol_list[vol] - ( mu_old / atm.mu[idx] * atm.p_vol[vol][idx] / p_tot_pre_condensation )
            else:
                atm.x_cond[vol][idx] = atm.x_cond[vol][idx-1] + atm.x_gas[vol][idx-1] - ( 1-atm.xc[idx-1] ) * ( mu_old / atm.mu[idx] * atm.p_vol[vol][idx] / p_tot_pre_condensation )
            
        #print("BADGER24")
        #atm.x_cond[vol][idx] *= atm.alpha_cloud
        # Add to molar concentration of total condensed phase
        atm.xc[idx]     += atm.x_cond[vol][idx]
        #print("15. atm.x_cond[vol][idx] = ", atm.x_cond[vol][idx])
        
        # Gas phase molar concentration
        #atm.x_gas[vol][idx] = atm.p_vol[vol][idx] / atm.p[idx]
        if idx == 0:
            atm.x_gas[vol][idx] =  mu_old / atm.mu[idx] * atm.p_vol[vol][idx] /  p_tot_pre_condensation
        else:
            atm.x_gas[vol][idx] =  ( 1-atm.xc[idx-1] ) * mu_old / atm.mu[idx] * atm.p_vol[vol][idx] /  p_tot_pre_condensation
        #print("BADGER25")
        # Add to molar concentration of total gas (dry or moist) phase
        # ! REVISIT ! keeping xd == 0 leads to a bug, why?
        if atm.p_vol[vol][idx] < p_sat(vol, tmp):
            atm.xd[idx]          += atm.x_gas[vol][idx]
        else:
            atm.xv[idx]          += atm.x_gas[vol][idx]
        #print("BADGER26")
        # Mean cp of both gas phase and retained condensates
        atm.cp[idx]    += atm.x_gas[vol][idx] * cpv(vol, atm.tmp[idx]) + atm.x_cond[vol][idx] * cp_cond(vol, atm.tmp[idx]) * atm.alpha_cloud

    # Renormalize cp w/ molar concentration (= cp_hat, Eq. 1, Li+2018)
    atm.cp[idx]  = atm.cp[idx] / ( atm.xd[idx] + atm.xv[idx] )
    #print("16. atm.cp[idx] = ", atm.cp[idx])

    # Dry concentration floor
    atm.xd[idx]  = np.amax([atm.xd[idx], 1e-10])
    #print("BADGER27")
    #print("17. atm.xd[idx] = ", atm.xd[idx])
    
    return atm

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
    Tsurf = atm.ts
    alpha = atm.alpha_cloud
    for vol in atm.vol_list.keys():
        if atm.vol_list[vol] * atm.ps > p_sat(vol, atm.ts):
            new_psurf += p_sat(vol,atm.ts)
            new_p_vol[vol] = p_sat(vol,atm.ts)
        else:
            new_psurf += atm.vol_list[vol] * atm.ps
            new_p_vol[vol] = atm.vol_list[vol] * atm.ps
            
    if new_psurf != atm.ps:
        for vol in atm.vol_list.keys():
            atm.vol_list[vol] = new_p_vol[vol] / new_psurf
        atm = atmos(Tsurf, new_psurf, atm.vol_list)
        atm.alpha_cloud = alpha
        
    
    ### Initialization
    
    # Initialize the tuple solution
    moist_tuple     = [] #[tuple([np.log(atm.ps), atm.ts])] 
   
    # Negative increment to go from ps to ptop < ps       
    step            = -.001

    # Integration counter
    idx             = 0  
    
    # Surface condensation
    atm             = condensation(atm, idx, prs_reset=False)

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
        atm             = condensation(atm, idx, prs_reset=False)

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

# Plotting
def plot_adiabats(atm):

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
            Psat_array = [ p_sat(vol, T) for T in T_sat_array ]
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
    fig.suptitle(r'$\alpha$=%.1f'%atm.alpha_cloud)
    #plt.show()

    plt.savefig('./output/general_adiabat.pdf', bbox_inches='tight')
    plt.close(fig)  

    return


####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    # Surface pressure & temperature
    P_surf                  = 20e+5       # Pa
    T_surf                  = 300         # K

    # Volatile molar concentrations: ! must sum to one !
    vol_list = { 
                  "H2O" : 0,    # 300e+5/P_surf --> specific p_surf
                  "CO2" : 1,     # 100e+5/P_surf
                  "H2"  : .0, 
                  "N2"  : .0,     # 1e+5/P_surf
                  "CH4" : .0, 
                  "O2"  : .0, 
                  "CO"  : .0, 
                  "He"  : .0,
                  "NH3" : .0, 
                }

    # Create atmosphere object
    atm                     = atmos(T_surf, P_surf, vol_list)
    
    # Set fraction of condensate retained in column
    atm.alpha_cloud         = 0.
    
    # Calculate moist adiabat + condensation
    atm                     = general_adiabat(atm)

    # Plot adiabat
    plot_adiabats(atm)
    
