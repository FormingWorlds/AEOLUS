'''
Created 28/01/19

@authors:
Mark Hammond (MH)
Tim Lichtenberg (TL)

SOCRATES radiative-convective model
'''

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
# import seaborn as sns
import copy
import pathlib
import pickle as pkl
import json
import glob, re, os

try:
    import phys
    import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
    import SocRadModel
    from atmosphere_column import atmos
except:
    import atm_rad_conv.phys as phys
    import atm_rad_conv.GeneralAdiabat as ga
    import atm_rad_conv.SocRadModel as SocRadModel
    from atm_rad_conv.atmosphere_column import atmos

# String sorting not based on natsorted package
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def surf_Planck_nu(atm):
    h   = 6.63e-34
    c   = 3.0e8
    kb  = 1.38e-23
    B   = np.zeros(len(atm.band_centres))
    c1  = 1.191042e-5
    c2  = 1.4387752
    for i in range(len(atm.band_centres)):
        nu      = atm.band_centres[i]
        B[i]    = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))
    B   = (1.-atm.albedo_s) * np.pi * B * atm.band_widths/1000.0
    return B

def RadConvEqm(dirs, time, atm, loop_counter, COUPLER_options, standalone, cp_dry, moist_timestep, trpp, rscatter):

    ### Moist/general adiabat
    atm_moist = compute_moist_adiabat(atm, dirs, standalone, trpp, rscatter)

    ### Dry adiabat
    if cp_dry == True:

        # Compute dry adiabat  w/ timestepping
        atm_dry   = compute_dry_adiabat(atm, dirs, standalone, rscatter)

        if standalone == True:
            print("Net, OLR => moist:", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print("| dry:", str(round(atm_dry.net_flux[0], 3)), str(round(atm_dry.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print()
    else: atm_dry = {}

    ### Moist/general adiabat w/ timestepping
    if moist_timestep == True:
        atm_moist_timestep = compute_moist_adiabat_timestep(atm, dirs, standalone, rscatter)

        if standalone == True:
            print("Net, OLR => moist:", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print("| moist timestepped:", str(round(atm_moist_timestep.net_flux[0], 3)), str(round(atm_moist_timestep.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print()
    else: atm_moist_timestep = {}

    # Plot
    if standalone == True:
        plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs)
        # Save to disk
        with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
            pkl.dump(atm_moist, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    return atm_dry, atm_moist, atm_moist_timestep

# Dry adiabat profile
def dry_adiabat_atm(atm):

    R_universal = 8.31446261815324              # Universal gas constant, J.K-1.mol-1

    # Calculate cp from molar concentrations
    cp_mix = 0.
    for vol in atm.vol_list.keys():
        cp_mix += atm.vol_list[vol] * ga.cpv(vol, atm.ts)

    # Calculate dry adiabat slope
    atm.Rcp = R_universal / cp_mix

    # Calculate dry adiabat temperature profile for staggered nodes (from ground up)
    for idx, prsl in enumerate(atm.pl):
        atm.tmpl[idx] = atm.ts * ( prsl / atm.ps ) ** ( atm.Rcp )

    # Interpolate temperature from staggered nodes
    atm.tmp = np.interp(atm.p, atm.pl, atm.tmpl)

    return atm
  
# Dry convective adjustment - changes tmp directly: bad
def DryAdj_old(atm):

    T   = atm.tmp
    p   = atm.p
    
    # Rcp is global
    # Downward pass
    for i in range(len(T)-1):
        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        
        # Adiabat slope
        pfact = (p1/p2)**atm.Rcp
        
        # If slope is shallower than adiabat (unstable), adjust to adiabat
        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) # Equal layer masses
                              # Not quite compatible with how
                              # heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.tmp[i]   = T1
            atm.tmp[i+1] = T2
    
    # Upward pass
    for i in range(len(T)-2,-1,-1):

        T1,p1 = T[i],p[i]
        T2,p2 = T[i+1],p[i+1]
        pfact = (p1/p2)**atm.Rcp

        if T1 < T2*pfact:
            Tbar = .5*(T1+T2) # Equal layer masses
                              # Not quite compatible with how
                              # heating is computed from flux
            T2 = 2.*Tbar/(1.+pfact)
            T1 = T2*pfact
            atm.tmp[i]   = T1
            atm.tmp[i+1] = T2 

    return atm      

# Dry convective adjustment - returns the tendency: good. Also conserves dry enthalpy.
def DryAdj(atm, itermax):

    Tl    = atm.tmpl
    pl    = atm.pl
    pe    = atm.p
    Tl_cc = Tl
    nlay  = len(Tl_cc)

    d_p   = np.ones(len(pl))

    small = 1.0e-6

    for i in range(1, nlay+1):
        d_p[i] = pe[i+1] - pe[i]

    for iteration in range(1, itermax+1):
        did_adj = False

        # Downward pass
        for i in range(1, nlay): # from layer 1 to layer nlay-1

            pfact = (pl[i]/pl[i+1])**atm.Rcp
            if (Tl_cc[i] < (Tl_cc[i+1]*pfact - small)):
                Tbar = (d_p[i]*Tl_cc[i] + d_p[i+1]*Tl_cc[i+1]) / (d_p[i] + d_p[i+1])
                Tl_cc[i+1] = (d_p[i] + d_p[i+1])*Tbar / (d_p[i+1]+pfact*d_p[i])
                Tl_cc[i] = Tl_cc[i+1] * pfact
                did_adj = True

        # Upward pass
        for i in range(nlay-1, 0, -1): # from layer nlay-1 to layer 1

            pfact = (pl[i]/pl[i+1])**atm.Rcp
            if (Tl_cc[i] < (Tl_cc[i+1]*pfact - small)):
                Tbar = (d_p[i]*Tl_cc[i] + d_p[i+1]*Tl_cc[i+1]) / (d_p[i] + d_p[i+1])
                Tl_cc[i+1] = (d_p[i] + d_p[i+1])*Tbar / (d_p[i+1]+pfact*d_p[i])
                Tl_cc[i] = Tl_cc[i+1] * pfact
                did_adj = True

        # If no adjustment required, exit the loop
        if (did_adj == False):
            break

    # Change in temperature is Tl_cc - Tl
    # adjust on timescale of 1 timestep (i.e. instant correction across one timestep)
    dT_conv_dry = (Tl_cc - Tl)/atm.dt
    return dT_conv_dry


# # Moist convective adjustment
# def MoistAdj(atm, dT):

#     # Apply heating
#     tmp_heated = atm.tmp + dT

#     # Reset to moist adiabat if convectively unstable
#     for idx, tmp_heated in enumerate(tmp_heated):
#         atm.tmp[idx] = np.amax([tmp_heated, atm.tmp[idx]])

#     return atm 

def MoistAdj(atm, itermax):

    Tl    = atm.tmpl
    pl    = atm.pl
    pe    = atm.p
    Tl_cc = Tl
    Tl_cc_H2O = Tl
    Tl_cc_CH4 = Tl
    Tl_cc_CO2 = Tl
    Tl_cc_CO  = Tl
    Tl_cc_N2  = Tl
    Tl_cc_O2  = Tl
    Tl_cc_H2  = Tl
    Tl_cc_He  = Tl
    Tl_cc_NH3 = Tl

    nlay  = len(Tl_cc)

    for iteration in range(1, itermax+1):
        did_adj = False
        #------------------- H2O -------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('H2O', pl[i])):
                Tl_cc_H2O[i] = ga.Tdew('H2O', pl[i]) # temperature stays the same during the phase change
                did_adj = True
        
        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('H2O', pl[i])):
                Tl_cc_H2O[i] = ga.Tdew('H2O', pl[i])
                did_adj = True

        #------------------- CH4 -------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('CH4', pl[i])):
                Tl_cc_CH4[i] = ga.Tdew('CH4', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('CH4', pl[i])):
                Tl_cc_CH4[i] = ga.Tdew('CH4', pl[i])
                did_adj = True

        #------------------- CO2 -------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('CO2', pl[i])):
                Tl_cc_CO2[i] = ga.Tdew('CO2', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('CO2', pl[i])):
                Tl_cc_CO2[i] = ga.Tdew('CO2', pl[i])
                did_adj = True

        #------------------- CO --------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('CO', pl[i])):
                Tl_cc_CO[i] = ga.Tdew('CO', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('CO', pl[i])):
                Tl_cc_CO[i] = ga.Tdew('CO', pl[i])
                did_adj = True

        #------------------- N2 --------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('N2', pl[i])):
                Tl_cc_N2[i] = ga.Tdew('N2', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('N2', pl[i])):
                Tl_cc_N2[i] = ga.Tdew('N2', pl[i])
                did_adj = True

        #------------------- O2 --------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('O2', pl[i])):
                Tl_cc_O2[i] = ga.Tdew('O2', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('O2', pl[i])):
                Tl_cc_O2[i] = ga.Tdew('O2', pl[i])
                did_adj = True

        #------------------- H2 --------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('H2', pl[i])):
                Tl_cc_H2[i] = ga.Tdew('H2', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('H2', pl[i])):
                Tl_cc_H2[i] = ga.Tdew('H2', pl[i])
                did_adj = True

        #------------------- He --------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('He', pl[i])):
                Tl_cc_He[i] = ga.Tdew('He', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('He', pl[i])):
                Tl_cc_He[i] = ga.Tdew('He', pl[i])
                did_adj = True

        #------------------- NH3 -------------------
        #Downward pass
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            if (Tl_cc[i] < ga.Tdew('NH3', pl[i])):
                Tl_cc_NH3[i] = ga.Tdew('NH3', pl[i]) # temperature stays the same during the phase change
                did_adj = True

        # Upward pass
        for i in range(nlay, 0, -1): # from layer nlay to layer 1
            if (Tl_cc[i] < ga.Tdew('NH3', pl[i])):
                Tl_cc_NH3[i] = ga.Tdew('NH3', pl[i])
                did_adj = True

        #----------------- Minimum -----------------
        for i in range(1, nlay+1): # from layer 1 to layer nlay
            Tl_cc[i] = np.min([Tl_cc_H2O[i], Tl_cc_CH4[i], Tl_cc_CO2[i], Tl_cc_CO[i], Tl_cc_N2[i], Tl_cc_O2[i], Tl_cc_H2[i], Tl_cc_He[i], Tl_cc_NH3[i]])

        # If no adjustment required, exit the loop
        if (did_adj == False):
            break

    # Change in temperature is Tmid_cc - Tmid
    dT_conv_moist = (Tl_cc - Tl)/atm.dt
    return dT_conv_moist


def plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs):

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13,10))
    # sns.set_style("ticks")
    # sns.despine()

    # Line settings
    col_idx  = 3
    col_vol1 = "H2O"
    col_vol2 = "N2"
    col_vol3 = "H2"
    col_vol4 = "O2"

    # Temperature vs. pressure
    ax1.semilogy(atm_moist.tmp,atm_moist.p, color=ga.vol_colors[col_vol1][col_idx+1], ls="-", label=r'Moist adiabat')
    if cp_dry == True: ax1.semilogy(atm_dry.tmp,atm_dry.p, color=ga.vol_colors[col_vol3][col_idx+1], ls="-", label=r'Dry adiabat')
    ax1.legend()
    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ (K)')
    ax1.set_ylabel(r'Pressure $P$ (Pa)')
    # ax1.set_ylim(bottom=atm_moist.ps*1.01)
    ax1.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)

    # Print active species
    active_species = r""
    for vol in atm_moist.vol_list:
        if atm_moist.vol_list[vol] > 1e-5:
            active_species = active_species + ga.vol_latex[vol] + ", "
    active_species = active_species[:-2]
    ax1.text(0.02, 0.02, r"Active species: "+active_species, va="bottom", ha="left", fontsize=10, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.5, boxstyle='round', pad=0.1), color=ga.vol_colors["black_1"] )

    # Fluxes vs. pressure

    # Zero line
    ax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
    # ax2.axvline(-1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)
    # ax2.axvline(1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)

    # SW down / instellation flux
    if cp_dry == True: ax2.semilogy(atm_dry.SW_flux_down*(-1),atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx], ls=":")
    ax2.semilogy(atm_moist.SW_flux_down*(-1),atm_moist.pl, color=ga.vol_colors[col_vol2][col_idx], ls=":", label=r'$F_{\odot}^{\downarrow}$')

    # LW down / thermal downward flux
    if cp_dry == True: ax2.semilogy(atm_dry.LW_flux_down*(-1),atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx], ls="--")
    ax2.semilogy(atm_moist.LW_flux_down*(-1),atm_moist.pl, color=ga.vol_colors[col_vol2][col_idx+1], ls="--", label=r'$F_\mathrm{t}^{\downarrow}$')
    # ls=(0, (3, 1, 1, 1))
    
    # Thermal upward flux, total
    if cp_dry == True: ax2.semilogy(atm_dry.flux_up_total,atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx], ls="--")
    ax2.semilogy(atm_moist.flux_up_total,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx], ls="--", label=r'$F_\mathrm{t}^{\uparrow}$')

    # Net flux
    if cp_dry == True: ax2.semilogy(atm_dry.net_flux,atm_dry.pl, color=ga.vol_colors[col_vol3][6], ls="-", lw=2)
    ax2.semilogy(atm_moist.net_flux,atm_moist.pl, color=ga.vol_colors[col_vol1][6], ls="-", lw=2, label=r'$F_\mathrm{net}^{\uparrow}$')

    # # SW up
    # if cp_dry == True: ax2.semilogy(atm_dry.SW_flux_up,atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx], ls=":")
    # ax2.semilogy(atm_moist.SW_flux_up,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx], ls=":", label=r'$F_\mathrm{SW}^{\uparrow}$')
    
    # # LW up
    # if cp_dry == True: ax2.semilogy(atm_dry.LW_flux_up,atm_dry.pl, color=ga.vol_colors[col_vol3][col_idx], ls=(0, (5, 1)))
    # ax2.semilogy(atm_moist.LW_flux_up,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx], ls=(0, (5, 1)), label=r'$F_\mathrm{LW}^{\uparrow}$')
    
    ax2.legend(ncol=6, fontsize=10, loc=3)
    ax2.invert_yaxis()
    ax2.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
    ax2.set_xlabel(r'Outgoing flux $F^{\uparrow}$ (W m$^{-2}$)')
    ax2.set_ylabel(r'Pressure $P$ (Pa)')
    ax2.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)

    # Wavenumber vs. OLR
    ax3.plot(atm_moist.band_centres, surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label=str(round(atm_moist.ts))+' K blackbody')
    if cp_dry == True: ax3.plot(atm_dry.band_centres, atm_dry.net_spectral_flux[:,0]/atm_dry.band_widths, color=ga.vol_colors[col_vol3][col_idx])
    ax3.plot(atm_moist.band_centres, atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
    ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
    ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
    ax3.legend(loc=1)
    ymax_plot = 1.2*np.max(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths)
    ax3.set_ylim(bottom=0, top=ymax_plot)
    ax3.set_xlim(left=0, right=np.max(np.where(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
    # ax3.set_xlim(left=0, right=30000)
    

    # Heating versus pressure
    ax4.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

    if cp_dry == True: 
        ax4.plot(atm_dry.LW_heating, atm_dry.p, ls="--", color=ga.vol_colors[col_vol3][col_idx+1])
        ax4.plot(atm_dry.net_heating, atm_dry.p, lw=2, color=ga.vol_colors[col_vol3][col_idx+1])
        ax4.plot(atm_dry.SW_heating, atm_dry.p, ls=":", color=ga.vol_colors[col_vol3][col_idx+1])

    ax4.plot(atm_moist.LW_heating, atm_moist.p, ls="--", color=ga.vol_colors[col_vol1][col_idx+1], label=r'LW')
    ax4.plot(atm_moist.net_heating, atm_moist.p, lw=2, color=ga.vol_colors[col_vol1][col_idx+1], label=r'Net')
    ax4.plot(atm_moist.SW_heating, atm_moist.p, ls=":", color=ga.vol_colors[col_vol1][col_idx+1], label=r'SW')

    # Plot tropopause
    trpp_idx = int(atm_moist.trpp[0])
    if trpp_idx > 0:
        ax4.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol2][col_idx], lw=1.0, ls="-.", label=r'Tropopause')
    
    ax4.invert_yaxis()
    ax4.legend(ncol=4, fontsize=10, loc=3)
    ax4.set_ylabel(r'Pressure $P$ (Pa)')
    ax4.set_xlabel(r'Heating (K/day)')
    # ax4.set_xscale("log")
    ax4.set_yscale("log") 
    x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
    x_minmax = np.max([ 20, x_minmax ])
    if not math.isnan(x_minmax):
        ax4.set_xlim(left=-x_minmax, right=x_minmax)
    ax4.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)
    # ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
    # ax4.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
    # ax4.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    
    # # Wavelength versus OLR log plot
    # OLR_cm_moist = atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths
    # wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ]          # microns
    # OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ]                    # microns
    # if cp_dry == True: 
    #     OLR_cm_dry = atm_dry.LW_spectral_flux_up[:,0]/atm_dry.band_widths
    #     wavelength_dry  = [ 1e+4/i for i in atm_dry.band_centres ]              # microns
    #     OLR_micron_dry  = [ 1e+4*i for i in OLR_cm_dry ]                        # microns
    #     ax4.plot(wavelength_dry, OLR_micron_dry, color=ga.vol_colors[col_vol3][col_idx+1])
    
    # ax4.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol1][col_idx+1])
    # ax4.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')
    # ax4.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
    # ax4.set_xscale("log")
    # ax4.set_yscale("log") 
    # ax4.set_xlim(left=0.1, right=100)
    # ax4.set_ylim(bottom=1e-20, top=1e5)
    # # ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
    # ax4.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
    # ax4.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])

    plt.savefig(dirs["output"]+"/"+"TP_"+str(round(time["planet"]))+'.pdf', bbox_inches="tight")
    plt.close(fig)

    # with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
    #     pkl.dump(atm, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    # # Save atm object to .json file
    # json_atm = json.dumps(atm.__dict__)
    # with open(dirs["output"]+"/"+str(int(time_current))+"_atm.json", "wb") as atm_file:
    #     json.dump(json_atm, atm_file)

# Time integration for n steps
def compute_dry_adiabat(atm, dirs, standalone, rscatter):

    # Dry adiabat settings 
    rad_steps   = 1000  # Maximum number of radiation steps
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

    # Time stepping
    for i in range(0, rad_steps):

        # atm.toa_heating = 0

        # Compute radiation, midpoint method time stepping
        try:
            atm_dry         = SocRadModel.radCompSoc(atm_dry, dirs, recalc=False, calc_cf=False, rscatter=rscatter)
            
            dT_dry          = atm_dry.net_heating * atm_dry.dt

            # Dry convective adjustment
            dT_conv_dry  = DryAdj(atm_dry, conv_steps)
            # Apply temperature tendency due to dry convection
            dT_dry += dT_conv_dry

            # Limit the temperature change per step
            dT_dry          = np.where(dT_dry > dT_max, dT_max, dT_dry)
            dT_dry          = np.where(dT_dry < -dT_max, -dT_max, dT_dry)

            # Apply heating
            atm_dry.tmp     += dT_dry

            # # Do the surface balance
            # kturb       = .1
            # atm.tmp[-1] += -atm.dt * kturb * (atm.tmp[-1] - atm.ts)
            
            # Dry convective adjustment (old, changed tmp directly)
            #for iadj in range(conv_steps):
                #atm_dry     = DryAdj(atm_dry)

            # Temperature floor to prevent SOCRATES crash
            if np.min(atm_dry.tmp) < T_floor:
                atm_dry.tmp = np.where(atm_dry.tmp < T_floor, T_floor, atm_dry.tmp)

            # Convergence criteria
            dTglobal_dry    = abs(round(np.max(atm_dry.tmp-PrevTemp_dry[:]), 4))
            dTtop_dry       = abs(round(atm_dry.tmp[0]-atm_dry.tmp[1], 4))

            # Break criteria
            dOLR_dry        = abs(round(atm_dry.LW_flux_up[0]-PrevOLR_dry, 6))
            dbreak_dry      = (0.01*(5.67e-8*atm_dry.ts**4)**0.5)

            # Inform during runtime
            if i % 2 == 1 and standalone == True:
                print("Dry adjustment step", i+1, end=": ")
                print("OLR = " + str(atm_dry.LW_flux_up[0]) + " W/m^2,", "dT_max = " + str(dTglobal_dry) + " K, dT_top = " + str(dTtop_dry) + " K, dOLR = " + str(dOLR_dry) + " W/m^2,")

            # Reduce timestep if heating is not converging
            if dTglobal_dry < 0.05 or dTtop_dry > dT_max:
                atm_dry.dt  = atm_dry.dt*0.99
                if standalone == True:
                    print("Dry adiabat not converging -> dt_new =", round(atm_dry.dt,5), "days")

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


def compute_moist_adiabat(atm, dirs, standalone, trpp, rscatter): # no time-stepping

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(copy.deepcopy(atm))

    # Run SOCRATES
    atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=False, rscatter=rscatter)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    if trpp == True:
        
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=False, rscatter=rscatter)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist

def compute_moist_adiabat_timestep(atm, dirs, standalone, rscatter): # with time-stepping

    # Moist adiabat settings
    rad_steps   = 10  # Maximum number of radiation steps
    conv_steps  = 30   # Number of convective adjustment steps (per radiation step)
    dT_max      = 20.  # K, Maximum temperature change per radiation step
    T_floor     = 10.  # K, Temperature floor to prevent SOCRATES crash

    # Build general adiabat structure before the time-stepping
    atm_moist = ga.general_adiabat(copy.deepcopy(atm))

    # Initialise previous OLR and TOA heating to zero
    PrevOLR_dry         = 0.
    PrevMaxHeat_dry     = 0.
    PrevTemp_dry        = atm.tmp * 0.

    # Time stepping
    for i in range(0, rad_steps):

        # Compute radiation, midpoint method time stepping

        try:
            # Run SOCRATES
            atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=False, rscatter=rscatter)

            dT_moist  = atm_moist.net_heating * atm_moist.dt

            # Dry convective adjustment
            dT_conv_dry  = DryAdj(atm_dry, conv_steps)
            # Apply temperature tendency due to dry convection
            dT_moist += dT_conv_dry

            # Moist convective adjustment
            dT_conv_moist  = MoistAdj(atm_moist, conv_steps)

            # Apply temperature tendency due to moist convection
            dT_moist += dT_conv_moist

            # Limit the temperature change per step
            dT_moist          = np.where(dT_moist > dT_max, dT_max, dT_moist)
            dT_moist          = np.where(dT_moist < -dT_max, -dT_max, dT_moist)

            # Apply heating
            atm_moist.tmp     += dT_moist

            # The stratosphere scheme is redundant here, a stratosphere will form on its own with timestepping and stellar heating
            # -----------------------------------------------------------------------------------------
            #if standalone == True:
            #    print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

            #if trpp == True:
        
            #    # Find tropopause index
            #    atm_moist = find_tropopause(atm_moist)

                # Reset stratosphere temperature and abundance levels
            #    atm_moist = set_stratosphere(atm_moist)

                # Recalculate fluxes w/ new atmosphere structure
            #    atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=False, rscatter=rscatter)

            #    if standalone == True:
            #        print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")
            # -----------------------------------------------------------------------------------------

            # Temperature floor to prevent SOCRATES crash
            if np.min(atm_moist.tmp) < T_floor:
                atm_moist.tmp = np.where(atm_moist.tmp < T_floor, T_floor, atm_moist.tmp)

            # Convergence criteria
            dTglobal_moist    = abs(round(np.max(atm_moist.tmp-PrevTemp_moist[:]), 4))
            dTtop_moist       = abs(round(atm_moist.tmp[0]-atm_moist.tmp[1], 4))

            # Break criteria
            dOLR_moist        = abs(round(atm_moist.LW_flux_up[0]-PrevOLR_moist, 6))
            dbreak_moist      = (0.01*(5.67e-8*atm_moist.ts**4)**0.5)

            # Inform during runtime
            if i % 2 == 1 and standalone == True:
                print("Moist adjustment step", i+1, end=": ")
                print("OLR = " + str(atm_moist.LW_flux_up[0]) + " W/m^2,", "dT_max = " + str(dTglobal_moist) + " K, dT_top = " + str(dTtop_moist) + " K, dOLR = " + str(dOLR_moist) + " W/m^2,")

            # Reduce timestep if heating is not converging
            if dTglobal_moist < 0.05 or dTtop_moist > dT_max:
                atm_moist.dt  = atm_moist.dt*0.99
                if standalone == True:
                    print("Moist adiabat not converging -> dt_new =", round(atm_moist.dt,5), "days")

            # Sensitivity break condition
            if (dOLR_moist < dbreak_moist) and i > 5:
                if standalone == True:
                    print("Timestepping break ->", end=" ")
                    print("dOLR/step =", dOLR_moist, "W/m^2, dTglobal_moist =", dTglobal_moist)
                break    # break here

        except:
            if standalone == True:
                print("Socrates cannot be executed properly, T profile:", atm_moist.tmp)
            break    # break here

        PrevOLR_moist       = atm_moist.LW_flux_up[0]
        PrevMaxHeat_moist   = abs(np.max(atm_moist.net_heating))
        PrevTemp_moist[:]   = atm_moist.tmp[:]

    return atm_moist


def find_tropopause(atm_moist):

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
    
        atm_moist.trpp[0] = trpp_idx                  # index
        atm_moist.trpp[1] = atm_moist.pl[trpp_idx]    # pressure 
        atm_moist.trpp[2] = atm_moist.tmpl[trpp_idx]  # temperature

    return atm_moist


def set_stratosphere(atm):

    trpp_idx = int(atm.trpp[0])
    trpp_prs = atm.trpp[1]
    trpp_tmp = atm.trpp[2]

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

def InterpolateStellarLuminosity(star_mass, time, mean_distance, albedo, Sfrac):

    # Constants
    L_sun                   = 3.828e+26        # W, IAU definition
    AU                      = 1.495978707e+11  # m

    # File name
    fname = str(pathlib.Path(__file__).parent.absolute())+"/luminosity_tracks/Lum_m"+str(star_mass)+".txt"

    # If file exists, just interpolate needed time
    if os.path.isfile(fname):

        luminosity_df           = pd.read_csv(fname)
        star_age                = (time["star"])/1e+6   # Myr
        ages                    = luminosity_df["age"]*1e+3         # Myr
        luminosities            = luminosity_df["lum"]              # L_sol

        # Interpolate luminosity for current time
        interpolate_luminosity  = interpolate.interp1d(ages, luminosities)
        interpolated_luminosity = interpolate_luminosity([star_age])
        interpolated_luminosity = interpolated_luminosity[0]

    # Else: interpolate from 2D grid
    else:

        # Find all luminosity tracks
        lum_tracks = natural_sort(glob.glob(str(pathlib.Path(__file__).parent.absolute())+"/luminosity_tracks/"+"Lum_m*.txt"))

        # Define data arrays for interpolation later on
        xy_age_mass = []
        z_lum       = []

        # Fill the arrays
        for lum_track in lum_tracks:

            # Read the specific data file
            luminosity_df    = pd.read_csv(lum_track)

            # Cut the string to get the mass of the star
            star_mass_lst    = float(lum_track[-7:-4])

            # Read out age and luminosity
            age_list         = list(luminosity_df["age"]*1e+3)
            luminosity_list  = list(luminosity_df["lum"])

            mass_list        = np.ones(len(age_list))*star_mass_lst

            # Fill the arrays
            zip_list = list(zip(age_list, mass_list))
            xy_age_mass.extend(zip_list)
            z_lum.extend(luminosity_list)

        # Bring arrays in numpy shape
        xy_age_mass = np.array(xy_age_mass)
        z_lum       = np.array(z_lum)

        # Define the interpolation grids
        grid_x, grid_y = np.mgrid[0.7:10000:100j, 0.1:1.4:100j]

        # Interpolate the luminosity from the 2D grid
        interpolated_luminosity = interpolate.griddata(xy_age_mass, z_lum, (time["star"]/1e+6, star_mass), method='linear', rescale=True)

    # Stellar constant, W m-2
    S_0    = interpolated_luminosity * L_sun / ( 4. * np.pi * (mean_distance*AU)**2. )

    # Scale instellation by fixed fraction
    S_0    = S_0 * Sfrac

    # Mean flux averaged over surface area, W m-2
    toa_heating             = ( 1. - albedo ) * S_0 / 4.
   
    return toa_heating

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    ##### Settings

    # Planet age and orbit
    time = { "planet": 0., "star": 4567e+6 } # yr,
    # time_current  = 0                 # yr, time after start of MO
    # time_offset   = 4567e+6           # yr, time relative to star formation
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 0.723                 # au, orbital distance

    # Surface pressure & temperature
    
    T_surf        = 800.                # K

    # # Volatile molar concentrations: must sum to ~1 !
    # P_surf        = 210e+5              # Pa
    # vol_list = { 
    #               "H2O"  : 100e5/P_surf,
    #               "CO2"  : 100e5/P_surf,
    #               "H2"   : 0., 
    #               "NH3"  : 100e5/P_surf,
    #               "N2"   : 10e5/P_surf,  
    #               "CH4"  : 0., 
    #               "O2"   : 0., 
    #               "CO"   : 0., 
    #               # # No thermodynamic data, RT only
    #               # "O3"   : 0.01, 
    #               # "N2O"  : 0.01, 
    #               # "NO"   : 0.01, 
    #               # "SO2"  : 0.01, 
    #               # "NO2"  : 0.01, 
    #               # "HNO3" : 0.01, 
    #               # "He"   : 0.01, 
    #               # "OCS"  : 0.01,
    #             }
    
    # Partial pressure guesses
    P_surf      = "calc"   
     # Volatiles considered
    vol_list    = { 
                          "H2O" :  1e+6,
                          "NH3" :  0.,
                          "CO2" :  0.,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  1e+5,
                          "H2"  :  0.,
                        }

    # Stellar heating on/off
    stellar_heating = True

    # Rayleigh scattering on/off
    rscatter = True

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, vol_list)

    # Compute stellar heating
    atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.toa_heating = 0.
    else:
        print("TOA heating:", round(atm.toa_heating), "W/m^2")

    # Compute heat flux
    atm_dry, atm_moist, atm_moist_timestep = RadConvEqm({"output": os.getcwd()+"/output", "rad_conv": os.getcwd()}, time, atm, [], [], standalone=True, cp_dry=False, moist_timestep=True, trpp=True, rscatter=rscatter) 

    print(len(atm_moist_timestep.p))

    # Plot abundances w/ TP structure
    ga.plot_adiabats(atm_moist)
