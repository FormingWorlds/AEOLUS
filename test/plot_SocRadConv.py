#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 09:47:56 2023

@author: x_anfey
"""

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
    import SocRadConv
    from atmosphere_column import atmos
except:
    import atm_rad_conv.phys as phys
    import atm_rad_conv.GeneralAdiabat as ga
    import atm_rad_conv.SocRadModel as SocRadModel
    import atm_rad_conv.SocRadConv as SocRadConv
    from atm_rad_conv.atmosphere_column import atmos


################################ SETTINGS ###################################################

# Planet age and orbit
time = { "planet": 0., "star": 4567e+6 } # yr,
# time_current  = 0                 # yr, time after start of MO
# time_offset   = 4567e+6           # yr, time relative to star formation
star_mass     = 1.0                 # M_sun, mass of star
mean_distance = 1.0                 # au, orbital distance


# Volatiles considered
## partial pressures
pp_water = 1e+3  # Water partial pressure for 1% contribution in a Earth-like atmosphere

vol_list_pp    = { 
                         "H2O" :  pp_water,
                         "NH3" :  0.,
                         "CO2" :  3.5e+1,
                         "CH4" :  0.,
                         "CO"  :  0.,
                         "O2"  :  2.07e+4,
                         "N2"  :  7.73e+4,
                         "H2"  :  0.,
                       }

## molar concentrations
vol_list_mc    = { 
                         "H2O" :  1.0e-2,
                         "NH3" :  0.,
                         "CO2" :  3.5e-4,
                         "CH4" :  0.,
                         "CO"  :  0.,
                         "O2"  :  2.07e-1,
                         "N2"  :  7.73e-1,
                         "H2"  :  0.,
                       }

vol_list_mc_bis    = { 
                         "H2O" :  1.,
                         "NH3" :  0.,
                         "CO2" :  0.,
                         "CH4" :  0.,
                         "CO"  :  0.,
                         "O2"  :  0.,
                         "N2"  :  0.,
                         "H2"  :  0.,
                       }

# Stellar heating on/off
stellar_heating = True

# Rayleigh scattering on/off
rscatter = True

# Compute contribution function
calc_cf = False

# Instellation scaling | 1.0 == no scaling
Sfrac = 1.0

# surface temperatures and pressures for wich to run the code
#T_grid = [288,500,1000,2000,3000]        # K
#P_grid = [1e+5,1e+6,5e+6,1e+7,2e+7]      # Pa

T_grid=[500,1000,1500,2000]
P_grid=[1e+6]

# Increase in Total Pressure : in proportion to molar concentrations (p) or by increasing only water contribution (w)?
ITP = 'p'


#####################################  FUNCTION DEFINITION  ############################################

# Adapted version of SocRadConv.RadConvEqm
def RadConvEqm(dirs, T_surf,P_surf, atm, loop_counter, COUPLER_options, standalone, cp_dry, trpp, calc_cf, rscatter):

    ### Moist/general adiabat
    atm_moist = SocRadConv.compute_moist_adiabat(atm, dirs, standalone, trpp, calc_cf, rscatter)

    ### Dry adiabat
    if cp_dry == True:

        # Compute dry adiabat  w/ timestepping
        atm_dry   = SocRadConv.compute_dry_adiabat(atm, dirs, standalone, calc_cf, rscatter)

    else: atm_dry = {}
    
    plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs)
    # Save to disk
    with open(dirs["subfolder"]+"/T"+str(T_surf)+"_P"+str(P_surf/1e+5)+"_atm.pkl", "wb") as atm_file: 
        pkl.dump(atm_moist, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    return atm_dry, atm_moist

# Adapted version of SocRadConv.plot_flux_balance
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
    ax1.set_xlabel(r'Temperature $T$ [K]')
    ax1.set_ylabel(r'Pressure $P$ [Pa]')
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
    ax2.set_xlabel(r'Outgoing flux $F^{\uparrow}$ [W m$^{-2}$]')
    ax2.set_ylabel(r'Pressure $P$ [Pa]')
    ax2.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)

    # Wavenumber vs. OLR
    #ax3.plot(atm_moist.band_centres, SocRadConv.surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label=str(round(atm_moist.ts))+' K blackbody')
    if cp_dry == True: ax3.plot(atm_dry.band_centres, atm_dry.net_spectral_flux[:,0]/atm_dry.band_widths, color=ga.vol_colors[col_vol3][col_idx])
    #ax3.plot(atm_moist.band_centres, atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
    
    ax3.plot(atm_moist.band_centres, (atm_moist.LW_spectral_flux_up[:,0]+atm_moist.SW_spectral_flux_up[:,0])/atm_moist.band_widths, color="red",label='TOA')
    #ax3.plot(atm_moist.band_centres, (atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths, color="red",ls=':',label='BOA')
    
    ax3.set_ylabel(r'Spectral flux density [W m$^{-2}$ [cm$^{-1}]^{-1}$]')
    ax3.set_xlabel(r'Wavenumber [cm$^{-1}$]')
    ax3.legend(loc=1)
    ymax_plot = 1.2*np.max((atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths)
    #ax3.set_ylim(bottom=0, top=ymax_plot)
    ax3.set_ylim(bottom=0, top=7.5)
    #ax3.set_xlim(left=0, right=np.max(np.where((atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
    ax3.set_xlim(left=0, right=10000)
    

    # Heating vs. pressure
    ax4.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

    if cp_dry == True: 
        ax4.plot(atm_dry.LW_heating, atm_dry.p, ls="--", color=ga.vol_colors[col_vol3][col_idx+1])
        ax4.plot(atm_dry.net_heating, atm_dry.p, lw=2, color=ga.vol_colors[col_vol3][col_idx+1])
        ax4.plot(atm_dry.SW_heating, atm_dry.p, ls=":", color=ga.vol_colors[col_vol3][col_idx+1])

    ax4.plot(atm_moist.LW_heating, atm_moist.p, ls="--", color=ga.vol_colors[col_vol1][col_idx+1], label=r'LW')
    ax4.plot(atm_moist.net_heating, atm_moist.p, lw=2, color=ga.vol_colors[col_vol1][col_idx+1], label=r'Net')
    ax4.plot(atm_moist.SW_heating, atm_moist.p, ls=":", color=ga.vol_colors[col_vol1][col_idx+1], label=r'SW')

    # Plot tropopause
    trpp_idx = int(atm_moist.trppidx)
    if trpp_idx > 0:
        ax4.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol2][col_idx], lw=1.0, ls="-.", label=r'Tropopause')
    
    ax4.invert_yaxis()
    ax4.legend(ncol=4, fontsize=10, loc=3)
    ax4.set_ylabel(r'Pressure $P$ [Pa]')
    ax4.set_xlabel(r'Heating [K/day]')
    # ax4.set_xscale("log")
    ax4.set_yscale("log") 
    x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
    x_minmax = np.max([ 20, x_minmax ])
    if not math.isnan(x_minmax):
        ax4.set_xlim(left=-x_minmax, right=x_minmax)
    ax4.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)
    ax4.set_xlim(left=-200,right=200)
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

    # plt.show()
    plt.savefig(dirs["subfolder"]+"/flux_T"+str(T_surf)+"_P"+str(P_surf/1e+5)+'.pdf', bbox_inches="tight")
    plt.close(fig)

    # with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
    #     pkl.dump(atm, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    # # Save atm object to .json file
    # json_atm = json.dumps(atm.__dict__)
    # with open(dirs["output"]+"/"+str(int(time_current))+"_atm.json", "wb") as atm_file:
    #     json.dump(json_atm, atm_file)

# Adapted version of GeneralAdiabat.plot_adiabats
def plot_adiabats(dirs,atm):

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
            Psat_array = [ ga.p_sat(vol, T) for T in T_sat_array ]
            ax1.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+ga.vol_latex[vol], lw=ls_ind, ls=":", color=ga.vol_colors[vol][4])

            # Plot partial pressures
            ax1.semilogy(atm.tmp, atm.p_vol[vol], color=ga.vol_colors[vol][4], lw=ls_ind, ls="-", label=r'$p$'+ga.vol_latex[vol],alpha=0.99)

            # Sum up partial pressures
            p_partial_sum += atm.p_vol[vol]

            # Plot individual molar concentrations
            ax2.semilogy(atm.x_cond[vol],atm.p, color=ga.vol_colors[vol][4], lw=ls_ind, ls="--", label=ga.vol_latex[vol]+" cond.")
            ax2.semilogy(atm.x_gas[vol],atm.p, color=ga.vol_colors[vol][4], lw=ls_ind, ls="-", label=ga.vol_latex[vol]+" gas")
            
    # # Plot sum of partial pressures as check
    #ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="--", label=r'$\sum p^\mathrm{i}$',alpha=0.99)
    
    # # # Dry adiabat function from RTB book
    # # ! careful: non-T dependent CPs used, can lead to differences with general adiabat !
    ax1.semilogy( ga.dry_adiabat( atm.ts, atm.pl, atm.cp[-1]), atm.pl , color=ga.vol_colors["black_3"], ls="-.", lw=ls_dry, label=r'Dry adiabat') # Functional form

    # General moist adiabat
    ax1.semilogy(atm.tmpl, atm.pl, color=ga.vol_colors["black_1"], lw=ls_moist,label="General adiabat",alpha=0.99)

    # Phase molar concentrations
    ax2.semilogy(atm.xd+atm.xv,atm.p, color=ga.vol_colors["black_2"], lw=ls_ind, ls=":", label=r"Gas phase")

    fs_l = 16
    fs_m = 14
    fs_s = 12

    ax1.invert_yaxis()
    ax1.set_xlabel(r'Temperature $T$ [K]', fontsize=fs_l)
    ax1.set_ylabel(r'Pressure $P$ [Pa]', fontsize=fs_l)
    # ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes', fontsize=fs_l)
    ax1.legend(loc=1, ncol=np.min([len(atm.vol_list)+1,2]), fontsize=fs_s)
    ax1.set_xlim([0,np.max(atm.ts)])

    ax2.invert_yaxis()
    # ax2.set_title('Phase & species abundances', fontsize=fs_l)
    ax2.set_xlabel(r'Molar concentration $X^{\mathrm{i}}_{\mathrm{phase}}$', fontsize=fs_l)
    ax2.set_ylabel(r'Pressure $P$ [Pa]', fontsize=fs_l)
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
    fig.suptitle(r'T_$\rm{surf}$=%.0f K'%atm.ts)
    #plt.show()

    plt.savefig(dirs["subfolder"]+"/general_adiabat_T"+str(T_surf)+"_P"+str(P_surf/1e+5)+'.pdf', bbox_inches='tight')
    plt.close(fig)  

# New function
def plot_OLRs(dirs,OLRs):
    
    plt.plot(atm_moist.band_centres, SocRadConv.surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label=str(round(atm_moist.ts))+' K blackbody')
    #if cp_dry == True: plt.plot(atm_dry.band_centres, atm_dry.net_spectral_flux[:,0]/atm_dry.band_widths, color=ga.vol_colors[col_vol3][col_idx])
    #ax3.plot(atm_moist.band_centres, atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
    
    for OLR in OLRs:
        plt.plot(atm_moist.band_centres, (OLR[1][0]+OLR[1][1])/atm_moist.band_widths,label=str(OLR[0]/1e+5)+'bar')
    
    plt.ylabel(r'Spectral flux density [W m$^{-2}$ cm$^{-1}$]')
    plt.xlabel(r'Wavenumber [cm$^{-1}$]')
    plt.legend(loc=1)
    ymax_plot = 1.2*np.max(OLRs[0][1]/atm_moist.band_widths)
    plt.ylim(bottom=0, top=ymax_plot)
    plt.xlim(left=0, right=np.max(np.where(OLRs[0][1]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
    # plt.xlim(left=0, right=5000)
    
    plt.savefig(dirs["subfolder"]+'/OLRs_T'+str(T_surf)+'.pdf', bbox_inches="tight")
    plt.close()

##################################### FUNCTION CALLS & OTHERS ##########################################


folder = "T"+str(T_grid[0])+"_"+str(T_grid[-1])+"__P"+str(P_grid[0]/1e+5)+"_"+str(P_grid[-1]/1e+5)+"_"+ITP # pressure given in bar
"""
os.system("mkdir output/"+"HITRAN")
os.system("mkdir output/HITRAN/"+folder)

if ITP=='p':
    vol_list=vol_list_mc_bis
    for T_surf in T_grid:
        subfolder="T"+str(T_surf)+"_"+ITP
        os.system("mkdir output/HITRAN/"+folder+'/'+subfolder)
        OLRs=[]
        for P_surf in P_grid:
            atm = atmos(T_surf, P_surf, vol_list, calc_cf=calc_cf)
            atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)
            
            # Set stellar heating on or off
            if stellar_heating == False: 
                atm.toa_heating = 0.
            else:
                print("TOA heating:", round(atm.toa_heating), "W/m^2")

            # Compute heat flux
            atm_dry, atm_moist = RadConvEqm({"output": os.getcwd()+"/output/HITRAN","subfolder": os.getcwd()+"/output/HITRAN/"+folder+"/"+subfolder, "rad_conv": os.getcwd()}, T_surf, P_surf, atm, [], [], standalone=True, cp_dry=False, trpp=False, calc_cf=calc_cf, rscatter=rscatter) 
            
            spectral_flux_up=[atm_moist.LW_spectral_flux_up[:,0],atm_moist.SW_spectral_flux_up[:,0]]
            OLRs.append([P_surf,spectral_flux_up])
            
            # Plot abundances w/ TP structure
            plot_adiabats({"output": os.getcwd()+"/output/HITRAN","subfolder": os.getcwd()+"/output/HITRAN/"+folder+"/"+subfolder}, atm_moist)
        plot_OLRs({"output": os.getcwd()+"/output/HITRAN","subfolder": os.getcwd()+"/output/HITRAN/"+folder+"/"+subfolder}, OLRs)

"""

if ITP=='w':
    vol_list=vol_list_pp
    for T_surf in T_grid:
        subfolder="T"+str(T_surf)+"_"+ITP
        os.system("mkdir output/"+folder+'/'+subfolder)
        OLRs=[]
        for P_surf in P_grid:
            pp_water=1e+3+(P_surf-1e+5)
            P="calc"
            atm = atmos(T_surf, P, vol_list, calc_cf=calc_cf)
            atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)
            
            # Set stellar heating on or off
            if stellar_heating == False: 
                atm.toa_heating = 0.
            else:
                print("TOA heating:", round(atm.toa_heating), "W/m^2")

            # Compute heat flux
            atm_dry, atm_moist = RadConvEqm({"output": os.getcwd()+"/output","subfolder": os.getcwd()+"/output/"+folder+"/"+subfolder, "rad_conv": os.getcwd()}, T_surf, P_surf, atm, [], [], standalone=True, cp_dry=False, trpp=False, calc_cf=calc_cf, rscatter=rscatter) 
            
            spectral_flux_up=[atm_moist.LW_spectral_flux_up[:,0],atm_moist.SW_spectral_flux_up[:,0]]
            OLRs.append([P_surf,spectral_flux_up])
            
            # Plot abundances w/ TP structure
            plot_adiabats({"output": os.getcwd()+"/output","subfolder": os.getcwd()+"/output/"+folder+"/"+subfolder}, atm_moist)
        plot_OLRs({"output": os.getcwd()+"/output","subfolder": os.getcwd()+"/output/"+folder+"/"+subfolder}, OLRs)
           

def plot_OLR(T):       
    with open('/home/x_anfey/AEOLUS/output/HITRAN/T500_2000__P10.0_10.0_p/T'+str(T)+'_p/T'+str(T)+'_P10.0_atm.pkl','rb') as f: 
        content_HITRAN = pkl.load(f)  
    
    plt.plot(content_HITRAN.band_centres, (content_HITRAN.LW_spectral_flux_up[:,0]+content_HITRAN.SW_spectral_flux_up[:,0])/content_HITRAN.band_widths, color="black",label='HITRAN')
    #ax3.plot(atm_moist.band_centres, (atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths, color="red",ls=':',label='BOA')
    
    with open('/home/x_anfey/AEOLUS/output/ExoMol/T500_2000__P10.0_10.0_p/T'+str(T)+'_p/T'+str(T)+'_P10.0_atm.pkl','rb') as f: 
        content_ExoMol = pkl.load(f)  
    
    plt.plot(content_ExoMol.band_centres, (content_ExoMol.LW_spectral_flux_up[:,0]+content_ExoMol.SW_spectral_flux_up[:,0])/content_ExoMol.band_widths, color="red",label='ExoMol')
    
    
    plt.ylabel(r'Spectral flux density [W m$^{-2}$ [cm$^{-1}]^{-1}$]')
    plt.xlabel(r'Wavenumber [cm$^{-1}$]')
    plt.legend(loc=1)
    plt.title(r'$(P,T)=(10~bar,$'+str(T)+'$~K)$')
    #ymax_plot = 1.2*np.max((atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths)
    #ax3.set_ylim(bottom=0, top=ymax_plot)
    #plt.ylim(bottom=0, top=7.5)
    #ax3.set_xlim(left=0, right=np.max(np.where((atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
    plt.xlim(left=0, right=7000)
            
def plot_heating(T):       
    with open('/home/x_anfey/AEOLUS/output/HITRAN/T500_2000__P10.0_10.0_p/T'+str(T)+'_p/T'+str(T)+'_P10.0_atm.pkl','rb') as f: 
        content_HITRAN = pkl.load(f)  
    
    plt.plot(content_HITRAN.net_heating,content_HITRAN.p, color="black",label='HITRAN')
    #ax3.plot(atm_moist.band_centres, (atm_moist.LW_spectral_flux_up[:,-1]+atm_moist.SW_spectral_flux_up[:,-1])/atm_moist.band_widths, color="red",ls=':',label='BOA')
    
    with open('/home/x_anfey/AEOLUS/output/ExoMol/T500_2000__P10.0_10.0_p/T'+str(T)+'_p/T'+str(T)+'_P10.0_atm.pkl','rb') as f: 
        content_ExoMol = pkl.load(f)  
    
    plt.plot(content_ExoMol.net_heating,content_ExoMol.p, color="red",label='ExoMol')
    
    
    plt.gca().invert_yaxis()
    plt.legend()
    plt.yscale('log')
    plt.axvline(x=0,color='grey',linestyle='--')
    plt.xlim(-200,200)
    plt.ylabel(r'Pressure $P$ [Pa]')
    plt.xlabel(r'Heating [K/day]')
    plt.title(r'$(P,T)=(10~bar,$'+str(T)+'$~K)$')
    plt.show()          
            
            
            