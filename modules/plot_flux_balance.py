#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:14:29 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import matplotlib.pyplot as plt
import numpy as np
import math 
#import json

from modules.spectral_planck_surface import surf_Planck_nu

try:
    import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
except:
    import atm_rad_conv.GeneralAdiabat as ga

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
    trpp_idx = int(atm_moist.trppidx)
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
