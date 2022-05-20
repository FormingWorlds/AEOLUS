import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
import copy
import pathlib
import pickle as pkl
import json
import glob, re, os
import seaborn as sns
import phys
import GeneralAdiabat as ga
import SocRadModel
import SocRadConv
from atmosphere_column import atmos


#### PLOT SETTINGS

ls_moist    = 2.5
ls_dry      = 2.0
ls_ind      = 1.5

# fig, ((Aax1, Aax2, Aax3), (Bax1, Bax2, Bax3), (Cax1, Cax2, Cax3)) = plt.subplots(3, 3, figsize=(13,9))
# fig, ((Aax1, Bax1, Cax1), (Aax3, Bax3, Cax3)) = plt.subplots(2, 3, figsize=(13,6))
# fig, ax = plt.subplots(1, 1, figsize=(8,6))
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13,6))
plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=0.2, wspace=0.2)
sns.despine()

# Line settings
col_idx  = 5
col_vol = "H2O"

ls_list = [ "-", "--", ":", "-." ]
lw      = 1.5
fs_l = 15
fs_m = 13
fs_s = 11

dirs =  {
               "output":   os.getcwd()+"/", 
               "rad_conv": "/Users/timlichtenberg/git/proteus/atm_rad_conv"
            }

#### LOOP OVER PARAMETERS
for set_idx, setting in enumerate([ "set1", "set2", "set3" ]): # "set1", "set2", "set3"

    # Retained condensate fraction
    for alpha_idx, alpha_cloud in enumerate([ 0.0 ]): # 0.0, 0.1, 1.0

        ls = ls_list[alpha_idx]

        ### Initial conditions

        # Magma ocean case
        if setting == "set1":
            name = "Bistable"
            time = { "planet": 0., "star": 4567e+06 } # yr,
            Mstar           = 1.0 
            mean_distance   = 1.0
            P_surf          = "calc"
            T_surf          = 280
            Sfrac=0.4
            vol_dict    = { 
                          "H2O" :  1e+5,
                          "NH3" :  0.,
                          "CO2" :  2e+5,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  0e+5,
                          "H2"  :  0e+5,
                        }

            col_vol = "CO2"

        # LATE VENEER case
        if setting == "set2":
            name = "Oscillatory"
            time = { "planet": 0., "star": 4567e+06 } # yr,
            Mstar           = 1.0 
            mean_distance   = 1.0
            
            P_surf          = "calc"
            T_surf          = 340
            Sfrac=0.45
 
            vol_dict        = { 
                              "H2O" :  1e+5,
                              "NH3" :  0e+5,
                              "CO2" :  20e+5,
                              "CH4" :  0e+5,
                              "CO"  :  0e+5,
                              "O2"  :  0e+5,
                              "N2"  :  0e+5,
                              "H2"  :  0e+5,
                            }

            col_vol = "H2O"

        # Hadean Earth / OHZ case
        if setting == "set3":

            name = r"Condensing"

            # Planet age and orbit
            time = { "planet": 0., "star": 4567e+06 } # yr,

            # Star mass, M_sun
            Mstar           = 1.0 

            # Planet-star distance, au
            mean_distance   = 1.0

            # Surface pressure range (Pa)
            P_surf          = "calc"

            # Surface temperature range (K)
            T_surf          = 270

            Sfrac=0.43

            # Volatiles considered
            vol_dict    = { 
                              "H2O" :  1e+5,
                              "NH3" :  0e+5,
                              "CO2" :  32e+5,
                              "CH4" :  0e+5,
                              "CO"  :  0e+5,
                              "O2"  :  0e+5,
                              "N2"  :  0e+5,
                              "H2"  :  0e+5,
                        }


            # Plot color
            col_vol = "H2"


        print("---------", setting, "case:", name, "alpha:", alpha_cloud)

        # Create atmosphere object
        atm                = atmos(T_surf, P_surf, vol_dict)

        # print(P_surf, atm.vol_list)

        # Set fraction of condensate retained in column
        atm.alpha_cloud       = alpha_cloud

        # Stellar heating
        atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, mean_distance, atm.albedo_pl, Sfrac=Sfrac)

        print("TOA heating:", round(atm.toa_heating), "W/m^2")

        atm.trppT = 150

        # Compute heat flux
        atm_dry, atm_moist = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=True, cp_dry=False, trpp=True, calc_cf=False, rscatter=True) 

        # Wavelength
        spectral_flux_up_total = atm_moist.LW_spectral_flux_up + atm_moist.SW_spectral_flux_up
        OLR_cm_moist = spectral_flux_up_total[:,0]/atm_moist.band_widths
        wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ]          # microns
        OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ]                    # microns
        ax.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol][col_idx], lw=1.0, label=name, ls=ls)
        ax2.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol][col_idx], lw=1.0, label=name, ls=ls)



# Zero line
# Aax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Bax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Cax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Aax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
# Bax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

ax.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')
ax.set_xscale("log")
ax.set_yscale("log") 
ax.set_xlim(left=0.5, right=1)
# ax.set_ylim(bottom=1e-20, top=1e5)
# ax.set_yticks([1e-10, 1e-5, 1e0, 1e5])
ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xticklabels(["0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
ax.legend(fontsize=10, loc=2)
ax.set_ylim(bottom=1e-22, top=1e-11)

ax2.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax2.set_xscale("log")
ax2.set_yscale("log") 
ax2.set_xlim(left=5, right=25)
ax2.set_xticks([6, 8, 10, 13, 17, 20, 25])
ax2.set_xticklabels(["6", "8", "10", "13", "17", "20", "25"])
ax2.set_ylim(bottom=9e+0, top=4e3)

# Subplot labels
lx = 0.98
ly = 0.015
ax.text(lx, ly, r'Shortwave', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=ax.transAxes)
ax2.text(lx, ly, r'IR', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=ax2.transAxes)
# Aax3.text(lx, ly, r'A2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax3.transAxes)
# Bax1.text(lx, ly, r'B1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax1.transAxes)
# # Bax2.text(lx, ly, r'B2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax2.transAxes)
# Bax3.text(lx, ly, r'B2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax3.transAxes)
# Cax1.text(lx, ly, r'C1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax1.transAxes)
# # Cax2.text(lx, ly, r'C2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax2.transAxes)
# Cax3.text(lx, ly, r'C2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax3.transAxes)

# Axes labels
# Aax1.set_ylabel(r'Pressure $P$ (bar)')


sns.despine()

plt.savefig(dirs["output"]+"RT.pdf", bbox_inches="tight")
plt.close(fig)
