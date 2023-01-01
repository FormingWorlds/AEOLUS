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
import time as t

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

def RadConvEqm(dirs, time, atm, loop_counter, COUPLER_options, standalone, cp_dry, trpp, calc_cf, rscatter):

    ### Moist/general adiabat
    atm_moist = compute_moist_adiabat(atm, dirs, standalone, trpp, calc_cf, rscatter)

    ### Dry adiabat
    if cp_dry == True:

        # Compute dry adiabat  w/ timestepping
        atm_dry   = compute_dry_adiabat(atm, dirs, standalone, calc_cf, rscatter)

        if standalone == True:
            print("Net, OLR => moist:", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print("| dry:", str(round(atm_dry.net_flux[0], 3)), str(round(atm_dry.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print()
    else: atm_dry = {}
    
    # Plot
    if standalone == True:
        #plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs)
        # Save to disk
        with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
            pkl.dump(atm_moist, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    return atm_dry, atm_moist

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
  
# Dry convective adjustment
def DryAdj(atm):

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

# # Moist convective adjustment
# def MoistAdj(atm, dT):

#     # Apply heating
#     tmp_heated = atm.tmp + dT

#     # Reset to moist adiabat if convectively unstable
#     for idx, tmp_heated in enumerate(tmp_heated):
#         atm.tmp[idx] = np.amax([tmp_heated, atm.tmp[idx]])

#     return atm 

#Moist adjustment routine.
def moist_adj(Tmid,pmid,pint,conv_timescale):
    dT_conv = np.zeros(len(pmid))
    Tmid_cc = np.zeros(len(pmid))
    for i in range(len(pmid)):
        Tmid_cc[i] = Tmid[i]
    nb_convsteps = 10
    for i_conv in range(nb_convsteps):
        did_adj = False
        #------------------- L = cst -------------------
        for i in range(len(Tmid_cc)-1): #Downward pass
            if (Tmid_cc[i] < ga.Tdew('H2O',pmid[i])):
                Tmid_cc[i]=ga.Tdew('H2O',pmid[i])
                did_adj = True
        for i in range(len(Tmid_cc)-2,-1,-1): #Upward pass
            if (Tmid_cc[i] < ga.Tdew('H2O',pmid[i])):
                Tmid_cc[i]=ga.Tdew('H2O',pmid[i])
                did_adj = True

        # If no adjustment required, exit the loop
        if (did_adj == False):
            break

    # Change in temperature is Tmid_cc - Tmid
    dT_conv[:] = (Tmid_cc[:] - Tmid[:])/conv_timescale
    return dT_conv

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

# Time integration for n steps
def compute_dry_adiabat(atm, dirs, standalone, calc_cf=False, rscatter=False):

    # Dry adiabat settings 
    rad_steps   = 100  # Maximum number of radiation steps
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
            atm_dry         = SocRadModel.radCompSoc(atm_dry, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter)
            
            dT_dry          = atm_dry.net_heating * atm_dry.dt

            # Limit the temperature change per step
            dT_dry          = np.where(dT_dry > dT_max, dT_max, dT_dry)
            dT_dry          = np.where(dT_dry < -dT_max, -dT_max, dT_dry)

            # Apply heating
            atm_dry.tmp     += dT_dry

            pure_steam_adj = True
            if pure_steam_adj:
                dT_moist = moist_adj(atm_dry.tmp,atm_dry.p,atm_dry.pl,1000.)
                atm_dry.tmp     += dT_moist

            # # Do the surface balance
            # kturb       = .1
            # atm.tmp[-1] += -atm.dt * kturb * (atm.tmp[-1] - atm.ts)
            
            # Dry convective adjustment
            for iadj in range(conv_steps):
                atm_dry     = DryAdj(atm_dry)

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


def compute_moist_adiabat(atm, dirs, standalone, trpp, calc_cf=False, rscatter=False):

    atm_moist = copy.deepcopy(atm)

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(atm_moist)

    # Run SOCRATES
    atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    if trpp == True:
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = SocRadModel.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=calc_cf, rscatter=rscatter)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist

def find_tropopause(atm_moist):

    # Heating criterion
    if atm_moist.trppT == 0:

        print("TROPOPAUSE SET BY HEATING")

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
        
            atm_moist.trppidx   = trpp_idx                  # index
            atm_moist.trppP     = atm_moist.pl[trpp_idx]    # pressure 
            atm_moist.trppT     = atm_moist.tmpl[trpp_idx]  # temperature

    # Temperature criterion
    else:

        print("TROPOPAUSE SET TO", atm_moist.trppT, "K")
        trpp_idx = (np.abs(atm_moist.tmp - atm_moist.trppT)).argmin()

        atm_moist.trppidx   = trpp_idx                  # index
        atm_moist.trppP     = atm_moist.pl[trpp_idx]    # pressure 
        atm_moist.trppT     = atm_moist.tmpl[trpp_idx]  # temperature


    

    return atm_moist


def set_stratosphere(atm):

    trpp_idx = int(atm.trppidx)
    trpp_prs = atm.trppP
    trpp_tmp = atm.trppT

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

    # If the age is not known with enough accuracy, define the ISR manually
    L_star                  = 0.000553*L_sun                                           # Trappist-1 luminosity
    toa_heating             = ( 1. - 0. ) * (L_star/(4.*np.pi*(mean_distance*AU)**2.)) # Tidally-locked, zero albedo
   
    return toa_heating

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    start = t.time()
    ##### Settings

    # Planet age and orbit
    time = { "planet": 0., "star": 4567e+6 } # yr,
    # time_current  = 0                 # yr, time after start of MO
    # time_offset   = 4567e+6           # yr, time relative to star formation
    star_mass     = 0.08    #1.0                 # M_sun, mass of star
    mean_distance = 0.01154 #1.0                 # au, orbital distance

    # Surface pressure & temperature
    
    T_surf        = 1225.0                # K

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
    #P_surf      = "calc" # does not do anything  
     # Volatiles considered
    vol_list    = { 
                          "H2O" :  1.00e+6,
                          "NH3" :  0.,
                          "CO2" :  0.,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  0.,
                          "H2"  :  0.
                        }

    P_surf      = sum(vol_list.values())

    # Stellar heating on/off
    stellar_heating = True

    # Rayleigh scattering on/off
    rscatter = True

    # Compute contribution function
    calc_cf = False

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, vol_list, calc_cf=calc_cf)

    # Compute stellar heating
    atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.toa_heating = 0.
    else:
        print("TOA heating:", round(atm.toa_heating), "W/m^2")

    # Compute heat flux
    atm_dry, atm_moist = RadConvEqm({"output": os.getcwd()+"/output", "rad_conv": os.getcwd()}, time, atm, [], [], standalone=True, cp_dry=True, trpp=False, calc_cf=calc_cf, rscatter=rscatter) 

    #print(len(atm_moist.p))

    # Plot abundances w/ TP structure
    #ga.plot_adiabats(atm_moist)

    end = t.time()
    print(end - start)

    ClausiusClapeyron = [ga.Tdew('H2O',p) for p in atm_dry.p]

    atm_print = atm_dry # atm_moist

    print("OPR = ", atm_print.flux_up_total[0], "ASR = ", atm_print.flux_down_total[0], "ASR-OPR = ", atm_print.flux_down_total[0]-atm_print.flux_up_total[0])

    np.savetxt(f'data_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.pl,atm_print.tmpl,atm_print.SW_flux_up,atm_print.LW_flux_up,atm_print.SW_flux_down,atm_print.LW_flux_down,atm_print.SW_flux_net,atm_print.LW_flux_net,atm_print.flux_up_total,atm_print.flux_down_total,atm_print.net_flux)))

    np.savetxt(f'data_heating_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.SW_heating,atm_print.LW_heating,atm_print.net_heating)))

    np.savetxt(f'data_spectral_TOA_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.LW_spectral_flux_up[:,0],atm_print.SW_spectral_flux_up[:,0],atm_print.net_spectral_flux[:,0])))

    np.savetxt(f'data_spectral_Surface_{int(P_surf)}_{int(T_surf)}.dat',np.column_stack((atm_print.LW_spectral_flux_up[:,-1],atm_print.SW_spectral_flux_up[:,-1],atm_print.net_spectral_flux[:,-1])))
