#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.solve_pt import *
from utils.socrates import CleanOutputDir

from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum

import utils.phys as phys

def run_once(sep, dirs, T_magma, P_surf, skin_d):

    # Planet 
    time = { "planet": 0., "star": 50e+6 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    P_top         = 0.1                  # Pa

    # Define volatiles by mole fractions
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }
    
    # Rayleigh scattering on/off
    rscatter = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?

    A_B = 0.175  # bond albedo
    inst_sf = 3.0/8.0
    
    ##### Function calls

    toa_heating = InterpolateStellarLuminosity(star_mass, time, sep)

    zenith_angle = 48.19 # cronin+14 (also for scaling by a factor of 3/8 ^^)

    # T_eqm = (inst * (1.0 - A_B) / (4.0 * phys.sigma))**(1.0/4.0)
    # T_trpp = T_eqm * (0.5**0.25)  # radiative skin temperature
    T_trpp = 0.01
    
    # print("T_trpp = %g K" % T_trpp)

    # Create atmosphere object
    atm = atmos(T_magma,  P_surf * 1e5, P_top, pl_radius, pl_mass, vol_mixing=vol_mixing, trppT=T_trpp)
    atm.albedo_pl = A_B
    atm.inst_sf = inst_sf
    atm.zenith_angle = zenith_angle
    atm.toa_heating = toa_heating
    atmos.skin_d = skin_d
    atm.tmp_magma = T_magma

    # Do rad trans
    atm = MCPA_CL(dirs, atm, trppD, rscatter, T_surf_max=9.0e99, T_surf_guess = T_trpp+100)
    # atm = MCPA(dirs, atm, False, trppD, rscatter)

    return [atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts]


if __name__=='__main__':

    print("Start")

    # Set up dirs
    dirs = {
            "aeolus": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+"/output/"
            }
    
    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        dirs["aeolus"]+"/spectral_files/Oak/Oak",
        dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]+"runtime_spectral_file"
    )

    # PARAMETERS
    P_surf  = 500.0   # surface pressure [bar]
    T_magma = 1500.0  # magma temperature [K]
    skin_d  = 1.0e-2  # conductive skin thickness [m]
    r_inner = 0.4     # inner orbital distane [AU]
    r_outer = 1.4     # outer orbital distance [AU]
    samples = 15       # number of samples
    logx    = False   # log x-axis?
    # /PARAMETERS
    
    # Run AEOLUS in a loop to generate data
    if logx:
        r_arr = np.logspace( np.log10(r_inner), np.log10(r_outer), samples)
    else:
        r_arr = np.linspace(r_inner, r_outer, samples)
    S0_arr  = []
    OLR_arr = []
    net_arr = []
    ts_arr  = []
    for r in r_arr:
        print("Orbital separation = %.2f AU" % r)
        out = run_once(r, dirs, T_magma, P_surf, skin_d)
        S0_arr.append(out[0] * -1.0)  # directed downwards => minus sign
        OLR_arr.append(out[1])
        net_arr.append(out[2])
        ts_arr.append(out[3])
        print(" ")
    
    # Setup plot
    lw = 2.5
    print("Making plots")
    
    fig,(ax1, ax2) = plt.subplots(2,1, figsize=(6,6))

    ax1.axvline(x=1.0,   color='deepskyblue', lw=lw/2, linestyle='dashed')  # earth
    ax1.axvline(x=0.723, color='goldenrod',   lw=lw/2, linestyle='dashed')  # venus

    ax1.text(r_arr[0], +1.0, "Cooling", verticalalignment='center', color='seagreen')
    ax1.text(r_arr[0], -1.0, "Warming", verticalalignment='center', color='seagreen')
    ax1.axhline(y=0,         color='black',     lw=0.7)  # zero flux

    ax1.plot(r_arr, S0_arr,  color='royalblue', lw=lw, label='ASF')    # absorbed stellar flux (sw down)
    ax1.plot(r_arr, OLR_arr, color='indianred', lw=lw, label='OLR')    # outgoing longwave radiation (lw up)
    ax1.plot(r_arr, net_arr, color='seagreen' , lw=lw, label='Net')    # net upward-directed flux
    
    ax1.legend(loc='center right')
    ax1.set_yscale("symlog")
    ax1.set_ylabel("Flux [W m$^{-2}$]")
    if logx:
        ax1.set_xscale("log")
    ax1.set_xticklabels([])

    yticks = ax1.get_yticks()
    yticks = np.delete(yticks, np.argwhere( yticks == 0.0 ))
    ax1.set_yticks(yticks)

    ax2.axvline(x=1.0,   color='deepskyblue', lw=lw/2, linestyle='dashed', label="Earth")  # earth
    ax2.axvline(x=0.723, color='goldenrod',   lw=lw/2, linestyle='dashed', label="Venus")  # venus

    ax2.plot(r_arr, np.ones(samples)*T_magma, color='orangered', lw=lw, label="$T_{surf}^{int}$") # magma temperature
    ax2.plot(r_arr, ts_arr,                   color='black',     lw=lw, label="$T_{surf}^{atm}$") # surface solution temperature 

    ax2.legend(loc='upper right')
    ax2.set_ylabel("Temperature [K]")
    if logx:
        ax2.set_xscale("log")

    ax2.set_xlabel("Orbital separation [AU]")
    fig.savefig(dirs["output"]+"instellation_demo.pdf", bbox_inches='tight')

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    print("Done!")

