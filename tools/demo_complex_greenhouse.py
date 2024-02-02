#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from AEOLUS.modules.solve_pt import RadConvEqm
from utils.socrates import CleanOutputDir

from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum



def run_once(T_surf, dirs):

    # Planet 
    time = { "planet": 0., "star": 456e+6 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 1.0                 # au, orbital distance
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    P_surf       = 127.0 * 1e5
    vol_mixing = {
                    "H2O" : 0.91805,
                    "CO2" : 5.98710,
                    "H2"  : 2.37994,
                    "CO"  : 115.89,
                    "N2"  : 1.77739
                }
    tot = np.sum(list(vol_mixing.values()))
    for key in vol_mixing.keys():
        vol_mixing[key] /= tot
    
    # Rayleigh scattering on/off
    rscatter = False

    # Compute contribution function
    calc_cf = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 12.0     # Fixed tropopause value if not calculated dynamically
    
    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    do_cloud = False
    re   = 1.0e-5 # Effective radius of the droplets [m] (drizzle forms above 20 microns)
    lwm  = 0.8    # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
    clfr = 0.0 

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, re, lwm, clfr, vol_mixing=vol_mixing, trppT=trppT, do_cloud=do_cloud)

    # Compute stellar heating
    _, atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Do rad trans
    _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, do_cloud=do_cloud) 

    return [T_surf, atm_moist.net_flux[0], atm_moist.net_flux[-1]]


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
        dirs["aeolus"]+"/spectral_files/Mallard/Mallard",
        dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]+"runtime_spectral_file"
    )
    print(" ")

    skin_k = 2.0
    skin_d = 0.02
    tmp_magma = 2800.0

    samples = 24
    
    # Run AEOLUS in a loop to generate runaway curve
    print("Running AEOLUS...")
    Ts_arr = []
    toa_arr = []
    boa_arr = []
    skn_arr = []
    for Ts in np.linspace(600, 3000, samples):
        print("T_surf = %d K" % Ts)
        out = run_once(Ts, dirs)
        Ts_arr.append(out[0])
        toa_arr.append(out[1])
        boa_arr.append(out[2])
        skn_arr.append(skin_k / skin_d * (tmp_magma - Ts))
        print(" ")
    
    # Setup plot
    print("Making plot")
    fig,ax = plt.subplots(1,1)

    # Plot data
    lw = 2
    ax.axvline(x=tmp_magma,  color='firebrick', lw=lw,  label="Magma")
    ax.plot(Ts_arr, toa_arr, color='gold',      lw=lw,  label='TOA')
    ax.plot(Ts_arr, boa_arr, color='orchid',    lw=lw,  label='BOA')
    ax.plot(Ts_arr, skn_arr, color='teal',      lw=lw,  label='Skin')

    # Setup figure and save
    fig.legend(loc='upper left')
    ax.set_xlabel("Surface temperature [K]")
    ax.set_ylabel("Upward-directed flux [W m-2]")
    ax.set_yscale("symlog")
    fig.savefig(dirs["output"]+"complex_runaway_demo.pdf", bbox_inches='tight')

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    print("Done!")

