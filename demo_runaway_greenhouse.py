#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.solve_pt import RadConvEqm
from utils.socrates import CleanOutputDir

from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum



def run_once(T_surf, dirs):

    # Planet 
    time = { "planet": 0., "star": 4.5e9 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 0.5                 # au, orbital distance
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    P_surf       = 300 * 1e5
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }
    
    # Rayleigh scattering on/off
    rscatter = False

    # Tropopause calculation
    trppT = 0.0     # Fixed tropopause value if not calculated dynamically
    
    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, vol_mixing=vol_mixing, trppT=trppT)

    # Compute stellar heating
    atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)

    # Do rad trans
    _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=False, calc_cf=False, rscatter=rscatter) 

    return [T_surf, atm_moist.LW_flux_up[0]]


if __name__=='__main__':

    print("Start")
    print(" ")

    # Set up dirs
    if os.environ.get('AEOLUS_DIR') == None:
        raise Exception("Environment variables not set! Have you sourced AEOLUS.env?")
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
    print(" ")

    
    # Run AEOLUS in a loop to generate runaway curve
    print("Running AEOLUS...")
    Ts_arr = []
    OLR_arr = []
    for Ts in np.linspace(200, 2200, 25):
        print("T_surf = %d K" % Ts)
        out = run_once(Ts, dirs)
        Ts_arr.append(out[0])
        OLR_arr.append(out[1])
        print(" ")
    
    # Get literature data
    g2013 = np.loadtxt(dirs["aeolus"]+"plotting_tools/comparison_data/Goldblatt13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    k2013 = np.loadtxt(dirs["aeolus"]+"plotting_tools/comparison_data/Kopparapu13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    h2015 = np.loadtxt(dirs["aeolus"]+"plotting_tools/comparison_data/Hamano15_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 

    # Setup plot
    print("Making plot")
    fig,ax = plt.subplots(1,1)

    # Plot data
    lw = 2
    ax.plot(k2013[0], k2013[1], color='tab:red',   lw=lw, label='Kopparapu+2013')
    ax.plot(g2013[0], g2013[1], color='tab:green', lw=lw, label='Goldblatt+2013')
    ax.plot(h2015[0], h2015[1], color='tab:blue',  lw=lw, label='Hamano+2015')
    ax.plot(Ts_arr, OLR_arr,    color='black',     lw=lw, label='AEOLUS')

    # Setup figure and save
    fig.legend(loc='upper center')
    ax.set_xlabel("Surface temperature [K]")
    ax.set_ylabel("OLR [W m-2]")
    fig.savefig(dirs["output"]+"runaway_demo.pdf")
    print(" ")

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    print("Done!")

