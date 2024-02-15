#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams['axes.formatter.useoffset'] = False
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

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
    time = { "planet": 0., "star": 150e+6 } # yr,
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
    rscatter = True

    # Tropopause calculation
    trppD = False   # Calculate dynamically?

    A_B = 0.2  # bond albedo
    inst_sf = 3.0/8.0
    
    ##### Function calls

    S_0 = InterpolateStellarLuminosity(star_mass, time, sep) 

    zenith_angle = 48.19 # cronin+14 (also for scaling by a factor of 3/8 ^^)

    T_eqm = (S_0 * inst_sf * (1.0 - A_B) /phys.sigma)**(1.0/4.0)
    T_trpp = T_eqm * (0.5**0.25)  # radiative skin temperature
    # T_trpp = 0.01
    
    # print("T_trpp = %g K" % T_trpp)

    # Create atmosphere object
    atm = atmos(T_magma,  P_surf * 1e5, P_top, pl_radius, pl_mass, vol_mixing=vol_mixing, trppT=T_trpp, req_levels=150)
    atm.albedo_pl = A_B
    atm.inst_sf = inst_sf
    atm.zenith_angle = zenith_angle
    atm.instellation = S_0
    atm.skin_d = skin_d
    atm.tmp_magma = T_magma

    # Do rad trans
    atm = MCPA_CBL(dirs, atm, trppD, rscatter, T_surf_max=9.0e99, T_surf_guess = T_trpp+100)
    # atm = MCPA(dirs, atm, False, trppD, rscatter)

    # Plot case 
    plt.ioff()
    fig,ax = plt.subplots(1,1)
    ax.plot(atm.tmpl, atm.pl, color='black', lw=2)
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_ylabel("Pressure [Pa]")
    ax.set_xlabel("Temperature")
    ax.set_title("a = %g AU" % sep)
    fig.savefig(dirs["output"]+"/recent.jpg",bbox_inches='tight', dpi=100)
    plt.close()

    # Save netcdf
    atm.write_ncdf(dirs["output"]+"/recent.nc")

    return [atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts, T_trpp]


if __name__=='__main__':

    print("Start")

    # Planet data
    x_planets = {
        "Earth": 1.0,
        "Venus": 0.723,
        "Mercury": 0.387
    }
    c_planets = {
        "Earth": 'deepskyblue',
        "Venus": 'goldenrod',
        "Mercury": 'violet'
    }

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

    # PARAMETERS
    P_surf  = 280.0   # surface pressure [bar]
    T_magma = 3000.0  # magma temperature [K]
    skin_d  = 1e-2  # conductive skin thickness [m]
    r_inner = 0.3     # inner orbital distane [AU]
    r_outer = 1.4     # outer orbital distance [AU]
    samples = 10       # number of samples
    logx    = False   # log x-axis?
    legend  = True    # make legend?
    dx_tick = 0.1     # x-tick spacing (set to 0 for automatic)
    # /PARAMETERS
    
    # Run AEOLUS in a loop to generate data
    if logx:
        r_arr = np.logspace( np.log10(r_inner), np.log10(r_outer), samples)
    else:
        r_arr = np.linspace(r_inner, r_outer, samples)
    asf_arr  = []
    OLR_arr = []
    net_arr = []
    ts_arr  = []
    tr_arr  = []
    for r in r_arr:
        print("Orbital separation = %.2f AU" % r)
        out = run_once(r, dirs, T_magma, P_surf, skin_d)
        asf_arr.append(out[0] * -1.0)  # directed downwards => minus sign
        OLR_arr.append(out[1])
        net_arr.append(out[2])
        ts_arr.append(out[3])
        tr_arr.append(out[4])
        print(" ")

    save_arr = [r_arr, asf_arr, OLR_arr, net_arr, ts_arr, tr_arr]
    np.savetxt(dirs["output"]+"/data_%dK.csv"%T_magma, 
               np.array(save_arr).T, fmt="%.5e", delimiter=",",
               header="r [AU],  S_0 [W m-2], OLR [W m-2], net [W m-2], ts [K], tr[K] ")
    
    # Setup plot
    lw = 2.5
    print("Making plots")

    plt.ioff()
    fig,(ax1, ax2) = plt.subplots(2,1, figsize=(6,6.4))

    ax1.set_title("%d K" % T_magma)

    ax1.text(r_arr[0], +1.0, "Cooling", verticalalignment='center', color='seagreen', weight='bold', zorder=8).set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
    ax1.text(r_arr[0], -1.0, "Warming", verticalalignment='center', color='seagreen', weight='bold', zorder=8).set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
    ax1.axhline(y=0, color='black', lw=0.7, zorder=0)  # zero flux

    ax1.plot(r_arr, asf_arr, color='royalblue', lw=lw, label='ASF')    # absorbed stellar flux (sw down)
    ax1.plot(r_arr, OLR_arr, color='crimson',   lw=lw, label='OLR')    # outgoing longwave radiation (lw up)
    ax1.plot(r_arr, net_arr, color='seagreen' , lw=lw, label='Net')    # net upward-directed flux

    if legend:
        ax1.legend(loc='center right', framealpha=1.0)
    ax1.set_yscale("symlog")
    ax1.set_ylabel(r"$F$ [W m$^{-2}$]")
    if logx:
        ax1.set_xscale("log")
    ax1.set_xticklabels([])
    if dx_tick > 1.0e-10:
        ax1.xaxis.set_minor_locator(MultipleLocator(dx_tick))

    yticks = ax1.get_yticks()
    yticks = np.delete(yticks, np.argwhere( yticks == 0.0 ))
    ax1.set_yticks(yticks)

    for p in x_planets.keys():
        for ax in (ax1,ax2):
            ax.axvline(x=x_planets[p], color=c_planets[p], lw=lw*0.5, linestyle='dashed', label=p , zorder=1)

    arr_magma = np.ones(len(r_arr))*T_magma
    ax2.plot(r_arr, np.zeros(len(r_arr)), zorder=6, color='silver', lw=lw, label=r"$\tilde{T}_s$") # magma temperature
    ax2.plot(r_arr, ts_arr - arr_magma,   zorder=6, color='black',  lw=lw, label=r"$T_s$") # surface solution temperature 

    if legend:
        ax2.legend(loc='center right', framealpha=1.0)
    ax2.set_ylabel(r"$T$ - $\tilde{T_s}$ [K]")
    if logx:
        ax2.set_xscale("log")
    if dx_tick > 1.0e-10:
        ax2.xaxis.set_minor_locator(MultipleLocator(dx_tick))

    ax2.set_xlabel("Orbital separation [AU]")
    fig.subplots_adjust(hspace=0.08)
    fig.savefig(dirs["output"]+"inst_%dK.pdf"%T_magma, bbox_inches='tight')
    

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    print("Done!")

