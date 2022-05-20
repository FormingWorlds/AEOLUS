'''
SocRadModel.py
Returns heating rates
MDH 25/01/19
'''

import os, glob, re, shutil
import netCDF4 as net
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import subprocess
import nctools
import RayleighSpectrum
from subprocess import call
from netCDF4 import Dataset

try:
    from atmosphere_column import atmos
except:
    from atm_rad_conv.atmosphere_column import atmos
# from natsort import natsorted # https://pypi.python.org/pypi/natsort

def radCompSoc(atm, dirs, recalc, calc_cf=False, rscatter=False):

    # Enable or disable calculating contribution function
    # ENABLE RIGHT ENVIRONMENT IN TERMINAL FIRST
    # calc_cf = False
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

    # Define stellar spectrum to use, options:
    # Sun_t456Myr_kurucz_95 F2V_hd128167 M45_ADLeo Sun_t0_0Ga_claire_12 (t: 0.0 – 4.55)
    star_name           = "Sun_t0_0Ga_claire_12"
    
    # Define path to spectral file
    spectral_base_name  = "sp_b318_HITRAN_a16"
    spectral_name       = "sp_b318_HITRAN_a16"+"_"+star_name
    spectral_dir        = dirs["rad_conv"]+"/spectral_files/"+spectral_base_name+"/"
    spectral_file       = spectral_dir+spectral_name

    # Rayleigh scattering for CO2
    if rscatter == True:

        # Generate new temp directory for scattering spectral file
        scatter_dir = dirs["output"]+"/"+"scatter_spectral_file/"
        if not os.path.isdir(scatter_dir):
            print(">>> Generate scatter spectral file dir:", scatter_dir)
            os.makedirs(scatter_dir)

        # New file
        spectral_file = scatter_dir+spectral_name

        # Copy if not there yet
        if not os.path.isfile(spectral_file):
            print(">>> Copy non-scatter spectral file to:", scatter_dir)
            for file in natural_sort(glob.glob(spectral_dir+"*")):
                shutil.copy(file, scatter_dir+os.path.basename(file))
                print(os.path.basename(file), end =" ")

        # Insert Rayleigh scattering into spectral file
        RayleighSpectrum.rayleigh_coeff_adder(species_list = ['co2', 'h2o', 'n2'], mixing_ratio_list = [atm.x_gas["CO2"][-1], atm.x_gas["H2O"][-1], atm.x_gas["N2"][-1]], spectral_file_path=spectral_file,wavelength_dummy_file_path=scatter_dir+'wavelength_band_file.txt')


    # # Enable cf SOCRATES environment
    # if calc_cf == True:
    #     path_to_cff_socrates = dirs["rad_conv"]+"/rad_trans/socrates_cff"
    #     call_sequence = [ "source", path_to_cff_socrates+"/set_rad_env" ]
    #     print(call_sequence)
    #     subprocess.call(call_sequence)

    # Remove auxiliary files from previous runs
    CleanOutputDir( os.getcwd() )

    # Solar zenith angle
    zenith_angle    = atm.zenith_angle  # Hamano+15: 54.7, Ranjan+18: 48.2, Katyal+19: 38

    # Surface albedo
    surface_albedo  = atm.albedo_s      # Hamano+15: 0.2, Schaefer+16: 0.75

    # Other parameters
    longitude       = 0
    latitude        = 0
    basis_function  = 1

    # Write values to netcdf: SOCRATES Userguide p. 45
    blockPrint()
    nctools.ncout_surf('profile.surf', longitude, latitude, basis_function, surface_albedo)
    nctools.ncout2d('profile.tstar', 0, 0, atm.ts, 'tstar', longname="Surface Temperature", units='K')
    nctools.ncout2d('profile.pstar', 0, 0, atm.ps, 'pstar', longname="Surface Pressure", units='PA')
    nctools.ncout2d('profile.szen', 0, 0, zenith_angle, 'szen', longname="Solar zenith angle", units='Degrees')
    nctools.ncout2d('profile.stoa', 0, 0, atm.toa_heating, 'stoa', longname="Solar Irradiance at TOA", units='WM-2')
    # T, P + volatiles
    nctools.ncout3d('profile.t', 0, 0,   atm.p,  atm.tmp, 't', longname="Temperature", units='K')
    nctools.ncout3d('profile.tl', 0, 0,  atm.pl, atm.tmpl, 'tl', longname="Temperature", units='K')
    nctools.ncout3d('profile.p', 0, 0,   atm.p,  atm.p, 'p', longname="Pressure", units='PA')
    nctools.ncout3d('profile.q', 0, 0,   atm.p,  molar_mass['H2O'] / atm.mu * atm.x_gas["H2O"], 'q', longname="q", units='kg/kg') 
    if "CO2" in atm.vol_list.keys():
        nctools.ncout3d('profile.co2', 0, 0, atm.p,  molar_mass['CO2'] / atm.mu * atm.x_gas["CO2"], 'co2', longname="CO2", units='kg/kg') 
    if "O3" in atm.vol_list.keys():
        nctools.ncout3d('profile.o3', 0, 0,  atm.p,  molar_mass['O3'] / atm.mu * atm.x_gas["O3"], 'o3', longname="O3", units='kg/kg') 
    if "N2O" in atm.vol_list.keys():
        nctools.ncout3d('profile.n2o', 0, 0,  atm.p,  molar_mass['N2O'] / atm.mu * atm.x_gas["N2O"], 'n2o', longname="N2O", units='kg/kg') 
    if "CO" in atm.vol_list.keys():
        nctools.ncout3d('profile.co', 0, 0,  atm.p,  molar_mass['CO'] / atm.mu * atm.x_gas["CO"], 'co', longname="CO", units='kg/kg') 
    if "CH4" in atm.vol_list.keys():
        nctools.ncout3d('profile.ch4', 0, 0, atm.p,  molar_mass['CH4'] / atm.mu * atm.x_gas["CH4"], 'ch4', longname="CH4", units='kg/kg') 
    if "O2" in atm.vol_list.keys():
        nctools.ncout3d('profile.o2', 0, 0,  atm.p,  molar_mass['O2'] / atm.mu * atm.x_gas["O2"], 'o2', longname="O2", units='kg/kg')
    if "NO" in atm.vol_list.keys():
        nctools.ncout3d('profile.no', 0, 0,  atm.p,  molar_mass['NO'] / atm.mu * atm.x_gas["NO"], 'no', longname="NO", units='kg/kg') 
    if "SO2" in atm.vol_list.keys():
        nctools.ncout3d('profile.so2', 0, 0,  atm.p,  molar_mass['SO2'] / atm.mu * atm.x_gas["SO2"], 'so2', longname="SO2", units='kg/kg') 
    if "NO2" in atm.vol_list.keys():
        nctools.ncout3d('profile.no2', 0, 0,  atm.p,  molar_mass['NO2'] / atm.mu * atm.x_gas["NO2"], 'no2', longname="NO2", units='kg/kg') 
    if "NH3" in atm.vol_list.keys():
        nctools.ncout3d('profile.nh3', 0, 0,  atm.p,  molar_mass['NH3'] / atm.mu * atm.x_gas["NH3"], 'nh3', longname="NH3", units='kg/kg') 
    if "HNO3" in atm.vol_list.keys():
        nctools.ncout3d('profile.hno3', 0, 0,  atm.p,  molar_mass['HNO3'] / atm.mu * atm.x_gas["HNO3"], 'hno3', longname="HNO3", units='kg/kg')
    if "N2" in atm.vol_list.keys():
        nctools.ncout3d('profile.n2', 0, 0,  atm.p,  molar_mass['N2'] / atm.mu * atm.x_gas["N2"], 'n2', longname="N2", units='kg/kg') 
    if "H2" in atm.vol_list.keys():
        nctools.ncout3d('profile.h2', 0, 0,  atm.p,  molar_mass['H2'] / atm.mu * atm.x_gas["H2"], 'h2', longname="H2", units='kg/kg')
    if "He" in atm.vol_list.keys():
        nctools.ncout3d('profile.he', 0, 0,  atm.p,  molar_mass['He'] / atm.mu * atm.x_gas["He"], 'he', longname="He", units='kg/kg')
    if "OCS" in atm.vol_list.keys():
        nctools.ncout3d('profile.ocs', 0, 0,  atm.p,  molar_mass['OCS'] / atm.mu * atm.x_gas["OCS"], 'ocs', longname="OCS", units='kg/kg')
    
    enablePrint()

    basename = 'profile'
    s = " "

    if rscatter == True:
        scatter_flag = " -r"
    else:
        scatter_flag = ""

    # Anchor spectral files and run SOCRATES
    seq4 = ("Cl_run_cdf -B", basename,"-s", spectral_file, "-R 1", str(atm.nbands), " -ch ", str(atm.nbands), " -S -g 2 -C 5 -u", scatter_flag)
    seq5 = ("fmove", basename,"currentsw")
    seq6 = ("Cl_run_cdf -B", basename,"-s", spectral_file, "-R 1 ", str(atm.nbands), " -ch ", str(atm.nbands), " -I -g 2 -C 5 -u", scatter_flag)
    seq7 = ("fmove", basename,"currentlw")

    if calc_cf == True:
        seq8 = ("Cl_run_cdf -B", basename,"-s", spectral_file, "-R 1 ", str(atm.nbands), " -ch ", str(atm.nbands), " -I -g 2 -C 5 -u -ch 1", scatter_flag)
        seq9 = ("fmove", basename, "currentlw_cff")

    comline1 = s.join(seq4)
    comline2 = s.join(seq5)
    comline3 = s.join(seq6)
    comline4 = s.join(seq7)
    if calc_cf == True:
        comline5 = s.join(seq8)
        comline6 = s.join(seq9)

    os.system(comline1)
    os.system(comline2)
    os.system(comline3)
    os.system(comline4)
    if calc_cf == True:
        os.system(comline5)
        os.system(comline6)


    # Open netCDF files produced by SOCRATES
    ncfile1  = net.Dataset('currentsw.vflx')
    ncfile2  = net.Dataset('currentsw.sflx')
    ncfile3  = net.Dataset('currentsw.dflx')
    ncfile4  = net.Dataset('currentsw.uflx')
    ncfile5  = net.Dataset('currentsw.nflx')
    ncfile6  = net.Dataset('currentsw.hrts')
    ncfile7  = net.Dataset('currentlw.dflx')
    ncfile8  = net.Dataset('currentlw.nflx')
    ncfile9  = net.Dataset('currentlw.uflx')
    ncfile10 = net.Dataset('currentlw.hrts')
    if calc_cf == True:
        ncfile11 = net.Dataset('currentlw_cff.cff')
        ncfile12 = net.Dataset('currentlw.cff')


    # Loop through netCDF variables and populate arrays
    vflxsw   = ncfile1.variables['vflx']  # SW downward flux (direct + diffuse)
    uflxsw   = ncfile4.variables['uflx']  # SW upward flux 
    nflxsw   = ncfile5.variables['nflx']  # SW net flux 
    hrtssw   = ncfile6.variables['hrts']  # SW heating rate (K/day)
    dflxlw   = ncfile7.variables['dflx']  # LW downward flux (diffuse)
    uflxlw   = ncfile9.variables['uflx']  # LW upward flux 
    nflxlw   = ncfile8.variables['nflx']  # LW net flux 
    hrtslw   = ncfile10.variables['hrts'] # LW heating rate (K/day)
    if calc_cf == True:
        cff    = ncfile11.variables['cff'] # Contribution function (flux, sum)
        cff_i  = ncfile12.variables['cff'] # Contribution function (channel, plev, lat, lon)


    ##### Fluxes

    if calc_cf == False:

        # Upward SW + LW flux summed over all bands (W/m^2)
        atm.SW_flux_up          = np.sum(uflxsw[:,:],axis=0)[:,0,0]
        atm.LW_flux_up          = np.sum(uflxlw[:,:],axis=0)[:,0,0]

        # Downward SW + LW flux summed over all bands (W/m^2)
        atm.SW_flux_down        = np.sum(vflxsw[:,:],axis=0)[:,0,0]
        atm.LW_flux_down        = np.sum(dflxlw[:,:],axis=0)[:,0,0]

        # Net SW + LW flux summed over all bands (W/m^2)
        atm.SW_flux_net         = np.squeeze(np.sum(uflxsw[:,:],axis=0)[:,0,0] - np.sum(vflxsw[:,:],axis=0)[:,0,0])
        atm.LW_flux_net         = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] - np.sum(dflxlw[:,:],axis=0)[:,0,0])

        # Upward SW + LW flux per band, W/m^2/(band)
        atm.LW_spectral_flux_up = uflxlw[:,:,0,0]
        atm.SW_spectral_flux_up = uflxsw[:,:,0,0]

        # Total up- and downward fluxes, (W/m^2)
        atm.flux_up_total       = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0]) 
        atm.flux_down_total     = np.squeeze(np.sum(vflxsw[:,:],axis=0)[:,0,0] + np.sum(dflxlw[:,:],axis=0)[:,0,0])

        # Total net flux (W/m^2)
        atm.net_flux            = np.squeeze(np.sum(uflxlw[:,:],axis=0)[:,0,0] - np.sum(dflxlw[:,:],axis=0)[:,0,0] + np.sum(uflxsw[:,:],axis=0)[:,0,0] -  np.sum(vflxsw[:,:],axis=0)[:,0,0])

        # Total net flux per band (W/m^2/(band))
        atm.net_spectral_flux   = uflxlw[:,:,0,0] + uflxsw[:,:,0,0] - dflxlw[:,:,0,0] - vflxsw[:,:,0,0]

        ##### Heating rates

        # Heating rates only for no recalc: recalc if tropopause was found
        if recalc == False:

            # Individual heating contributions (K/day)
            atm.SW_heating          = np.sum(hrtssw[:,:],axis=0)[:,0,0]
            atm.LW_heating          = np.sum(hrtslw[:,:],axis=0)[:,0,0]

            # Total heating (K/day)
            atm.net_heating       = np.squeeze(np.sum(hrtssw[:,:],axis=0) + np.sum(hrtslw[:,:],axis=0))
    
    ### Contribution function
    if calc_cf == True:
        # print(cff)
        atm.cff     = cff[:,0,0]
        atm.cff_i   = cff_i[:,:,0,0]
        atm.LW_flux_up_i = uflxlw[:,:,0,0]

    # Close netCDF files
    ncfile1.close()
    ncfile2.close()
    ncfile3.close()
    ncfile4.close()
    ncfile5.close()
    ncfile6.close()
    ncfile7.close()
    ncfile8.close()
    ncfile9.close()
    ncfile10.close()
    if calc_cf == True:
        ncfile11.close()
        # ncfile12.close()

    # Remove auxiliary files
    CleanOutputDir( os.getcwd() )

    return atm

# Sting sorting not based on natsorted package
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# Clean SOCRATES output files after run
def CleanOutputDir( output_dir ):

    types = ("current*.*", "profile.*") 
    files_to_delete = []
    for files in types:
        files_to_delete.extend(glob.glob(output_dir+"/"+files))

    # print("Remove old output files:")
    for file in natural_sort(files_to_delete):
        os.remove(file)
    #     print(os.path.basename(file), end =" ")
    # print("\n==> Done.")

# Disable and enable print: https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__
