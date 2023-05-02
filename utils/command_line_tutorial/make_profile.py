# This standalone code generates netCDF input files for SOCRATES including a modifiable PT profile corresponding to a pure steam adiabatic column. 
# It requires the SOCRATES environment to access the module nctools from socrates_main/python/nctools.py 

import numpy as np
import nctools
import os
import sys

# Disable and enable print: https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

# Solar zenith angle
zenith_angle    = 54.7  # Hamano+15: 54.7, Ranjan+18: 48.2, Katyal+19: 38

# Surface albedo
surface_albedo  = 0.2   # Hamano+15: 0.2, Schaefer+16: 0.75

# Stellar irradiance at the TOA
AU = 149597870700.
L_Sun = 3.828e26

# Semi-major axis
# Trappist-1b: 0.01154
# Trappist-1d: 0.02227

sigma = 5.670374419e-8
a_major = 0.01154
I0 = (1.0 - 0.0) * 0.000553*L_Sun / ( 4. * np.pi * (a_major*AU)**2. )

# Other parameters
longitude       = 0
latitude        = 0
basis_function  = 1

# Defining the hybrid-sigma pressure array
# Units: pref is in Pa, ak is in Pa, and bk is nondimensional (Pa/pref)

a47 = np.array([ 10.00000,    24.45365,    48.76776,    \
                 85.39458,    133.41983,   191.01402,   \
                 257.94919,   336.63306,   431.52741,   \
                 548.18995,   692.78825,   872.16512,   \
                 1094.18467,  1368.11917,  1704.99489,  \
                 2117.91945,  2622.42986,  3236.88281,  \
                 3982.89623,  4885.84733,  5975.43260,  \
                 7286.29500,  8858.72424,  10739.43477, \
                 12982.41110, 15649.68745, 18811.37629, \
                 22542.71275, 25724.93857, 27314.36781, \
                 27498.59474, 26501.79312, 24605.92991, \
                 22130.51655, 19381.30274, 16601.56419, \
                 13952.53231, 11522.93244, 9350.82303,  \
                 7443.47723,  5790.77434,  4373.32696,  \
                 3167.47008,  2148.51663,  1293.15510,  \
                 581.62429,   0.00000,     0.00000 ])

b47 = np.array([ 0.0000,        0.0000,        0.0000,  \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.00000,       0.00000, \
                 0.00000,       0.01188,       0.04650, \
                 0.10170,       0.17401,       0.25832, \
                 0.34850,       0.43872,       0.52448, \
                 0.60307,       0.67328,       0.73492, \
                 0.78834,       0.83418,       0.87320, \
                 0.90622,       0.93399,       0.95723, \
                 0.97650,       0.99223,       1.00000 ])

# Define the slope exponent
k = 1.38065812e-23            # Boltzman thermodynamic constant
N_avogadro = 6.022136736e+23
Rstar = 1000.*k*N_avogadro    # Universal gas constant
H2O_MolecularWeight = 1.8e+01
R = Rstar/H2O_MolecularWeight # J/(kg*K) specific gas constant of water vapor
H2O_cp = 1.847000e+03
Rcp = R/H2O_cp

psurf = 1e6
Tsurf = 1225.0
#Interface levels
pint = a47 + b47*psurf
Tint = Tsurf*(pint/psurf)**Rcp
#Midpoint levels
pmid = 0.5*(pint[1:len(a47)]+pint[0:len(a47)-1])
Tmid = Tsurf*(pmid/psurf)**Rcp # Initialize on an adiabat
#Tmid = Tsurf*(pmid/pmid)      # Initialize on an isotherm

# Write values to netcdf: SOCRATES Userguide p. 45
blockPrint()
nctools.ncout_surf('profile.surf', longitude, latitude, basis_function, surface_albedo)
nctools.ncout2d('profile.tstar', 0, 0, Tsurf, 'tstar', longname="Surface Temperature", units='K')
nctools.ncout2d('profile.pstar', 0, 0, psurf, 'pstar', longname="Surface Pressure", units='PA')
nctools.ncout2d('profile.szen', 0, 0, zenith_angle, 'szen', longname="Solar zenith angle", units='Degrees')
nctools.ncout2d('profile.stoa', 0, 0, I0, 'stoa', longname="Solar Irradiance at TOA", units='WM-2')
# T, P + volatiles
nctools.ncout3d('profile.t', 0, 0,   pmid,  Tmid, 't', longname="Temperature", units='K')
nctools.ncout3d('profile.tl', 0, 0,  pint, Tint, 'tl', longname="Temperature", units='K')
nctools.ncout3d('profile.p', 0, 0,   pmid,  pmid, 'p', longname="Pressure", units='PA')
nctools.ncout3d('profile.q', 0, 0,   pmid, 1.0, 'q', longname="q", units='kg/kg') # q = 1 if pure steam
enablePrint()

basename = 'profile'
