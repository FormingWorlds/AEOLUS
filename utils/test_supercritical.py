import utils.GeneralAdiabat as ga
from utils.atmosphere_column import atmos
import copy

T_surf                  = 700.0                         # K
pH2O                    = ga.p_sat('H2O',T_surf)        # Pa
print("pH2O = ", pH2O)
pCO2                    = 0.                            # Pa
pH2                     = 0.                            # Pa
pN2                     = 1e+5                          # Pa
pCH4                    = 0.                            # Pa
pO2                     = 0.                            # Pa
pHe                     = 0.                            # Pa
pNH3                    = 0.                            # Pa
P_surf                  = pH2O + pCO2 + pH2 + pN2 + pCH4 + pO2 + pHe + pNH3  # Pa
print("P_surf = ", P_surf)

alpha_cloud             = 0.0

pl_radius     = 6.371e6             # m, planet radius
pl_mass       = 5.972e24            # kg, planet mass
P_top         = 1.0                 # Pa

re   = 1.0e-5 
lwm  = 0.8    
clfr = 0.0 

vol_partial = {}
vol_mixing = { 
                "H2O" : pH2O / P_surf,   
                "CO2" : pCO2 / P_surf,
                "H2"  : pH2  / P_surf,
                "N2"  : pN2  / P_surf,
                "CH4" : pCH4 / P_surf,
                "O2"  : pO2  / P_surf,
                "CO"  : pN2  / P_surf,
                "He"  : pHe  / P_surf,
                "NH3" : pNH3 / P_surf,
            }

# P_surf = 0.0
# vol_mixing = {}
# vol_partial = {
#     "H2O" : ga.p_sat('H2O',T_surf), #4e3,
#     "NH3" : 0.0,
#     "CO2" : 0.0,#4.0e2,
#     "CH4" : 0.0,
#     "CO" : 0.0,#3.0e2,
#     "O2" : 0.0,#0.0e4,
#     "N2" : 1.0e5,
#     "H2" : 0.0
#     }

atm = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, re, lwm, clfr,
                    vol_mixing=vol_mixing, vol_partial=vol_partial, calc_cf=False, trppT=False, req_levels=100, water_lookup=False, do_cloud=False)
#atm                     = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, re, lwm, clfr, vol_mixing=vol_mixing)

atm_moist = copy.deepcopy(atm)
    
atm_moist = ga.general_adiabat(atm_moist)