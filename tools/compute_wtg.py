import os
import numpy as np
import utils.GeneralAdiabat as ga
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import LogLocator
import seaborn as sns

aeolus_dir = os.getenv('AEOLUS_DIR')+"/"

plot_output_dir = aeolus_dir+'plots/'

if not os.path.exists(plot_output_dir):
        os.makedirs(plot_output_dir)
		
def wtg(rdgas,T,omega,radius):
    """Computing the Weak Temperature Gradient parameter of a planet.

    Parameters
    ----------
        rdgas : float
            Specific gas constant of the dominant species in the atmosphere [J.kg^-1.K^-1]
        T : float
            Characteristic temperature of the atmosphere [K]
        omega : float
            Angular velocity of the planet [rad.s^-1]
        radius : float
            Radius of the planet [m]
    """

    return np.sqrt(rdgas*T)/(omega*radius)

####################################
##### Global constants
####################################

# AU definition [m]
AU = 149597870700.
# Universal gravitational constant [m^3.kg^-1.s^-2]
G_universal = 6.67430e-11
# Universal gas constant [J.mol^-1.K^-1]
R = 8.31446261815324 
# Number of seconds in a day
day_s = 24.0*3600.0

####################################
##### Sun/Earth parameters
####################################

L_Sun        = 3.828e26    # Solar luminosity [W]
M_Sun        = 1.988500e30 # Solar mass [kg]
radius_Earth = 6.371e6     # Earth radius [m]
mass_Earth   = 5.9722e24   # Earth mass [kg]
grav_Earth   = 9.807       # Earth surface gravity [m.s^-2]

def surface_gravity(mass,radius):
	g = G_universal*np.array(mass)*mass_Earth/(np.array(radius)*radius_Earth)**2
	if g.ndim > 0:
		g = np.sort(g)
	return g

def angular_velocity(orbital_period):
	return 2.0*np.pi/(orbital_period*day_s)

def weighted_mean(mix_ratios,y):
	return sum(np.multiply(np.array(mix_ratios),np.array(y))) 

# Sources: TRAPPIST-1: doi:10.3847/psj/abd022 ; Proxima-b: doi:10.3847/2041-8205/831/2/l16 (Figure 2) ; Teegarden: doi:10.1051/0004-6361/201935460 (Figure 12). 
planets = { 'Earth': {'Stellar Mass': M_Sun, 'Stellar Luminosity': L_Sun, 'Radius': radius_Earth, 'Gravity': grav_Earth, 'Semi-major axis': AU, 'Angular velocity': angular_velocity(1.0)},
            'TRAPPIST-1b': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 1.116*radius_Earth, 'Gravity': 1.102*grav_Earth, 'Semi-major axis': 0.01154*AU, 'Angular velocity': angular_velocity(1.510826)},
			'TRAPPIST-1d': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 0.788*radius_Earth, 'Gravity': 0.624*grav_Earth, 'Semi-major axis': 0.02227*AU, 'Angular velocity': angular_velocity(4.049219)},
			'TRAPPIST-1e': {'Stellar Mass': 0.0898*M_Sun, 'Stellar Luminosity': 0.000553*L_Sun, 'Radius': 0.920*radius_Earth, 'Gravity': 0.817*grav_Earth, 'Semi-major axis': 0.02925*AU, 'Angular velocity': angular_velocity(6.101013)},
			'Proxima-b':   {'Stellar Mass': 0.1221*M_Sun, 'Stellar Luminosity': 0.001567*L_Sun, 'Radius': np.array([0.94,1.4])*radius_Earth, 'Gravity': surface_gravity(1.07,[0.94,1.4]), 'Semi-major axis': 0.04856*AU, 'Angular velocity': angular_velocity(11.1868)},
			'Teegarden-b':   {'Stellar Mass': 0.09455120*M_Sun, 'Stellar Luminosity': 0.00073*L_Sun, 'Radius': 1.02*radius_Earth, 'Mass': 1.05*mass_Earth, 'Gravity': surface_gravity(1.05,1.02), 'Semi-major axis': 0.0252*AU, 'Angular velocity': angular_velocity(4.9100)},
			'Teegarden-c':   {'Stellar Mass': 0.09455120*M_Sun, 'Stellar Luminosity': 0.00073*L_Sun, 'Radius': 1.04*radius_Earth, 'Mass': 1.11*mass_Earth, 'Gravity': surface_gravity(1.11,1.04), 'Semi-major axis': 0.0443*AU, 'Angular velocity': angular_velocity(11.409)}}

atmospheres = { 'Earth': {'H2O': 1e-3, 'CO2': 3.50e-4, 'O3': 0.07e-6, 'N2O': 0.31e-6, 'CO': 0.10e-6, 'CH4': 1.70e-6, 'O2': 0.20947, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.02e-6, 'NH3': 1.0e-7, 'HNO3': 0.0, 'N2': 0.78084, 'H2': 0.53e-6, 'He': 5.24e-6, 'OCS': 0.0},
                'Pure Steam': {'H2O': 1.0, 'CO2': 0.0, 'O3': 0.0, 'N2O': 0.0, 'CO': 0.0, 'CH4': 0.0, 'O2': 0.0, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.0, 'NH3': 0.0, 'HNO3': 0.0, 'N2': 0.0, 'H2': 0.0, 'He': 0.0, 'OCS': 0.0},
				'300K N2-H2O': {'H2O': ga.p_sat('H2O',300.0)/(1e5+ga.p_sat('H2O',300.0)), 'CO2': 0.0, 'O3': 0.0, 'N2O': 0.0, 'CO': 0.0, 'CH4': 0.0, 'O2': 0.0, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.0, 'NH3': 0.0, 'HNO3': 0.0, 'N2': 1e5/(1e5+ga.p_sat('H2O',300.0)), 'H2': 0.0, 'He': 0.0, 'OCS': 0.0},
				'Post-runaway N2-H2O': {'H2O': 260.0/261.0, 'CO2': 0.0, 'O3': 0.0, 'N2O': 0.0, 'CO': 0.0, 'CH4': 0.0, 'O2': 0.0, 'NO': 0.0, 'SO2': 0.0, 'NO2': 0.0, 'NH3': 0.0, 'HNO3': 0.0, 'N2': 1.0/261.0, 'H2': 0.0, 'He': 0.0, 'OCS': 0.0}}

molecules = {'Molecular Weight': {'H2O': 0.018, 'CO2': 0.044, 'O3': 0.048, 'N2O': 0.044, 'CO': 0.028, 'CH4': 0.016, 'O2': 0.032, 'NO': 0.030, 'SO2': 0.064, 'NO2': 0.046, 'NH3': 0.017, 'HNO3': 0.063, 'N2': 0.028, 'H2': 0.002, 'He': 0.004, 'OCS': 0.060},
			 'Isobaric Heat Capacity': {'H2O': 1864.0, 'CO2': 849.0, 'O3': 819.375, 'N2O': 877.364, 'CO': 1040.0, 'CH4': 2232.0, 'O2': 918.0, 'NO': 995.0, 'SO2': 624.0625, 'NO2': 805.0, 'NH3': 2175.0, 'HNO3': 849.365, 'N2': 1040.0, 'H2': 14310.0, 'He': 5197.5, 'OCS': 41.592}}

molecular_weights = list(molecules['Molecular Weight'].values())

print("psat 300K = ", ga.p_sat('H2O',300.0))
T_array = np.arange(200.0, 3000.0+50.0, 50)

ls_ind = 1.5
fs_l = 16
fs_m = 14
fs_s = 10

colors = cm.RdYlBu(np.linspace(0, 1, len(planets)+4)) 

fig, ax = plt.subplots(figsize=(8,8))

# ===== SELECT PLANET AND ATMOSPHERE =====

planet     = 'Earth'
atmosphere = 'Earth'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[0], label=planet)

planet     = 'Earth'
atmosphere = '300K N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="--", color=colors[1], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}=p_{\mathrm{sat}}$(300K) Pa')

planet     = 'Earth'
atmosphere = 'Post-runaway N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-.", color=colors[2], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}$=260 bar')

planet     = 'TRAPPIST-1b'
atmosphere = 'Pure Steam'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[3], label=planet+' - pure steam')

planet     = 'TRAPPIST-1d'
atmosphere = 'Pure Steam'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[4], label=planet+' - pure steam')

planet     = 'TRAPPIST-1e'
atmosphere = 'Earth'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color='y', label=planet+' - Earth composition')

# planet     = 'Proxima-b'
# atmosphere = 'Earth'
# atmospheric_composition = list(atmospheres[atmosphere].values())
# rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
# omega = planets[planet]['Angular velocity']
# radius = planets[planet]['Radius'] 
# wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
# ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[6], label=planet+' - Earth composition')

planet     = 'Teegarden-b'
atmosphere = '300K N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[7], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}=p_{\mathrm{sat}}$(300K) Pa')

planet     = 'Teegarden-b'
atmosphere = 'Post-runaway N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="--", color=colors[8], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}$=260 bar')

planet     = 'Teegarden-c'
atmosphere = '300K N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="-", color=colors[9], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}=p_{\mathrm{sat}}$(300K) Pa')

planet     = 'Teegarden-c'
atmosphere = 'Post-runaway N2-H2O'
atmospheric_composition = list(atmospheres[atmosphere].values())
rdgas = R/weighted_mean(atmospheric_composition,molecular_weights) 
omega = planets[planet]['Angular velocity']
radius = planets[planet]['Radius'] 
wtg_param = np.array([wtg(rdgas,T,omega,radius) for T in T_array])
ax.plot( T_array, wtg_param, lw=ls_ind, ls="--", color=colors[10], label=planet+r' $p_{\mathrm{N2}}$=1 bar, $p_{\mathrm{H2O}}$=260 bar')

ax.axhline(y=1.0, lw=0.5, ls="-", color='k')

ax.legend(fontsize=fs_s, loc='upper right')

ax.set_xlabel(r'Characteristic temperature, $T$ (K)', fontsize=fs_l)
ax.set_ylabel(r'WTG parameter, $\Lambda$', fontsize=fs_l)
ax.set_xlim([min(T_array),max(T_array)])
ax.set_yscale("log")
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)

ax.tick_params(axis='both', which='major', labelsize=fs_m)
ax.tick_params(axis='both', which='minor', labelsize=fs_m)
current_ticks = ax.get_xticks()
new_ticks = np.append(current_ticks, 0)
ax.set_xticks(new_ticks)

plt.show()

fig.savefig(plot_output_dir+'wtg.pdf', bbox_inches='tight')
fig.savefig(plot_output_dir+'wtg.png', bbox_inches='tight')
#plt.close()