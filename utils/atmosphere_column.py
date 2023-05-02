# atmosphere_column.py
# class for atmospheric column data

import numpy as np

R_universal = 8.31446261815324 # Universal gas constant, J.K-1.mol-1

class atmos:
    '''
    Atmosphere class
    '''
    def __init__(self, T_surf, P_surf, vol_list, calc_cf=False, trppT=0):
        self.alpha_cloud 	= 0.0 	    	# The fraction of condensate retained in the column; 1 -> Li et al 2018; 0 -> full rainout

        # If vol_list is given in partial pressures, calculate mixing ratios
        if (type(P_surf) == str) or (type(P_surf) == float and P_surf <= 0.):
            P_surf          = sum(vol_list.values())
            print("Calculate mixing ratios from partial pressures.")
            print("P_surf:", P_surf)
            print("p_i:", vol_list)
            for vol in vol_list.keys():
                vol_list[vol] = vol_list[vol]/P_surf
            print("x_i:", vol_list)

        self.ps 			= P_surf 	 	# Surface pressure, Pa
        self.ts 			= T_surf		# Surface temperature, K
        self.vol_list 		= vol_list		# Names + mixing ratios dict

        self.ptop 			= 10 			# Top pressure in Pa
        self.nlev 			= 10000  	   	# Number of vertical levels for adiabat integration
        self.step    		= 0.01  		# Adjust to match self.nlev
        self.nlev_save		= 100   		# Number of levels to save object
        self.p 				= np.zeros(self.nlev) 	   		# np.ones(self.nlev)
        self.pl 			= np.zeros(self.nlev+1)    		# np.ones(self.nlev+1)

        self.trppidx		= 0 				 	# Tropopause: idx
        self.trppP 			= 0 				 	# Tropopause: prs

        # Nominal tropopause T in K; if set to zero dynamically calculated in SocRadConv.py
        self.trppT 			= 290 				 	# Tropopause: tmp
        
        self.dt 			= 0.5 							# days

        self.toa_heating    = 0. 							# W/m^2
        self.star_lum       = 0. 							# L_sun

        self.albedo_s   	= 0.0 							# surface albedo
        self.albedo_pl   	= 0.175 						# Bond albedo (scattering)
        self.zenith_angle  	= 54.55							# solar zenith angle, Hamano+15 (arccos(1/sqrt(3) = 54.74), Wordsworth+ 10: 48.19 (arccos(2/3)), see Cronin 14 (mu = 0.58 -> theta = arccos(0.58) = 54.55) for definitions

        self.planet_mass 	= 5.972e+24 					# kg
        self.planet_radius 	= 6.3781e+6 					# m
        self.grav_s 		= 6.67408e-11*self.planet_mass/(self.planet_radius**2) # m s-2
        
        self.tmp 			= np.zeros(self.nlev)      		# self.ts*np.ones(self.nlev)
        self.tmpl 			= np.zeros(self.nlev+1)
        self.Rcp    		= 2./7. 						# standard earth air
        self.n_species 		= 7
        self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])
        
        # self.bands 			= np.concatenate((np.arange(0,3000,20),np.arange(3000,9000,50),np.arange(9000,24500,500))) # cm
        self.bands 			= np.concatenate((np.arange(0,3000,25),np.arange(3000,11000,50),np.arange(11000,30500,500))) # cm, 318 bands: HITEMP-compatible spacing
        #self.bands 			= np.concatenate((np.arange(0,3000,500),np.arange(3000,11000,1000),np.arange(11000,30500,3000))) # cm, 20 bands
        
        self.band_centres 	= (self.bands[1:] + self.bands[:-1]) / 2
        self.band_widths 	= np.diff(self.bands)
        self.nbands 	    = np.size(self.bands)-1

        # Level-dependent quantities
        self.p_vol 			= {} # Gas phase partial pressures
        self.pl_vol 		= {} # Gas phase partial pressures
        self.x_gas 			= {} # Gas phase molar concentration
        self.x_cond         = {} # Condensed phase molar concentration
        self.grav_z			= np.zeros(self.nlev) # Local gravity
        self.z 				= np.zeros(self.nlev) # Atmospheric height
        self.mu 			= np.zeros(self.nlev) # Mean molar mass in level
        self.xd 			= np.zeros(self.nlev) # Molar concentration of dry gas
        self.xv 			= np.zeros(self.nlev) # Molar concentration of moist gas
        self.xc 			= np.zeros(self.nlev) # Molar concentration of condensed phase
        self.rho            = np.zeros(self.nlev) # Density of atmosphere at a given level
        self.ifatm 			= np.zeros(self.nlev) # Defines nth level to which atmosphere is calculated
        self.cp      		= np.zeros(self.nlev) # Mean heat capacity

        # Define T and P arrays from surface up
        self.tmp[0]         = T_surf         		# K
        self.p[0]           = P_surf         		# Pa
        self.z[0]           = 0         			# m
        self.grav_z[0]      = self.grav_s 			# m s-2


        # H2O floor to prevent NaNs
        self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-30 ] )

        # Instantiate object dicts and arrays
        for vol in self.vol_list.keys():
            # Instantiate as zero
            self.p_vol[vol]      = np.zeros(self.nlev)
            self.pl_vol[vol]     = np.zeros(self.nlev+1)
            self.x_gas[vol]      = np.zeros(self.nlev)
            self.x_cond[vol]     = np.zeros(self.nlev)
            
            #self.x_ocean[vol]    = 0.

            # Surface partial pressures
            self.p_vol[vol][0]   = self.ps * vol_list[vol]

        # Radiation heating and fluxes
        self.LW_flux_up 			= np.zeros(self.nlev)				# W/m^2
        self.LW_flux_down 			= np.zeros(self.nlev)				# W/m^2
        self.LW_flux_net			= np.zeros(self.nlev)				# W/m^2
        self.LW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.LW_heating				= np.zeros(self.nlev)				# K/day
        self.SW_flux_up 			= np.zeros(self.nlev)				# W/m^2
        self.SW_flux_down 			= np.zeros(self.nlev)				# W/m^2
        self.SW_flux_net			= np.zeros(self.nlev)				# W/m^2
        self.SW_spectral_flux_up 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.SW_heating				= np.zeros(self.nlev)				# K/day
        self.flux_up_total			= np.zeros(self.nlev)				# W/m^2
        self.flux_down_total		= np.zeros(self.nlev)				# W/m^2
        self.net_flux				= np.zeros(self.nlev)				# W/m^2
        self.net_spectral_flux	 	= np.zeros([self.nbands,self.nlev])	# W/m^2/(band)
        self.net_heating 			= np.zeros(self.nlev) 				# K/day
        # Contribution function arrays
        if calc_cf == True:
            self.cff 					= np.zeros(self.nlev) 				# normalised
            self.cff_i					= np.zeros([self.nbands,self.nlev]) # cf per band
            self.LW_flux_up_i 			= np.zeros([self.nbands,self.nlev])


        
