# atmosphere_column.py
# class for atmospheric column data

import numpy as np
import netCDF4 as nc
from utils import phys
from utils.height import AtmosphericHeight
import os

class atmos:
    
    def __init__(self, T_surf: float, P_surf: float, P_top: float, pl_radius: float, pl_mass: float, re: float, lwm: float, clfr: float,
                 vol_mixing: dict = {}, vol_partial: dict = {}, 
                 calc_cf: bool=False, req_levels: int = 100, water_lookup: bool=False,
                 trppT: float = 290.0, minT: float = 1.0, maxT: float = 9000.0, do_cloud: bool=False):
        
        """Atmosphere class    
    
        Stores compositional and thermodynamic information for the column.    
        Also stores planetary parameters.   

        One of either vol_mixing or vol_partial must be passed in. 
        If vol_partial is passed, then the value of P_surf is recalculated as the sum of partial pressures.

        Parameters
        ----------
            T_surf : float
                Surface temperature [K]
            P_surf : float
                Surface pressure [Pa]
            P_top : float
                Pressure at top of column [Pa]
            pl_radius : float
                Radius of rocky part of planet [m]
            pl_mass : float
                Mass of rocky part of planet [kg]
            re : float
                Effective radius of cloud droplets [m]
            lwm : float
                Liquid water mass fraction [kg/kg]
            clfr : float
                Water cloud fraction [adimensional]

            vol_mixing : dict
                Dictionary of volatiles (keys) and mixing ratios (values)
            vol_partial: dict
                Dictionary of volatiles (keys) and partial pressures (values)

            calc_cf : bool
                Calculate contribution function?
            req_levels : int
                Requested number of vertical levels
            water_lookup : bool
                Use lookup table for water thermodynamic values (e.g. L, c_p)
            trppT : float
                Tropopause temperature
            minT : float
                Temperature floor
            maxT : float
                Temperature ceiling
                
        """

        # Parse volatiles 
        if (len(vol_mixing) == 0) and (len(vol_partial) == 0):
            raise Exception("Either vol_mixing OR vol_partial must be passed to atmos.__init__ function!\nNeither were.")
        if (len(vol_mixing) > 0) and (len(vol_partial) > 0):
            raise Exception("Either vol_mixing OR vol_partial must be passed to atmos.__init__ function!\nBoth were.")
        
        if len(vol_mixing) > 0:  # Set by mixing ratio
            self.ps = P_surf
            if (P_surf <= 0.0):
                raise Exception("Surface pressure passed to atmos.__init__ must be positive!\nValue passed = '%g'" % P_surf)
            
            tot_mixing =  float(sum(vol_mixing.values()))  # Ensure mixing ratios add up to unity
            self.vol_list = {}
            for key in vol_mixing.keys():
                self.vol_list[key] = vol_mixing[key]/tot_mixing

        
        if len(vol_partial) > 0: # Set by partial pressure
            self.ps = float(sum(vol_partial.values()))
            self.vol_list = {}
            for key in vol_partial.keys():
                self.vol_list[key] = vol_partial[key]/self.ps

        # Required volatiles
        required_vols = {"H2O","CO2","N2"}
        if len(required_vols.intersection(self.vol_list.keys())) < len(required_vols):
            raise Exception("Missing required volatiles!\nRequired vols = %s" % str(required_vols))

        # H2O floor to prevent NaNs
        self.vol_list["H2O"] = np.max( [ self.vol_list["H2O"], 1e-30 ] )

        # Initialise other variables
        self.alpha_cloud 	= 0.0 	    	# The fraction of condensate retained in the column; 1 -> Li et al 2018; 0 -> full rainout
        
        self.ts 			= T_surf		# Surface temperature, K

        self.ptop 			= P_top 			# Top pressure in Pa
        self.nlev 			= 10000  	   	# Number of vertical levels for adiabat integration
        self.step    		= 0.01  		# Adjust to match self.nlev
        self.nlev_save		= int(max(req_levels, 10))  		# Number of levels to save object
        self.p 				= np.zeros(self.nlev) 	   		# np.ones(self.nlev)
        self.pl 			= np.zeros(self.nlev+1)    		# np.ones(self.nlev+1)

        self.trppT          = trppT                 # Fixed value [K]
        self.trppidx		= 0 				 	# Tropopause: idx
        self.trppP 			= 0 				 	# Tropopause: prs
        self.minT           = minT                  # Minimum temperature allowed [K]
        self.maxT           = maxT                  # Maximum ^

        self.dt 			= 0.5 							# days

        self.toa_heating    = 0. 							# W/m^2
        self.star_lum       = 0.0							# L_sun

        self.albedo_s   	= 0.0 							# surface albedo
        self.albedo_pl   	= 0.175 						# Bond albedo (scattering)
        self.zenith_angle  	= 54.74							# solar zenith angle, Hamano+15 (arccos(1/sqrt(3) = 54.74), Wordsworth+10: 48.19 (arccos(2/3)), see Cronin+14 for definitions

        self.planet_mass    = pl_mass
        self.planet_radius  = pl_radius
        self.grav_s 		= phys.G*self.planet_mass/(self.planet_radius**2) # m s-2
        
        self.tmp 			= np.zeros(self.nlev)      		# self.ts*np.ones(self.nlev)
        self.tmpl 			= np.zeros(self.nlev+1)
        self.Rcp    		= 2./7. 						# standard earth air
        self.n_species 		= 7
        self.mixing_ratios 	= np.zeros([self.n_species,self.nlev])

        self.water_lookup   = water_lookup
        
        # self.bands 			= np.concatenate((np.arange(0,3000,20),np.arange(3000,9000,50),np.arange(9000,24500,500))) # cm
        self.bands 			= np.concatenate((np.arange(0,3000,25),np.arange(3000,11000,50),np.arange(11000,30500,500))) # cm, 318 bands: HITEMP-compatible spacing
        
        self.band_centres 	= (self.bands[1:] + self.bands[:-1]) / 2
        self.band_widths 	= np.diff(self.bands)
        self.nbands 	    = np.size(self.bands)-1

        self.tmp_magma      = 3000.0
        self.skin_d         = 0.01 # m
        self.skin_k         = 2.0  # W m-1 K-1

        # Level-dependent quantities
        self.p_vol 			= {} # Gas phase partial pressures
        self.pl_vol 		= {} # Gas phase partial pressures
        self.x_gas 			= {} # Gas phase molar concentration
        self.x_cond         = {} # Condensed phase molar concentration
        self.grav_z			= np.zeros(self.nlev) # Local gravity
        self.z 				= np.zeros(self.nlev)   # Atmospheric height (centres)
        self.zl 	    	= np.zeros(self.nlev+1) # Atmospheric height (edges)
        self.mu 			= np.zeros(self.nlev) # Mean molar mass in level
        self.xd 			= np.zeros(self.nlev) # Molar concentration of dry gas
        self.xv 			= np.zeros(self.nlev) # Molar concentration of moist gas
        self.xc 			= np.zeros(self.nlev) # Molar concentration of condensed phase
        self.rho            = np.zeros(self.nlev) # Density of atmosphere at a given level
        self.ifatm 			= np.zeros(self.nlev) # Defines nth level to which atmosphere is calculated
        self.cp      		= np.zeros(self.nlev) # Mean heat capacity
        self.rh      		= np.zeros(self.nlev) # Relative humidity

        # Define T and P arrays from surface up
        self.tmp[0]         = self.ts         		# K
        self.tmpl[0]         = self.ts         		# K
        self.p[0]            = self.ps         		# Pa
        self.pl[0]           = self.ps         		# Pa
        self.z[0]           = 0         			# m
        self.zl[0]
        self.grav_z[0]      = self.grav_s 			# m s-2


        # Instantiate object dicts and arrays
        for vol in self.vol_list.keys():
            # Instantiate as zero
            self.p_vol[vol]      = np.zeros(self.nlev)
            self.pl_vol[vol]     = np.zeros(self.nlev+1)
            self.x_gas[vol]      = np.zeros(self.nlev)
            self.x_cond[vol]     = np.zeros(self.nlev)
        
            # Surface partial pressures
            self.p_vol[vol][0]   = self.ps * self.vol_list[vol]

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
        self.net_heating 			= np.zeros(self.nlev) 				# K/day from socrates
        self.heat                   = np.zeros(self.nlev)               # K/day from *

        # Contribution function arrays
        if calc_cf == True:
            self.cff 					= np.zeros(self.nlev) 				# normalised
            self.cff_i					= np.zeros([self.nbands,self.nlev]) # cf per band
            self.LW_flux_up_i 			= np.zeros([self.nbands,self.nlev])

        # Cloud flags (socrates/bin/rad_pcf.f90) and input arrays
        if do_cloud == True:            
            self.cloud_scheme 			= 2	# -C 2 = ip_cloud_mix_max : maximum/random overlap in a mixed column 
            self.cloud_representation   = 1 # -K 1 = ip_cloud_homogen : ice and water mixed homogeneously 
            self.droplet_type           = 5 # -d 5 
            self.solver                 = 16 # -v 16 = i_solver : solver used for the two-stream calculations, chosen from those defined in solver_pcf.f90.
            # The below three variables should be zero until a cloud scheme forms clouds
            self.re                     = np.zeros(self.nlev_save) # Effective radius of the droplets [m]
            self.lwm                    = np.zeros(self.nlev_save) # Liquid water mass fraction [kg/kg]
            self.clfr                   = np.zeros(self.nlev_save) # Water cloud fraction
            # Floats defined in SocRadConv and used in the cloud scheme
            self.effective_radius       = re
            self.liquid_water_fraction  = lwm
            self.cloud_fraction         = clfr
        else:
            self.cloud_scheme 			= 5 # -C 5 = ip_cloud_off : clear sky. Other flags will be ignored. In principle we shouldn't need to define the variables below in this case.
            self.cloud_representation   = 1 
            self.droplet_type           = 5  
            self.solver                 = 13
            self.re                     = np.zeros(self.nlev_save)
            self.lwm                    = np.zeros(self.nlev_save)
            self.clfr                   = np.zeros(self.nlev_save) 
            self.effective_radius       = 0.0
            self.liquid_water_fraction  = 0.0
            self.cloud_fraction         = 0.0



    def write_PT(self,filename: str="output/PT.tsv", punit:str = "Pa"):
        """Write PT profile to file, with descending pressure.

        Useful for when run outside of PROTEUS and for debugging.
        Designed with VULCAN compatibility in mind.

        Parameters
        ----------
            filename : string
                Output filename
            punit : string
                Pressure unit to use. Options: Pa, bar, dyne/cm2, atm.

        """

        p_scalefactor = 1.0
        match punit:
            case "Pa":
                p_scalefactor = 1.0
            case "bar":
                p_scalefactor = 1.0e-5
            case "dyne/cm2":
                p_scalefactor = 1.0e5
            case "atm":
                p_scalefactor = 1.01325e5
            case _:
                raise Exception("Unrecognised pressure unit '%s'"%punit)

        p_save = np.array(self.p) * p_scalefactor
        T_save = np.array(self.tmp)

        X = np.array([p_save,T_save]).T[::-1]

        header = '# (%s)\t(K)\nPressure\tTemp' % punit

        np.savetxt(filename,X,fmt='%1.5e',header=header,comments='',delimiter='\t')


    def write_ncdf(self, fpath):
        """Write atmosphere arrays to a netCDF file

        Parameters
        ----------
            fpath : string
                Output filename
        """ 

        # ----------------------
        # Calculate gravity and height (in case it hasn't been done already)
        self.z = AtmosphericHeight(self, self.planet_mass, self.planet_radius)

        self.zl = np.zeros(len(self.z)+1)
        for i in range(1,len(self.z)):
            self.zl[i] = 0.5 * (self.z[i-1] + self.z[i])
        self.zl[0] = 2*self.z[0] - self.zl[1] # estimate TOA height

        # ----------------------
        # Prepare NetCDF
        if os.path.exists(fpath):
            os.remove(fpath)

        ds = nc.Dataset(fpath, 'w', format='NETCDF4')
        ds.description = 'AEOLUS atmosphere data'

        nlev_c = len(self.p)
        nlev_l = nlev_c + 1

        gas_list = [str(gas) for gas in self.vol_list.keys()]
        ngases = len(gas_list)

        nchars = 16

        # ----------------------
        # Create dimensions
        ds.createDimension('nlev_c', nlev_c)      # Cell centres
        ds.createDimension('nlev_l', nlev_l)      # Cell edges
        ds.createDimension('ngases', ngases)      # Gases
        ds.createDimension('nchars', nchars)      # Length of string containing gas names
        ds.createDimension('nbands', self.nbands) # Number of bands in the spectral file

        # ----------------------
        # Scalar quantities  
        #    Create variables
        var_tstar =     ds.createVariable('tstar',      'f4');  var_tstar.units = "K"       # BOA LW BC
        var_toah =      ds.createVariable('toa_heating','f4');  var_toah.units = "W m-2"    # TOA SW BC
        var_ps =        ds.createVariable('ps',         'f4');  var_ps.units = "Pa"
        var_ptop =      ds.createVariable('ptop',       'f4');  var_ptop.units = "Pa" 
        var_trppP =     ds.createVariable('trppP',      'f4');  var_trppP.units = "Pa" 

        #     Store data
        var_tstar.assignValue(self.ts)
        var_toah.assignValue(self.toa_heating)
        var_ps.assignValue(self.ps)
        var_ptop.assignValue(self.ptop)
        var_trppP.assignValue(self.trppP)

        # ----------------------
        # Layer quantities  
        #    Create variables
        var_p =         ds.createVariable('p',       'f4', dimensions=('nlev_c'));           var_p.units = "Pa"
        var_pl =        ds.createVariable('pl',      'f4', dimensions=('nlev_l'));           var_pl.units = "Pa"
        var_tmp =       ds.createVariable('tmp',     'f4', dimensions=('nlev_c'));           var_tmp.units = "K"
        var_tmpl =      ds.createVariable('tmpl',    'f4', dimensions=('nlev_l'));           var_tmpl.units = "K"
        var_z =         ds.createVariable('z',       'f4', dimensions=('nlev_c'));           var_z.units = "m"
        var_zl =        ds.createVariable('zl',      'f4', dimensions=('nlev_l'));           var_zl.units = "m"
        var_grav =      ds.createVariable('gravity', 'f4', dimensions=('nlev_c'));           var_grav.units = "m s-2"
        var_cp =        ds.createVariable('cp',      'f4', dimensions=('nlev_c'));           var_cp.units = "J/(kg K)"
        var_xd =        ds.createVariable('xd',      'f4', dimensions=('nlev_c'));           var_xd.units = "none"
        var_xv =        ds.createVariable('xv',      'f4', dimensions=('nlev_c'));           var_xv.units = "none"

        var_mmw =       ds.createVariable('mmw',     'f4', dimensions=('nlev_c'));           var_mmw.units = "kg mol-1"
        var_gases =     ds.createVariable('gases',   'S1', dimensions=('ngases', 'nchars'))  # Names of gases
        var_mr =        ds.createVariable('x_gas',   'f4', dimensions=('nlev_c', 'ngases'))  # Mixing ratios per level
        var_cmr =       ds.createVariable('x_cond',  'f4', dimensions=('nlev_c', 'ngases'))  # Condensate mixing ratios per level
        var_pvol =      ds.createVariable('p_vol',   'f4', dimensions=('nlev_c', 'ngases')); var_pvol.units = "Pa"  # Gas phase partial pressures
        var_plvol =     ds.createVariable('pl_vol',  'f4', dimensions=('nlev_l', 'ngases')); var_pvol.units = "Pa"
        var_fdl =       ds.createVariable('fl_D_LW', 'f4', dimensions=('nlev_l'));           var_fdl.units = "W m-2"
        var_ful =       ds.createVariable('fl_U_LW', 'f4', dimensions=('nlev_l'));           var_ful.units = "W m-2"
        var_fnl =       ds.createVariable('fl_N_LW', 'f4', dimensions=('nlev_l'));           var_fnl.units = "W m-2"
        var_fds =       ds.createVariable('fl_D_SW', 'f4', dimensions=('nlev_l'));           var_fds.units = "W m-2"
        var_fus =       ds.createVariable('fl_U_SW', 'f4', dimensions=('nlev_l'));           var_fus.units = "W m-2"
        var_fns =       ds.createVariable('fl_N_SW', 'f4', dimensions=('nlev_l'));           var_fns.units = "W m-2"
        var_fd =        ds.createVariable('fl_D',    'f4', dimensions=('nlev_l'));           var_fd.units = "W m-2"
        var_fu =        ds.createVariable('fl_U',    'f4', dimensions=('nlev_l'));           var_fu.units = "W m-2"
        var_fn =        ds.createVariable('fl_N',    'f4', dimensions=('nlev_l'));           var_fn.units = "W m-2"
        var_hr =        ds.createVariable('rad_hr',  'f4', dimensions=('nlev_c'));           var_hr.units = "K day-1"
        var_sful =      ds.createVariable('Sfl_U_LW','f4', dimensions=('nbands', 'nlev_l')); var_sful.units = "W m-2 m-1"
        var_sfus =      ds.createVariable('Sfl_U_SW','f4', dimensions=('nbands', 'nlev_l')); var_sfus.units = "W m-2 m-1"

        var_re =        ds.createVariable('re',      'f4', dimensions=('nlev_c'));           var_re.units = "m"
        var_lwm =       ds.createVariable('lwm',     'f4', dimensions=('nlev_c'));           var_lwm.units = "kg kg-1"
        var_clfr =      ds.createVariable('clfr',    'f4', dimensions=('nlev_c'));           var_clfr.units = "none"

        #     Store data
        var_p[:] =      self.p[:]
        var_pl[:] =     self.pl[:]
        var_tmp[:] =    self.tmp[:]
        var_tmpl[:] =   self.tmpl[:]
        var_z[:]    =   self.z[:]
        var_zl[:]    =  self.zl[:]
        var_mmw[:]  =   self.mu[:]
        var_grav[:]  =  self.grav_z[:]
        var_cp[:]  =  self.cp[:]
        var_xd[:] = self.xd[:]
        var_xv[:] = self.xv[:]

        var_gases[:] =  np.array([ [c for c in gas.ljust(nchars)[:nchars]] for gas in gas_list ] , dtype='S1')
        var_mr[:] =     np.array([ [ self.x_gas[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_cmr[:] =    np.array([ [ self.x_cond[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_pvol[:] =   np.array([ [ self.p_vol[gas][i] for i in range(nlev_c-1,-1,-1) ] for gas in gas_list  ]).T
        var_plvol[:] = np.array([ [ self.pl_vol[gas][i] for i in range(nlev_l-1,-1,-1) ] for gas in gas_list  ]).T
        var_fdl[:] =    self.LW_flux_down[:]
        var_ful[:] =    self.LW_flux_up[:]
        var_fnl[:] =    self.LW_flux_net[:]

        var_fds[:] =    self.SW_flux_down[:]
        var_fus[:] =    self.SW_flux_up[:]
        var_fns[:] =    self.SW_flux_net[:]

        var_fd[:] =     self.flux_down_total[:]
        var_fu[:] =     self.flux_up_total[:]
        var_fn[:] =     self.net_flux[:]

        var_hr[:] =     self.net_heating[:]

        var_sful[:,:] = self.LW_spectral_flux_up[:,:]
        var_sfus[:,:] = self.SW_spectral_flux_up[:,:]

        var_re[:]   =   self.re[:]
        var_lwm[:]  =   self.lwm[:]
        var_clfr[:] =   self.clfr[:]

        # ----------------------
        # Close
        ds.close()


