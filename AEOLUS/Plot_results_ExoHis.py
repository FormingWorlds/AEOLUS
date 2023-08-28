#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 22:32:14 2023

@author: x_anfey
"""

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import pickle as pkl
import os,sys


##################################################################################################
#################################### SECTION 0 : SETTINGS ########################################
##################################################################################################

from ExoHiS.Gen_files import settings

settings = settings()

path_ExoHiS = '/proj/bolinc/users/x_anfey/AEOLUS/ExoHiS/' 
path_output = '/proj/bolinc/users/x_anfey/AEOLUS/output/'



##################################################################################################
####################### SECTION 1 : k-coefficients distribution in CPS ###########################
##################################################################################################

def plot_kcoeff_CPS(T,P,band):
    
    t       = settings.Tg.index(T)
    p       = settings.Pg.index(P)
    pt_pair = len(settings.Tg)*p + t
    
    ds      = netCDF4.Dataset(path_ExoHiS+'block5/output_ExoMol_'+settings.resxf+'.nc')
    k_coeff = ds['kopt'][band,pt_pair,:20]
    weights = ds['w_k'][band,:20]
    
    k_coeff_sorted = []
    g_values       = [0]
    
    for k in range(len(k_coeff)):
        k_min = min(k_coeff)
        c     = 0
        while k_coeff[c] != k_min :
            c+=1
        k_coeff_sorted.append(k_min)
        k_coeff[c]=10**10
        g_values.append(g_values[-1]+weights[c])
    del(g_values[0])
    while g_values[0] == 0 and k_coeff_sorted[0] == 0 :
        del(g_values[0])
        del(k_coeff_sorted[0])
    
    plt.plot(g_values,k_coeff_sorted,label='ExoMol (corr_k)',color="red",linestyle="",marker="o",ms=5)
    
    
    k_coeff  = settings.tmp.kdata[p,t,band,:]
    g_values = settings.tmp.ggrid
    plt.plot(g_values,k_coeff,label='ExoMol (exo_k)',color="orange",linestyle="",marker="o",ms=5)
    
    
    dsH      = netCDF4.Dataset(path_ExoHiS+'block5/output_HITRAN_'+settings.resxf+'.nc')
    k_coeffH = dsH['kopt'][band,pt_pair,:20]
    weightsH = dsH['w_k'][band,:20]
    k_coeff_sortedH=[]
    g_valuesH=[0]
    
    for k in range(len(k_coeffH)):
        k_min = min(k_coeffH)
        c     = 0
        while k_coeffH[c] != k_min :
            c+=1
        k_coeff_sortedH.append(k_min)
        k_coeffH[c]=10**10
        g_valuesH.append(g_valuesH[-1]+weightsH[c])
    del(g_valuesH[0])
    while g_valuesH[0] == 0 and k_coeff_sortedH[0] == 0 :
        del(g_valuesH[0])
        del(k_coeff_sortedH[0])
    
    plt.plot(g_valuesH,k_coeff_sortedH,label='HITRAN',color="black",linestyle="",marker="o",ms=5)
    
    plt.yscale('log')
    plt.grid()
    plt.legend(fontsize=13)
    plt.title(r'$(P,T)=($'+str(P)+'$~bar,$'+str(T)+'$~K)$, band : '+str(band)+r', $\widetilde{\nu} = $'+str(settings.wnedges[band])+' $cm^{-1}$')
    plt.xlabel('Cumulative probability $g$',fontsize=13)
    plt.ylabel('Absorption coefficient $\kappa~[m^2.kg^{-1}]$',fontsize=13)
    plt.show()
    
##################################################################################################
############################## SECTION 2 : transmission function #################################
##################################################################################################


def plot_trans_func(T,P):
    
    t       = settings.Tg.index(T)
    p       = settings.Pg.index(P)
    pt_pair = len(settings.Tg)*p + t
    
    mean_TEC,mean_TEE,mean_TH = [],[],[]
    inv = [-1 for k in range(20)]
    for band in range(len(settings.wnedges)-1):
        
        dsEC        = netCDF4.Dataset(path_ExoHiS+'block5/output_ExoMol_'+settings.resxf+'.nc')
        k_coeffEC   = dsEC['kopt'][band,pt_pair,:20]
        weightsEC   = dsEC['w_k'][band,:20]
        mean_TEC.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEC)),weightsEC)))
        
         
        dsH = netCDF4.Dataset(path_ExoHiS+'block5/output_HITRAN_'+settings.resxf+'.nc')
        k_coeffH = dsH['kopt'][band,pt_pair,:20]
        weightsH = dsH['w_k'][band,:20]
        mean_TH.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffH)),weightsH)))
        
        k_coeffEE = settings.tmp.kdata[p,t,band,:]
        gvaluesEE = settings.tmp.ggrid
        weightsEE = [gvaluesEE[0]]
        for i in range(1,len(gvaluesEE)):
            weightsEE.append(gvaluesEE[i]-gvaluesEE[i-1])
        mean_TEE.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEE)),weightsEE)))
    
    
    diffEC_H  = [mean_TEC[k]-mean_TH[k] for k in range(len(mean_TEC))]
    diffEE_H  = [mean_TEE[k]-mean_TH[k] for k in range(len(mean_TEE))]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
     
    ax1.plot(settings.wnedges[:-1],mean_TEC,label='ExoMol (corr_k)',color='red')
    ax1.plot(settings.wnedges[:-1],mean_TEE,label='ExoMol (exo_k)',color='orange')
    ax1.plot(settings.wnedges[:-1],mean_TH,label='HITRAN',color='black')
    ax1.set_xlim(-500,settings.wnedges[-1])
    
    ax2.plot(settings.wnedges[:-1],diffEC_H,label='HITRAN vs ExoMol (corr_k)',color='red')
    ax2.plot(settings.wnedges[:-1],diffEE_H,label='HITRAN vs ExoMol (exo_k)',color='orange')
    ax2.set_xlim(-500,settings.wnedges[-1])

    plt.grid(visible=True)
    ax1.legend() 
    ax2.legend()
    ax1.set_ylim(-2e-2,1.02)
        
    ax2.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    ax2.set_ylabel('difference',fontsize=13)
     
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    ax1.set_ylabel('transmission coefficient $T$',fontsize=13)
     
    ax3 = ax1.twiny()
    ax3.set_xlabel('wavelength '+r'$\lambda~[\mu m]$',fontsize=13)
    
    ax1.set_xticks([1]+list(np.arange(5000,35000,5000)))    
    ax2.set_xticks([1]+list(np.arange(5000,35000,5000)))
    wl_ticks=[round(10**4/wn,9) for wn in ax1.get_xticks()]
    wl_ticks_ploted = [round(wl,3) for wl in wl_ticks]
     
    ax3.set_xticks([10**4/x for x in wl_ticks])
    ax3.set_xbound(ax1.get_xbound())
    ax3.set_xticklabels(wl_ticks_ploted,)
     
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax3.tick_params(axis='both', which='major', labelsize=12)
     
    plt.title(r'$(P,T)=($'+str(P)+'$~bar,$'+str(T)+'$~K)$')
    plt.show()

    
def plot_comp_grid():
    dsEC = netCDF4.Dataset(path_ExoHiS+'block5/output_ExoMol_'+settings.resxf+'.nc')
    dsH  = netCDF4.Dataset(path_ExoHiS+'block5/output_HITRAN_'+settings.resxf+'.nc')
    lenP, lenT = len(settings.Pg), len(settings.Tg)
    mean_diffEC, mean_diffEE = np.zeros((lenP,lenT)),np.zeros((lenP,lenT))
    for pt_pair in range(lenP*lenT):
        t_index, p_index = pt_pair%lenT, pt_pair//lenT
        TEC,TEE,TH = [],[],[]
        inv = [-1 for k in range(20)]
        for band in range(len(settings.wnedges)-1):
             
            k_coeffEC = dsEC['kopt'][band,pt_pair,:20]
            weightsEC = dsEC['w_k'][band,:20]
            TEC.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEC)),weightsEC)))
          
         
            k_coeffH = dsH['kopt'][band,pt_pair,:20]
            weightsH = dsH['w_k'][band,:20]
            TH.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffH)),weightsH)))
         
            k_coeffEE = settings.tmp.kdata[p_index,t_index,band,:]
            gvaluesEE = settings.tmp.ggrid
            weightsEE = [gvaluesEE[0]]
            for i in range(1,len(gvaluesEE)):
                weightsEE.append(gvaluesEE[i]-gvaluesEE[i-1])
            TEE.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEE)),weightsEE)))
            
        diffEC_H  = np.abs([TEC[k]-TH[k] for k in range(len(TEC))])
        diffEE_H  = np.abs([TEE[k]-TH[k] for k in range(len(TEE))])
         
        mean_diffEC[p_index,t_index] = np.sum([diffEC_H[k]*(settings.wnedges[k+1]-settings.wnedges[k]) for k in range(len(TEC))])*100/settings.wnedges[-1]
        mean_diffEE[p_index,t_index] = np.sum([diffEE_H[k]*(settings.wnedges[k+1]-settings.wnedges[k]) for k in range(len(TEE))])*100/settings.wnedges[-1]
      
    p   = np.array(settings.Pg)
    t   = np.array(settings.Tg)
    
    
    
    with open(path_output+'ExoMol_'+settings.resxf+'/300.0K_1.0bar_moist_pt.tsv', 'r') as f:
        lignes1 = f.readlines()
    lignes_modifiees1 = lignes1[2:]
    with open(path_output+'ExoMol_'+settings.resxf+'/300.0K_1.0bar_moist_pt_rpz.tsv', 'w') as f:
        f.writelines(lignes_modifiees1)
    
    p_and_t1 = np.fromfile(path_output+'ExoMol_'+settings.resxf+'/300.0K_1.0bar_moist_pt_rpz.tsv',count=-1, sep=' ')
    pprofile1, tprofile1 = [], []
    for i in range(len(p_and_t1)//2):
        pprofile1.append(1e-5*p_and_t1[2*i])
        tprofile1.append(p_and_t1[2*i+1])
    
    
    
    
    with open(path_output+'ExoMol_'+settings.resxf+'/2500.0K_200.0bar_moist_pt.tsv', 'r') as f:
        lignes1 = f.readlines()
    lignes_modifiees1 = lignes1[2:]
    with open(path_output+'ExoMol_'+settings.resxf+'/2500.0K_200.0bar_moist_pt_rpz.tsv', 'w') as f:
        f.writelines(lignes_modifiees1)
    
    p_and_t2 = np.fromfile(path_output+'ExoMol_'+settings.resxf+'/2500.0K_200.0bar_moist_pt_rpz.tsv',count=-1, sep=' ')
    pprofile2, tprofile2 = [], []
    for i in range(len(p_and_t2)//2):
        pprofile2.append(1e-5*p_and_t2[2*i])
        tprofile2.append(p_and_t2[2*i+1])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    vmin = min(np.array(mean_diffEC).min(),np.array(mean_diffEE).min())
    vmax = max(np.array(mean_diffEC).max(),np.array(mean_diffEE).max())
    
    p1 = ax1.pcolormesh(t,p,mean_diffEC,vmin=vmin,vmax=vmax)
    p2 = ax2.pcolormesh(t,p,mean_diffEE,vmin=vmin,vmax=vmax)
        
    fig.colorbar(p1, ax=ax1, orientation='horizontal', label='transmission absolute difference [%]')
    fig.colorbar(p2, ax=ax2, orientation='horizontal', label='transmission absolute difference [%]')
    
    ax1.plot(tprofile1, pprofile1, color='red', label=r'$(T_{surf},P_{surf})=(300~K,1~bar)$')
    ax2.plot(tprofile1, pprofile1, color='red', label=r'$(T_{surf},P_{surf})=(300~K,1~bar)$')
    
    ax1.plot(tprofile2, pprofile2, color='red', linestyle='--', label=r'$(T_{surf},P_{surf})=(2500~K,200~bar)$')
    ax2.plot(tprofile2, pprofile2, color='red', linestyle='--', label=r'$(T_{surf},P_{surf})=(2500~K,200~bar)$')

    ax1.set_yscale('log')
    ax1.set_ylim(5e-5,2e2)
    ax1.set_ylabel('Pressure $P~[bar]$',fontsize=13)
    ax1.set_xlabel('Temperature $T~[K]$',fontsize=13)
    ax1.legend(loc='lower center')
    ax1.set_title('HITRAN vs ExoMol (corr_k)')
    
    ax2.set_yscale('log')
    ax2.set_ylim(5e-5,2e2)
    ax2.set_ylabel('Pressure $P~[bar]$',fontsize=13)
    ax2.set_xlabel('Temperature $T~[K]$',fontsize=13)
    ax2.legend(loc='lower center')
    ax2.set_title('HITRAN vs ExoMol (exo_k)')
    
    fig.tight_layout()
    plt.show()
    

##################################################################################################
############################## SECTION 3 : results from SOCRATES #################################
##################################################################################################


def plot_OLR(T,P):       
    with open(path_output+'HITRAN_'+settings.resxf+'/'+str(T)+'.0K_'+str(P)+'bar_atm.pkl','rb') as f: 
        content_HITRAN = pkl.load(f)  
    
    plt.plot(content_HITRAN.band_centres, (content_HITRAN.LW_spectral_flux_up[:,-1]+content_HITRAN.SW_spectral_flux_up[:,-1])/content_HITRAN.band_widths, color="purple",linestyle='--',label='Blackbody emission $(300~K)$')
    plt.plot(content_HITRAN.band_centres, (content_HITRAN.LW_spectral_flux_up[:,0]+content_HITRAN.SW_spectral_flux_up[:,0])/content_HITRAN.band_widths, color="black",label='HITRAN TOA')
    
    
    with open(path_output+'ExoMol_'+settings.resxf+'/'+str(T)+'.0K_'+str(P)+'bar_atm.pkl','rb') as f: 
        content_ExoMol = pkl.load(f)  
    
    plt.plot(content_ExoMol.band_centres, (content_ExoMol.LW_spectral_flux_up[:,0]+content_ExoMol.SW_spectral_flux_up[:,0])/content_ExoMol.band_widths, color="red",label='ExoMol TOA')
       
        
    plt.ylabel(r'Spectral flux density $F_{\widetilde{\nu}}$ [W m$^{-2}$ [cm$^{-1}]^{-1}$]',fontsize=13)
    plt.xlabel(r'Wavenumber $\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    plt.legend(loc=1)
    plt.title(r'$(P,T)=('+str(P)+'~bar,$'+str(T)+'$~K)$')
    plt.xlim(left=0, right=7000)
            
def plot_heating(T,P):   
       
    with open(path_output+'HITRAN_'+settings.resxf+'/'+str(T)+'.0K_'+str(P)+'bar_atm.pkl','rb') as f: 
        content_HITRAN = pkl.load(f)  
    
    plt.plot(content_HITRAN.net_heating,content_HITRAN.p, color="black",label='HITRAN')
    
    with open(path_output+'ExoMol_'+settings.resxf+'/'+str(T)+'.0K_'+str(P)+'bar_atm.pkl','rb') as f: 
        content_ExoMol = pkl.load(f)  
    
    plt.plot(content_ExoMol.net_heating,content_ExoMol.p, color="red",label='ExoMol')
    
    
    plt.gca().invert_yaxis()
    plt.legend()
    plt.yscale('log')
    plt.axvline(x=0,color='grey',linestyle='--')
    plt.xlim(-200,200)
    plt.ylabel(r'Pressure $P$ [Pa]',fontsize=13)
    plt.xlabel(r'Heating [K/day]',fontsize=13)
    plt.title(r'$(P,T)=(10~bar,$'+str(T)+'$~K)$')
    plt.show()