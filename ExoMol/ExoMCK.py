#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 20:12:27 2023

@author: x_anfey
"""

import exo_k as xk
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
# import seaborn as sns
import copy
import shutil
import pathlib
import pickle as pkl
import json
import glob, re, os
import netCDF4



"""
################## Downloading all the necessary files from ExoMol database ######################


os.system("wget https://www.exomol.com/db//H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL.states.bz2")
os.system("wget https://www.exomol.com/db//H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL.pf")
os.system("wget https://www.exomol.com/db/H2O/1H2-16O/1H2-16O__air_a0.broad")

for k in range (301):
    lwn=k*100
    uwn=(k+1)*100
    if uwn<1000:
        address="https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__00" +str(lwn)+ "-00"+str(uwn)+".trans.bz2"
    elif uwn==1000:
        address="https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__00900-01000.trans.bz2"
    elif 1000<uwn<10000:
        address="https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__0" +str(lwn)+ "-0"+str(uwn)+".trans.bz2"
    elif uwn==10000:
        address="https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__09900-10000.trans.bz2"
    else :
        address="https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__" +str(lwn)+ "-"+str(uwn)+".trans.bz2"
    os.system("wget "+address)

os.system("bunzip2 *.bz2")

######################## Converting xsec files (ASCII) to hdf5 format #############################
"""
path_in='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_10_files/'
path_out='/proj/bolinc/users/x_anfey/ExoMol/hdf5_cstep_10_files/'

T_grid=[100,300,500,700,1000,1300,1500,2000,2500]
P_grid=[1e-4,1e-3,1e-2,1e-1,5e-1,1e0,5e0,1e1,2.5e1,5e1,7.5e1,1e2,1.25e2,1.5e2,1.75e2,2e2]

T_grid_bis=[300,500,1000,1500,2000,2500]
P_grid_bis=[1e-4,1e-3,1e-2,1e-1,1e0,1e1,5e1,1e2,1.5e2,2e2]

# T & P grids names of files

T_grid_str=[]
P_grid_str=[]

for T in T_grid_bis:
    if T<1000:
        T_grid_str.append('0'+str(T))
    else :
        T_grid_str.append(str(T))

for P in P_grid_bis:
    if P<1:
        P_grid_str.append(str(P))
    elif 1<=P<10:
        P_grid_str.append('00'+str(int(P)))
    elif 10<=P<100:
        P_grid_str.append('0'+str(int(P)))
    else:
        P_grid_str.append(str(int(P)))


# conversion from cm²/molecule to m²/kg : m²/kg = mult_factor * cm²/molecule

Na = 6.0221408e+23 # Avogadro number
M = 18.01528 # Molar masse of molecule considered [g/mol] (water)

mult_factor = Na/(10*M)
"""
for T_str in T_grid_str:
    for P_str in P_grid_str:
        file='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_10hires'
        
        xk.convert_kspectrum_to_hdf5(path_in+file+'.xsec', path_out+file+'.h5',wn_column=0, kdata_column=1, file_kdata_unit='cm^2', data_type='xsec',mult_factor=mult_factor)

"""
########################### Creating a corr-k matrix with the spectra ####################################

path_out_01 = '/proj/bolinc/users/x_anfey/ExoMol/hdf5_cstep_files/'
path_out_02 = '/proj/bolinc/users/x_anfey/ExoMol/hdf5_cstep_02_files/'
path_out_10 = '/proj/bolinc/users/x_anfey/ExoMol/hdf5_cstep_10_files/'

file_grid_01=[]
for P_str in P_grid_str:
    P_line=[]
    for T_str in T_grid_str:
        P_line.append('1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_hires'+'.h5')
    file_grid_01.append(P_line)
file_grid_01=np.array(file_grid_01)

file_grid_02=[]
for P_str in P_grid_str:
    P_line=[]
    for T_str in T_grid_str:
        P_line.append('1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_02hires'+'.h5')
    file_grid_02.append(P_line)
file_grid_02=np.array(file_grid_02)

file_grid_10=[]
for P_str in P_grid_str:
    P_line=[]
    for T_str in T_grid_str:
        P_line.append('1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_10hires'+'.h5')
    file_grid_10.append(P_line)
file_grid_10=np.array(file_grid_10)

Tgrid=np.array(T_grid_bis)
logPgrid=np.array(np.log10(P_grid_bis))
wnedges=np.array(np.concatenate((np.arange(0,3000,25),np.arange(3000,11000,50),np.arange(11000,30500,500))))


#tmp_01=xk.hires_to_ktable(path=path_out_01,filename_grid=file_grid_01,logpgrid=logPgrid, grid_p_unit='bar', p_unit='Pa', tgrid=Tgrid, wnedges=wnedges, mol='H2O')
#tmp.remove_zeros()
#print(tmp)
#print(tmp.kdata[5,1,280,:])
#tmp_02=xk.hires_to_ktable(path=path_out_02,filename_grid=file_grid_02,logpgrid=logPgrid, grid_p_unit='bar', p_unit='Pa', tgrid=Tgrid, wnedges=wnedges, mol='H2O')
#tmp_10=xk.hires_to_ktable(path=path_out_10,filename_grid=file_grid_10,logpgrid=logPgrid, grid_p_unit='bar', p_unit='Pa', tgrid=Tgrid, wnedges=wnedges, mol='H2O')


###################################### Ploting cross-sections ###########################################

def plot_cs_T(Pressure,T_grid):
    path='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_10_files/'
    
    if Pressure<1:
        P_str=str(Pressure)
    elif 1<=Pressure<10:
        P_str='00'+str(Pressure)
    elif 10<=Pressure<1000:
        P_str='0'+str(Pressure)
    else:
        P_str=str(Pressure)
    
    fig = plt.figure()
    ax1= fig.add_subplot(111)
    
    for T in T_grid :
        if 100<=T<1000:
            T_str='0'+str(T)
        else :
            T_str=str(T)
        
        filename='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_10hires.xsec'
        wn_and_xs=np.fromfile(path+filename,count=-1, sep=' ')
        wn,xs=[],[]
        for i in range(len(wn_and_xs)//2):
            wn.append(wn_and_xs[2*i])
            xs.append(mult_factor*wn_and_xs[2*i+1])
        ax1.plot(wn,xs,label=r'$(P,T)=($'+str(Pressure)+' $bar,$'+str(T)+' $K)$',color='black')
    
    plt.yscale('log')
    plt.grid()
    plt.legend()
    #ax1.set_xlim(11098.2,11099.8)
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$')
    ax1.set_ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('wavelength '+r'$\lambda~[\mu m]$')
    
    print(ax1.get_xticks())
    wl_ticks=[round(10**4/wn,2) for wn in ax1.get_xticks()]
    
    ax2.set_xticks([10**4/x for x in wl_ticks])
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(wl_ticks)
    
    plt.title(r'$P~=~$'+str(Pressure)+' atm')
    plt.show()

def plot_cs_T_bis():
    path='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_10_files/'
    
    fig = plt.figure()
    ax1= fig.add_subplot(111)
    
        
    filename='1H2-16O_'+'0300'+'K_Voigt_P'+'001'+'_cstep_10hires.xsec'
    wn_and_xs=np.fromfile(path+filename,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label='1'+' bar',color='black')
    
    filename='1H2-16O_'+'0300'+'K_Voigt_P'+'010'+'_cstep_10hires.xsec'
    wn_and_xs=np.fromfile(path+filename,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label='10'+' bar',color='yellow')
    
    filename='1H2-16O_'+'0300'+'K_Voigt_P'+'100'+'_cstep_10hires.xsec'
    wn_and_xs=np.fromfile(path+filename,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label='100'+' bar',color='red')
    
    plt.yscale('log')
    plt.grid()
    plt.legend()
    #ax1.set_xlim(400,1600)
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$')
    ax1.set_ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('wavelength '+r'$\lambda~[\mu m]$')
    
    print(ax1.get_xticks())
    wl_ticks=[round(10**4/wn,2) for wn in ax1.get_xticks()]
    
    ax2.set_xticks([10**4/x for x in wl_ticks])
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(wl_ticks)
    
    plt.title(r'$P~=~$'+'1'+' atm')
    plt.show()

def plot_cs_P(P_grid,Temperature):
    path='/proj/bolinc/users/x_anfey/ExoMol/xsec_files/'
      
    if 100<=Temperature<1000:
        T_str='0'+str(Temperature)
    else :
        T_str=str(Temperature)
      
    for P in P_grid :
                    
        if P<1:
            P_str=str(P)
        elif 1<=P<10:
            P_str='00'+str(P)
        elif 10<=P<100:
            P_str='0'+str(P)
        else:
            P_str=str(P)
          
        filename='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'.xsec'
        wn_and_xs=np.fromfile(path+filename,count=-1, sep=' ')
        wn,xs=[],[]
        for i in range(len(wn_and_xs)//2):
            wn.append(wn_and_xs[2*i])
            xs.append(mult_factor*wn_and_xs[2*i+1])
        plt.plot(wn,xs,label=str(P)+' atm')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$')
    plt.ylabel('absorption cross section $\kappa~[m^2.kg^{-1}]$')
    plt.title(r'$T~=~$'+str(Temperature)+' K')
    plt.show()
    
def plot_cs_res(Pressure,Temperature): 
    if Pressure<1:
        P_str=str(Pressure)
    elif 1<=Pressure<10:
        P_str='00'+str(Pressure)
    elif 10<=Pressure<1000:
        P_str='0'+str(Pressure)
    else:
        P_str=str(Pressure)
    
    if 100<=Temperature<1000:
        T_str='0'+str(Temperature)
    else :
        T_str=str(Temperature)
    
    fig = plt.figure()
    ax1= fig.add_subplot(111)
    
    filename_01='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_hires.xsec'
    path_01='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_files/'
    wn_and_xs=np.fromfile(path_01+filename_01,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label=r'$\Delta \widetilde{\nu}=1.0~cm^{-1}$',color='yellow')
    
    
    
    filename_02='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_02hires.xsec'
    path_02='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_02_files/'
    wn_and_xs=np.fromfile(path_02+filename_02,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label=r'$\Delta \widetilde{\nu}=0.5~cm^{-1}$',color='orange')
    
    filename_10='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_10hires.xsec'
    path_10='/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_10_files/'
    wn_and_xs=np.fromfile(path_10+filename_10,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label=r'$\Delta \widetilde{\nu}=0.1~cm^{-1}$',color='red')
    
    filename_100='1H2-16O_'+T_str+'K_Voigt_P'+P_str+'_cstep_100hires.xsec'
    path_100='/proj/bolinc/users/x_anfey/ExoMol/'
    wn_and_xs=np.fromfile(path_100+filename_100,count=-1, sep=' ')
    wn,xs=[],[]
    for i in range(len(wn_and_xs)//2):
        wn.append(wn_and_xs[2*i])
        xs.append(mult_factor*wn_and_xs[2*i+1])
    ax1.plot(wn,xs,label=r'$\Delta \widetilde{\nu}=0.01~cm^{-1}$',color='black')
    
    plt.yscale('log')
    plt.grid(visible=True)
    plt.legend()
    ax1.set_xlim(11098.2,11099.8)
    #ax1.set_xlim(11095,11100)
    
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=14)
    ax1.set_ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$',fontsize=14)
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('wavelength '+r'$\lambda~[\mu m]$',fontsize=14)
    
    
    #ax1.set_xlabel('''nombre d'onde '''+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=14)
    #ax1.set_ylabel('''coefficient d'absorption $\kappa~[m^2.kg^{-1}]$''',fontsize=14)
    
    #ax2 = ax1.twiny()
    #ax2.set_xlabel('''longueur d'onde '''+r'$\lambda~[\mu m]$',fontsize=14)
    
    
    print(ax1.get_xticks())
    wl_ticks=[round(10**4/wn,4) for wn in ax1.get_xticks()]
    
    ax2.set_xticks([10**4/x for x in wl_ticks])
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(wl_ticks,)
    
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    
    #plt.title(r'$(P,T)~=~($'+str(Pressure)+' atm,'+str(Temperature)+' K$)$')
    plt.show()

"""
####################### Writing lookup table for SOCRATES ##############################################

lookup_name = 'sp_b318_ExoMol_k'
with open(lookup_name,"a") as file:
    file.write('*BLOCK: k-table\n\n')
    file.write('Lookup table:\t'+str(len(P_grid))+' pressures,\t'+str(len(T_grid))+' temperatures.\n')
    for p in range(len(P_grid)):
        file.write(' '+str(format(tmp.pgrid[p],'.6E'))+' ')
        c=1 
        for t in range(len(T_grid)):
            if c%6==0:
                file.write('\n '+str(format(tmp.tgrid[t],'.6E'))+' ')
                c+=1
            else :
                file.write(str(format(tmp.tgrid[t],'.6E'))+' ')
                c+=1
        file.write('\n')
    file.write('\n')
    for b in range(len(tmp.wns)):
        file.write('Band:\t'+str(b+1)+', gas:\t'+str(1)+', k-terms:\t'+str(tmp.Ng)+'\n')
        for w in range(tmp.Ng):
            for p in range(len(P_grid)):
                file.write(' '+str(format(tmp.kdata[p,0,b,w],'.6E'))+' ')
                c=1
                for t in range(1,len(T_grid)):
                    if c%6==0:
                        file.write('\n '+str(format(tmp.kdata[p,t,b,w],'.6E'))+' ')
                        c+=1
                    else :
                        file.write(str(format(tmp.kdata[p,t,b,w],'.6E'))+' ')
                        c+=1
                file.write('\n')
        file.write('\n')
    file.write('\n*END')

############################# Writing spectral files with prep_spec ######################################

#os.system('rm sp_b318_ExoMol')
exec_file_name = "sp_exec_b318_ExoMol.sh"

f = open(exec_file_name, "w+")

f.write('prep_spec <<EOF'+ '\n')
f.write('sp_b318_ExoMol_CIA_bis\n')
f.write('318\n')                                  # number of bands
f.write('1\n')                                    # number of absorbing gases
f.write('1\n')                                    # absorber ids (1 for water)
f.write('1\n')                                    # number of generalised continua to be included
f.write('1 1\n')                                  # gas identifiers for continuum 1
f.write('0\n')                                    # number of aerosols
f.write('c\n')                                    # band units (c for inverse cm)

f.write('1'+' ')                                    # write band edges one by one
for band in wnedges[1:-1]:                        #
    f.write(str(band)+'\n')                       #
    f.write(str(band)+' ')                        #
f.write(str(wnedges[-1])+'\n')                    #


for band in wnedges[:-1]:                         # absorber ids in each band 
    f.write('1 \n')                               #

for band in wnedges[:-1]:                         # continua ids in each band
    f.write('0 \n')                               #
    
f.write('n\n')                                    # exclude no bands

f.write('-1\n')                                   # close prep_spec
f.write('EOF\n')                                  #

f.close()
os.chmod(exec_file_name,0o777)

os.system('./'+exec_file_name)                    # run file script

########################################## Writing block 5 ###########################################

block5 = "spectral_file_block5"
with open(block5,"a") as file:
    file.write('*BLOCK: TYPE ='+'\t'+'5: SUBTYPE ='+'\t'+'0: version ='+'\t'+'1'+'\n')
    file.write('Exponential sum fiting coefficients: (exponents: m2/kg)'+'\n')
    file.write('Band'+'\t'+'Gas, Number of k-terms, Scaling type and scaling function,'+'\n')
    file.write('\t\t'+'followed by reference pressure and temperature,'+'\n')
    file.write('\t\t\t'+'k-terms, weights and scaling parameters.'+'\n')
    for b in range(len(tmp.wns)):
        values=str(b+1)+'\t'+'1'+'\t'+'20'+'\t'+'2'+'\t'+'9'+'\t'+str(len(P_grid))+'\t'+str(len(T_grid))+'\n'
        for iw in range(len(tmp.weights)):
            values+='    '+str(format(tmp.kdata[-1,-1,b,iw],'.9E'))+'    '+str(format(tmp.weights[iw],'.9E'))+'\n'
        file.write(values)
    file.write('*END')

####################################### Creating a netCDF file ########################################

os.system("rm 1H2-16O_ExoMol.nc")

# create a netCDF file
ncfile = netCDF4.Dataset("1H2-16O_ExoMol_cstep.nc","w",format='NETCDF3_CLASSIC',clobber='true')

# create the dimensions for the netCDF file
##ncfile.createDimension('scalar',1)
ncfile.createDimension('nu',30501)
ncfile.createDimension('pt_pair',len(tmp.pgrid)*len(tmp.tgrid))

# create the variables
p_calc = ncfile.createVariable('p_calc','f4',('pt_pair',))
t_calc = ncfile.createVariable('t_calc','f4',('pt_pair',))
nu = ncfile.createVariable('nu', 'f8', ('nu',))
kabs = ncfile.createVariable('kabs','f8',('pt_pair','nu',))

# assign the data
with open(path_in+'1H2-16O_'+T_grid_str[0]+'K_Voigt_P'+P_grid_str[0]+'.xsec',"r") as f:
    WN = list()
    for line in f:
        data = line.split()
        WN.append(data[0])
nu[:] = [float(wn) for wn in WN]

t_calc[:] = (tmp.tgrid.tolist())*len(tmp.pgrid)
p_calc[:] = sorted(tmp.pgrid.tolist()*len(tmp.tgrid))


for i in range(len(tmp.pgrid)*len(tmp.tgrid)):
    n = len(tmp.tgrid)
    p = i//n
    t = i%n
    with open(path_in+'1H2-16O_'+T_grid_str[t]+'K_Voigt_P'+P_grid_str[p]+'.xsec',"r") as f:
        ACS = list()
        for line in f:
            data = line.split()
            ACS.append(data[1])
        kabs[i,:] = [float(acs) for acs in ACS]
        
ncfile.close()


####################################### Creating a netCDF file BIS ########################################

os.system("rm 1H2-16O_ExoMol.nc")

# create a netCDF file
ncfile = netCDF4.Dataset("1H2-16O_ExoMol.nc","w",format='NETCDF3_CLASSIC',clobber='true')
##, format='NETCDF3_CLASSIC'

# create the dimensions for the netCDF file
ncfile.createDimension('scalar',1)
ncfile.createDimension('nu',305001)
ncfile.createDimension('pt_pair',len(tmp_10.pgrid)*len(tmp_10.tgrid))

# create the variables
nu = ncfile.createVariable('nu', 'f8', ('nu',))
kabs = ncfile.createVariable('kabs','f4',('pt_pair','nu',))
t_calc = ncfile.createVariable('t_calc','f8',('pt_pair',))
p_calc = ncfile.createVariable('p_calc','f8',('pt_pair',))



# assign the data
with open(path_in+'1H2-16O_'+T_grid_str[0]+'K_Voigt_P'+P_grid_str[0]+'_cstep_10hires'+'.xsec',"r") as f:
    WN = list()
    for line in f:
        data = line.split()
        WN.append(data[0])
nu[:] = [float(WN[i])*100 for i in range(len(WN))]

t_calc[:] = (T_grid_bis)*len(tmp_10.pgrid)
p_calc[:] = sorted([p*1e+5 for p in P_grid_bis]*len(tmp_10.tgrid))


for i in range(len(tmp_10.pgrid)*len(tmp_10.tgrid)):
    n = len(tmp_10.tgrid)
    p = i//n
    t = i%n
    with open(path_in+'1H2-16O_'+T_grid_str[t]+'K_Voigt_P'+P_grid_str[p]+'_cstep_10hires'+'.xsec',"r") as f:
        ACS = list()
        for line in f:
            data = line.split()
            ACS.append(data[1])
        kabs[i,:] = [float(ACS[j])*mult_factor for j in range(len(ACS))]
        
# assign attributes
#nu.title = 'wavenumber'
#nu.long_name = 'wavenumber'
#nu.units = 'm-1'
nu.step = nu[-1]-nu[-2]

#p_calc.title = 'pressure'
#p_calc.long_name = 'pressure'
#p_calc.units = 'Pa'

#t_calc.title = 'temperature'
#t_calc.long_name = 'temperature'
#t_calc.units = 'K'

#kabs.title = 'absorption'
#kabs.long_name = 'absorption'
#kabs.units = 'm2 kg-1'
        
ncfile.close()

######################################## Using corr_k #################################################

os.system("rm output_ExoMol")
os.system("rm monitoring_ExoMol")
os.system('HDF5_USE_FILE_LOCKING=FALSE')

exec_file_name = "corr_k_ExoMol.sh"

f = open(exec_file_name, "w+")

f.write('corr_k <<EOF'+ '\n')

#f.write('/proj/bolinc/users/x_anfey/socrates_2211/socrates_main/examples/corr_k/h2o_lbl.nc'+'\n')
f.write('1H2-16O_ExoMol.nc'+'\n')

f.write('N'+'\n')
f.write('N'+'\n')
f.write('sp_b318_ExoMol'+'\n')
f.write('N'+'\n')
f.write('Y'+'\n')
f.write('N'+'\n')
f.write('N'+'\n')
f.write('1'+'\n')
f.write('1 318'+'\n')
f.write('i'+'\n')

for P in P_grid_bis:
    f.write(str(P*1e+5))
    for T in T_grid_bis:
        f.write(' '+str(T))
    f.write('\n')
f.write('*END'+'\n')

f.write('Y'+'\n')
f.write('2'+'\n')
f.write('i'+'\n')

for k in range(318):
    f.write('1e+5 298'+'\n')
    
f.write('N'+'\n')
f.write('1'+'\n')
f.write('1.0e-2'+'\n')
f.write('10.0'+'\n')
f.write('N'+'\n')
f.write('N'+'\n')
f.write('Y'+'\n')
f.write('mapping_ExoMol'+'\n')
f.write('1'+'\n')
#f.write('/home/x_anfey/AEOLUS/profile.stoa'+'\n')
f.write('output_ExoMol'+'\n')
f.write('monitoring_ExoMol'+'\n')
f.write('8'+'\n')
f.write('EOF\n')

f.close()
os.chmod(exec_file_name,0o777)

os.system('./'+exec_file_name)

######################################### Trying Ccorr_k #################################################

#pt file
os.system('rm pt_file')
with open('pt_file',"a") as file:
    file.write('*PTVAL'+'\n')
    for P in P_grid_bis:
        file.write(str(P*1e+5))
        for T in T_grid_bis:
            file.write(' '+str(T))
        file.write('\n')
    file.write('*END')

#ref pt_file
os.system('rm ref_pt_file')
with open('ref_pt_file',"a") as file:
    file.write('*REF 1 318 1 1e+5 298')

#solar pathname
solar_path = '/proj/bolinc/users/x_anfey/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt'

#output pathname
output_path = 'output_ExoMol_k20_2306'
#os.system('rm '+output_path)
#os.system('rm '+output_path+'.nc')

#monitoring pathname
mon_path = 'monitoring_ExoMol_k20_2306'
#os.system('rm '+mon_path)

#LbL pathname
LbL_path = '1H2-16O_ExoMol.nc'

#self-broadening pathname
dir_continua_296 = '/proj/bolinc/users/x_anfey/spectral_files/dat_continua/mt_ckd_v3.0_s296'
dir_continua_260 = '/proj/bolinc/users/x_anfey/spectral_files/dat_continua/mt_ckd_v3.0_s260'

exec_file_name = "corr_k_ExoMol.sh"

f = open(exec_file_name, "w+")

f.write('Ccorr_k'+' ')
f.write('-s sp_b318_ExoMol'+' ')                 # spectral file
f.write('-R 1 318'+' ')                          # band limits 
f.write('-F pt_file'+' ')                        # pressures and temperatures at which to calculate coefficients
f.write('-r ref_pt_file'+' ')                    # reference conditions for scaling
f.write('-b 1.5e-3'+' ')                         # tolerance
f.write('-l 1 1.0e1'+' ')                        # generate line absorption data. gas id then maximum absorptive pathlength for the gas (kg/m2)
f.write('-n 20'+' ')                             # Number of k-terms the correlated-k fit should use
f.write('-lk'+' ')                               # a look-up table will be used for the pressure/temperature scaling
#f.write('-ct 1 1 10.0'+' ')
#f.write('-t 1.0e-3'+' ')
#f.write('-e '+dir_continua_296+' '+dir_continua_260+' ')
#f.write('+S '+solar_path+' ')                    # Solar Weighting
f.write('+p'+' ')                                # Planckian Weighting
f.write('-o '+output_path+' ')                   # Pathname of output file
f.write('-m '+mon_path+' ')                      # Pathname of monitoring file
f.write('-L '+LbL_path+' ')                      # Pathname of LbL file (.nc)
f.write('-np 1')                                 # Number of OpenMP threads

f.close()
os.chmod(exec_file_name,0o777)

os.system('./'+exec_file_name)


################################### Ploting CKD-method Exo_k ###############################################

def plot_kcoeff_CPS(nP,nT,nWN):
    k_coeff = tmp_10.kdata[nP,nT,nWN,:]
    g_values = tmp_10.ggrid
    plt.plot(g_values,k_coeff)
    plt.yscale('log')
    plt.grid()
    plt.xlabel('Cumulative probability $g$')
    plt.ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    plt.show()
    
def plot_kcoeff_CPS_bis(nP,nT,nWN):
    plt.plot(tmp_01.ggrid,tmp_01.kdata[nP,nT,nWN,:],label=r'$\widetilde{\nu}=1~cm^{-1}$')
    plt.plot(tmp_02.ggrid,tmp_02.kdata[nP,nT,nWN,:],label=r'$\widetilde{\nu}=0.5~cm^{-1}$')
    plt.plot(tmp_10.ggrid,tmp_10.kdata[nP,nT,nWN,:],label=r'$\widetilde{\nu}=0.1~cm^{-1}$')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.xlabel('Cumulative probability $g$')
    plt.ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    plt.show()
    
def plot_kcoeff_CPS_EH(band):
    pt_pair = 34
    ds = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k10.nc')
    k_coeff = ds['kopt'][band,pt_pair,:10]
    weights = ds['w_k'][band,:10]
    
    k_coeff_sorted=[]
    g_values=[0]
    
    for k in range(len(k_coeff)):
        k_min = min(k_coeff)
        c = 0
        while k_coeff[c] != k_min :
            c+=1
        k_coeff_sorted.append(k_min)
        k_coeff[c]=10**10
        g_values.append(g_values[-1]+weights[c])
    del(g_values[0])
    plt.plot(g_values,k_coeff_sorted,label='ExoMol',color="red")
    
    k_coeff = tmp_10.kdata[5,4,band,:]
    g_values = tmp_10.ggrid
    plt.plot(g_values,k_coeff,'--',label='ExoMol (Exo_k)',color="orange")
    
    
    dsH = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_HITRAN_k10.nc')
    k_coeffH = dsH['kopt'][band,pt_pair,:10]
    weightsH = dsH['w_k'][band,:10]
    k_coeff_sortedH=[]
    g_valuesH=[0]
    
    for k in range(len(k_coeffH)):
        k_min = min(k_coeffH)
        c = 0
        while k_coeffH[c] != k_min :
            c+=1
        k_coeff_sortedH.append(k_min)
        k_coeffH[c]=10**10
        g_valuesH.append(g_valuesH[-1]+weightsH[c])
    del(g_valuesH[0])
    plt.plot(g_valuesH,k_coeff_sortedH,label='HITRAN',color="black")
    
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.title(r'$(P,T)=(10~bar,2000~K)$, band : '+str(band)+r', $\widetilde{\nu} = $'+str(wnedges[band])+' $cm^{-1}$')
    plt.xlabel('Cumulative probability $g$')
    plt.ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    plt.show() 
    
def plot_kcoeff_CPS_EH_20(band):
    pt_pair = 24
    ds = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k20.nc')
    k_coeff = ds['kopt'][band,pt_pair,:20]
    weights = ds['w_k'][band,:20]
    
    k_coeff_sorted=[]
    g_values=[0]
    
    for k in range(len(k_coeff)):
        k_min = min(k_coeff)
        c = 0
        while k_coeff[c] != k_min :
            c+=1
        k_coeff_sorted.append(k_min)
        k_coeff[c]=10**10
        g_values.append(g_values[-1]+weights[c])
    del(g_values[0])
    plt.plot(g_values,k_coeff_sorted,label='ExoMol',color="red")
    print('g_values :',g_values,'\n\n','k_coeff :',k_coeff_sorted)
    
    k_coeff = tmp_10.kdata[4,0,band,:]
    g_values = tmp_10.ggrid
    plt.plot(g_values,k_coeff,'--',label='ExoMol (Exo_k)',color="orange")
    
    
    dsH = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_HITRAN_k20.nc')
    k_coeffH = dsH['kopt'][band,pt_pair,:20]
    weightsH = dsH['w_k'][band,:20]
    k_coeff_sortedH=[]
    g_valuesH=[0]
    
    for k in range(len(k_coeffH)):
        k_min = min(k_coeffH)
        c = 0
        while k_coeffH[c] != k_min :
            c+=1
        k_coeff_sortedH.append(k_min)
        k_coeffH[c]=10**10
        g_valuesH.append(g_valuesH[-1]+weightsH[c])
    del(g_valuesH[0])
    plt.plot(g_valuesH,k_coeff_sortedH,label='HITRAN',color="black")
    
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.title(r'$(P,T)=(1~bar,300~K)$, band : '+str(band)+r', $\widetilde{\nu} = $'+str(wnedges[band])+' $cm^{-1}$')
    plt.xlabel('Cumulative probability $g$')
    plt.ylabel('absorption coefficient $\kappa~[m^2.kg^{-1}]$')
    plt.show()

def plot_comp_spectra():
    pt_pair = 24
    mean_kE,mean_kH = [],[]
    for band in range(len(wnedges)-1):
        dsE = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k10.nc')
        k_coeffE = dsE['kopt'][band,pt_pair,:10]
        weightsE = dsE['w_k'][band,:10]
        mean_kE.append(np.sum(np.multiply(k_coeffE,weightsE)))
        
        dsH = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_HITRAN_k10.nc')
        k_coeffH = dsH['kopt'][band,pt_pair,:10]
        weightsH = dsH['w_k'][band,:10]
        mean_kH.append(np.sum(np.multiply(k_coeffH,weightsH)))
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    ax1.plot(wnedges[:-1],mean_kE,label='ExoMol',color='red')
    ax1.plot(wnedges[:-1],mean_kH,label='HITRAN',color='black')
    
    plt.yscale('log')
    plt.grid(visible=True)
    plt.legend() 
    ax1.set_xlim(-50,550)
    ax1.set_ylim(2e-4,200)
    
    y = np.arange(2e-4,2e+2,1e-4)
    #ax1.fill_betweenx(y,x1=500,x2=1300,color='yellow')
    #ax1.fill_betweenx(y,x1=2200,x2=2900,color='yellow')
    
    
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    ax1.set_ylabel('mean value absorption coefficient $\kappa~[m^2.kg^{-1}]$',fontsize=13)
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('wavelength '+r'$\lambda~[\mu m]$',fontsize=13)
        
    wl_ticks=[round(10**4/wn,4) for wn in ax1.get_xticks()]
    
    ax2.set_xticks([10**4/x for x in wl_ticks])
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(wl_ticks,)
    
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    
    plt.title(r'$(P,T)=(10~bar,500~K)$')
    plt.show()
    
def plot_comp_spectra_trans(pt_pair):
    mean_TEC,mean_TEE,mean_TH = [],[],[]
    inv = [-1 for k in range(20)]
    for band in range(len(wnedges)-1):
        dsEC = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k20_2306.nc')
        k_coeffEC = dsEC['kopt'][band,pt_pair,:20]
        weightsEC = dsEC['w_k'][band,:20]
        mean_TEC.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEC)),weightsEC)))
         
        dsH = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_HITRAN_k20.nc')
        k_coeffH = dsH['kopt'][band,pt_pair,:20]
        weightsH = dsH['w_k'][band,:20]
        mean_TH.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffH)),weightsH)))
        
        k_coeffEE = tmp_10.kdata[4,0,band,:]
        gvaluesEE = tmp_10.ggrid
        weightsEE = [gvaluesEE[0]]
        for i in range(1,len(gvaluesEE)):
            weightsEE.append(gvaluesEE[i]-gvaluesEE[i-1])
        mean_TEE.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEE)),weightsEE)))
    
    
    diffEC_H  = [mean_TEC[k]-mean_TH[k] for k in range(len(mean_TEC))]
    diffEC_EE = [mean_TEE[k]-mean_TH[k] for k in range(len(mean_TEC))]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
     
    ax1.plot(wnedges[:-1],mean_TEC,label='ExoMol Ccorr_k',color='red')
    ax1.plot(wnedges[:-1],mean_TEE,label='ExoMol Exo_k',color='orange')
    ax1.plot(wnedges[:-1],mean_TH,label='HITRAN',color='black')
    ax1.set_xlim(-10,500)
    
    ax2.plot(wnedges[:-1],diffEC_H,label='Ccorr_k',color='red')
    ax2.plot(wnedges[:-1],diffEC_EE,label='ExoMol',color='orange')
    ax2.set_xlim(-10,500)

    plt.grid(visible=True)
    ax1.legend() 
    ax2.legend()
    #ax1.set_xlim(400,3000)
    ax1.set_ylim(-1e-2,1.01)
    
    y = np.arange(0,1,1e-2)
    #ax1.fill_betweenx(y,x1=500,x2=1300,color='yellow')
    #ax1.fill_betweenx(y,x1=2200,x2=2900,color='yellow')
    
    ax2.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    ax2.set_ylabel('difference',fontsize=13)
     
    ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
    ax1.set_ylabel('transmission coefficient $T$',fontsize=13)
     
    ax3 = ax1.twiny()
    ax3.set_xlabel('wavelength '+r'$\lambda~[\mu m]$',fontsize=13)
         
    wl_ticks=[round(10**4/wn,4) for wn in ax1.get_xticks()]
     
    ax3.set_xticks([10**4/x for x in wl_ticks])
    ax3.set_xbound(ax1.get_xbound())
    ax3.set_xticklabels(wl_ticks,)
     
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax3.tick_params(axis='both', which='major', labelsize=12)
     
    plt.title(r'$(P,T)=(1~bar,300~K)$')
    plt.show()   
    
def plot_comp_grid():
    dsEC = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k20_tol.nc')
    dsH  = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_HITRAN_k20.nc')
    mean_diffEC, mean_diffEE = np.zeros((10,6)),np.zeros((10,6))
    for pt_pair in range(60):
        t_index,p_index=pt_pair%6,pt_pair//6
        TEC,TEE,TH = [],[],[]
        inv = [-1 for k in range(20)]
        for band in range(len(wnedges)-1):
             
            k_coeffEC = dsEC['kopt'][band,pt_pair,:20]
            weightsEC = dsEC['w_k'][band,:20]
            TEC.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEC)),weightsEC)))
          
         
            k_coeffH = dsH['kopt'][band,pt_pair,:20]
            weightsH = dsH['w_k'][band,:20]
            TH.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffH)),weightsH)))
         
            k_coeffEE = tmp_10.kdata[p_index,t_index,band,:]
            gvaluesEE = tmp_10.ggrid
            weightsEE = [gvaluesEE[0]]
            for i in range(1,len(gvaluesEE)):
                weightsEE.append(gvaluesEE[i]-gvaluesEE[i-1])
            TEE.append(np.sum(np.multiply(np.exp(np.multiply(inv,k_coeffEE)),weightsEE)))
            
        diffEC_H  = np.abs([TEC[k]-TH[k] for k in range(len(TEC))])
        diffEE_H  = np.abs([TEE[k]-TH[k] for k in range(len(TEE))])
         
        mean_diffEC[p_index,t_index] = np.sum([diffEC_H[k]*(wnedges[k+1]-wnedges[k]) for k in range(len(TEC))])/wnedges[-1]
        mean_diffEE[p_index,t_index] = np.sum([diffEE_H[k]*(wnedges[k+1]-wnedges[k]) for k in range(len(TEE))])/wnedges[-1]
      
    p   = np.array(P_grid_bis)
    t   = np.array(T_grid_bis)
    T,P = np.meshgrid(t,p)
    #plt.pcolormesh(T,P,mean_diffEC)
    plt.pcolormesh(T,P,mean_diffEE)
    
    plt.yscale('log')
    plt.ylabel('Pressure $P~[bar]$')
    plt.xlabel('Temperature $T~[K]$')
    plt.ylim(5e-5,2e2)
    plt.colorbar()
    plt.show()   
    
    


####################################### HITRAN ############################################################

dat_home        = "/proj/bolinc/users/x_anfey/spectral_files/"
dir_hitran      = dat_home+"dat_hitran/"

exec_file_name = "corr_k_HITRAN.sh"
f = open(exec_file_name, "w+")
f.write('Ccorr_k -F '+'pt_file'+' -D '+dir_hitran+'h2o_data.par -R 1 318'+' -c 3000.0 -i '+'0.1'+' -l 1 1.0e1 -b 1.5e-3 -s '+'sp_b318_ExoMol'+' +p -n 20 -lk -k -o output_HITRAN_k20 -m h2o_l'+'318'+'_lm_k20 -L h2o_lbl_lwf_pt48_k20.nc -sm h2o_l'+'318'+'_1-'+'318'+'map_k20.nc -np '+'1'+'\n')
f.close()
os.chmod(exec_file_name,0o777)

os.system('./'+exec_file_name)

################################## CIA MT_CKD ##########################################################

with open('pt_cia',"a") as file:
    file.write('*PTVAL'+'\n')
    for T in T_grid_bis:
        file.write(' '+str(T))
    file.write('\n')
    file.write('*END')

dat_home     = "/proj/bolinc/users/x_anfey/spectral_files/"
dir_hitran   = dat_home+"dat_hitran/"
dir_continua = dat_home+"dat_continua/"

exec_file_name = "corr_k_CIA_ExoMol.sh"
f = open(exec_file_name, "w+")
f.write('Ccorr_k -F '+'pt_cia'+' -R 1 318'+' -c 2500.0 -i '+'10.0'+' -ct 1 1 10.0 -t 1.0e-3 -e '+dir_continua+'mt_ckd_v3.0_s296 '+dir_continua+'mt_ckd_v3.0_s260 -k -s '+'sp_b318_HITRAN'+' +p -lk -o h2o-h2o_lc_HITRAN -m h2o-h2o_lcm_HITRAN -L 1H2-16O_ExoMol_CIA.nc'+ '\n')
f.close()
os.chmod(exec_file_name,0o777)

os.system('./'+exec_file_name)



# create a netCDF file
ncfile = netCDF4.Dataset("1H2-16O_ExoMol_CIA.nc","w",format='NETCDF3_CLASSIC',clobber='true')

# create the dimensions for the netCDF file
ncfile.createDimension('scalar',1)
ncfile.createDimension('nu',305001)
ncfile.createDimension('pt_pair',len(tmp_10.tgrid))

# create the variables
nu = ncfile.createVariable('nu', 'f8', ('nu',))
kabs = ncfile.createVariable('kabs','f8',('pt_pair','nu',))
t_calc = ncfile.createVariable('t_calc','f4',('pt_pair',))
p_calc = ncfile.createVariable('p_calc','f4',('pt_pair',))

# assign the data
with open(path_in+'1H2-16O_'+T_grid_str[0]+'K_Voigt_P'+P_grid_str[0]+'_cstep_10hires.xsec',"r") as f:
    WN = list()
    for line in f:
        data = line.split()
        WN.append(data[0])
nu[:] = [float(wn)*100 for wn in WN]

t_calc[:] = (tmp_10.tgrid.tolist())
p_calc[:] = [1e5 for k in range(len(tmp_10.tgrid))]


for i in range(len(tmp_10.tgrid)):
    with open(path_in+'1H2-16O_'+T_grid_str[i]+'K_Voigt_P'+'001'+'_cstep_10hires.xsec',"r") as f:
        ACS = list()
        for line in f:
            data = line.split()
            ACS.append(data[1])
        kabs[i,:] = [float(acs) for acs in ACS]

# assign attributes
#nu.title = 'wavenumber'
#nu.long_name = 'wavenumber'
#nu.units = 'm-1'
nu.step = nu[-1]-nu[-2]

#p_calc.title = 'pressure'
#p_calc.long_name = 'pressure'
#p_calc.units = 'Pa'

#t_calc.title = 'temperature'
#t_calc.long_name = 'temperature'
#t_calc.units = 'K'

#kabs.title = 'absorption'
#kabs.long_name = 'absorption'
#kabs.units = 'm2 kg-1'
        
ncfile.close()

############################# modify a netCDF file with k-coefficient from ExoMol CCorr_k ###############

dsEC = netCDF4.Dataset('/proj/bolinc/users/x_anfey/ExoMol/output_ExoMol_k20_bis.nc','r+')
for band in range(13):
    dsEC['w_k'][band,1] = dsEC['w_k'][band,0] + dsEC['w_k'][band,1]
    dsEC['w_k'][band,0] = 0 
dsEC.close()


################################ Ray ####################################################################

"""
#Hitran field locations
fieldLengths = [2,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7]
Sum = 0
fieldStart = [0]
for length in fieldLengths:
    Sum += length
    fieldStart.append(Sum)   
iso = 1
waveNum = 2
lineStrength = 3
airWidth =  5
selfWidth = 6
Elow = 7
TExp = 8
#
#
#Total internal partition functions (or generic approximations).
#These are used in the temperature scaling of line strength.
#The generic partition functions are OK up to about 500K
#(somewhat less for CO2)
def QGenericLin(T): #Approx. for linear molecules like CO2
    return T
def QGenericNonLin(T): #Approx for nonlinear molecules like H2O
    return T**1.5

molecules = {} #Start an empty dictionary
molecules['H2O'] = [1,18.,QGenericNonLin] #Entry is [molecule number,mol wt,partition fn]
molecules['CO2'] = [2,44.,QGenericLin]
molecules['CO'] = [5,28.,QGenericLin]
molecules['CH4'] = [6,16.,QGenericNonLin]
molecules ['HF'] = [14,20.,QGenericLin]
molecules['HCl'] = [15,36.,QGenericLin]
molecules['SO2'] = [9,64.,QGenericNonLin]


def computeAbsorption(p,T,dWave,numWidths = 1000.):
    waveGrid = [waveStart + dWave*i for i in range(int((waveEnd-waveStart)/dWave))]
    waveGrid = np.array(waveGrid)
    absGrid = np.zeros(len(waveGrid),float)

    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        gam = gamList[i]*(p/1.013e5)*(296./T)**TExpList[i]
        #Temperature scaling of line strength
        Tfact = math.exp(-100.*(h*c/k)*ElowList[i]*(1/T-1/296.))
        #The following factor is usually pretty close to unity
        #for lines that aren't far from the peak of the Planck spectrum
        #for temperature T, but it can become important on the low frequency
        #side, and is easy to incorporate.
        Tfact1 = (1.- math.exp(-100.*(h*c/k)*n/T))/ \
              (1.- math.exp(-100.*(h*c/k)*n/296.))   
        #The following has the incorrect algebraic prefactor used in
        #  the original version of PyTran (see Errata/Updates document)
        #S = sList[i]*(T/296.)**TExpList[i]*Tfact
        #The following is the corrected version, including also the
        #  low frequency factor Tfact1
        #S = sList[i]*(Q(296.)/Q(T))*TExpList[i]*Tfact*Tfact1
        #Preceding line didn't delete "*TExpList" factor.  Results now
        #checked against LMD kspectrum code, for Lorentz line case
        #-->Corrected on 6/10/2013
        S = sList[i]*(Q(296.)/Q(T))*Tfact*Tfact1
        #
        iw = int(len(waveGrid)*(n-waveStart)/(waveEnd-waveStart))
        nsum = int(numWidths*gam/dWave)
        i1 = max(0,iw-nsum)
        i2 = min(len(waveGrid)-1,iw+nsum)
        if i2>0:
            dn = np.arange(i1-iw,i2-iw)*dWave
            abs = S*gam/(math.pi*( dn**2 + gam**2))
            absGrid[i1:i2] += abs
    return waveGrid,absGrid

def get(line,fieldNum):
    return line[fieldStart[fieldNum]:fieldStart[fieldNum]+fieldLengths[fieldNum]]


def loadSpectralLines(molName):
    global waveList,sList,gamList,gamAirList,gamSelfList,ElowList,EupList,TExpList,Q
    molNum = molecules[molName][0]
    Q = molecules[molName][2] #Partition function for this molecule
    file = '/proj/bolinc/users/x_anfey/spectral_files/dat_hitran/h2o_data.par'
    f = open(file)
    waveList = []
    sList = []
    gamList =[]
    gamAirList = []
    gamSelfList = []
    TExpList = []
    ElowList = [] #Lower state energy
    #
    count = 0
    line = f.readline()
    while len(line)>0:
        count += 1
        isoIndex = int(get(line,iso))
        n = float(get(line,waveNum))
        S = float(get(line,lineStrength))
        El = float(get(line,Elow))
        #Convert line strength to (m**2/kg)(cm**-1) units
        #The cm**-1 unit is there because we still use cm**-1
        #as the unit of wavenumber, in accord with standard
        #practice for IR
        #
        #**ToDo: Put in correct molecular weight for the
        #        isotope in question.  
        S = .1*(N_avogadro/molecules[molName][1])*S
        gamAir = float(get(line,airWidth))
        gamSelf = float(get(line,selfWidth))
        TemperatureExponent = float(get(line,TExp))
        if  isoIndex == 1:
            waveList.append(n)
            sList.append(S)
            #gamList.append(gamSelf) #Put in switch to control this
            gamAirList.append(gamAir)
            gamSelfList.append(gamSelf)
            ElowList.append(El)
            TExpList.append(TemperatureExponent)
        #Read the next line, if there is one
        line = f.readline()
    #Convert the lists to numpy/Numeric arrays
    waveList = np.array(waveList)
    sList = np.array(sList)
    #gamList = numpy.array(gamList)
    gamAirList = np.array(gamAirList)
    gamSelfList = np.array(gamSelfList)
    ElowList =  np.array(ElowList)
    EupList = ElowList+waveList
    TExpList = np.array(TExpList)
    f.close()


######################################################################################################
####################################### MAIN SCRIPT ##################################################
######################################################################################################

h         = 6.62607015e-34
c         = 2.88792458e+8
k         = 1.380649e-23
N_avogadro= 6.0221408e+23

mult_factor = 3.3427961152976805e+21

waveStart = 1
waveEnd   = 2500

molName   = 'H2O'
loadSpectralLines(molName)

p         = 1.0e5
T         = 300.
gamList   = gamSelfList
dWave     = .1*(p/1.0e5)
nWidths   = 1000.
waveGrid,absGrid = computeAbsorption(p,T,dWave,nWidths)
"""
############

filename_01='1H2-16O_0300K_Voigt_P001_cstep_10hires.xsec'
path_ExoMol = '/proj/bolinc/users/x_anfey/ExoMol/xsec_cstep_10_files/'
wn_and_xs=np.fromfile(path_ExoMol+filename_01,count=-1, sep=' ')
wn,xs=[],[]
for i in range(len(wn_and_xs)//2):
    wn.append(wn_and_xs[2*i])
    xs.append(mult_factor*wn_and_xs[2*i+1])

############


fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(waveGrid,absGrid,label='HITRAN',color='black')
ax1.plot(wn,xs,label='ExoMol',color='red')
ax1.grid(True)
ax1.legend()
ax1.set_xlim(0,2500)
ax1.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
ax1.set_ylabel(r'absorption coefficient $\kappa~[m^2.kg^{-1}]$',fontsize=13)

ax2.plot(waveGrid,np.exp(-1*absGrid),label='HITRAN',color='black')
ax2.plot(wn,np.exp(-1*np.array(xs)),label='ExoMol',color='red')
ax2.grid(True)
ax2.legend()
ax2.set_xlim(0,2500)
ax2.set_xlabel('wavenumber '+r'$\widetilde{\nu}~[cm^{-1}]$',fontsize=13)
ax2.set_ylabel(r'transmission coefficient T$',fontsize=13)

plt.show()

"""