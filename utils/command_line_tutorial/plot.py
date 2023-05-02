# This code plots the selected netCDF file. It does not require the SOCRATES environment
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import seaborn as sns
import sys

CONST = 287. * 250. / (9.80665 * 1000.0) #

if __name__ == '__main__':
    if (len(sys.argv) > 2):
        filename = sys.argv[1]
        filename2 = sys.argv[2]
        dgs2 = Dataset(filename2)
    elif (len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        raise RuntimeError('please enter a file name')
 
dgs = Dataset(filename)
name = filename[filename.find('.') + 1:]

lon = dgs.variables['lon'][:]
lat = dgs.variables['lat'][:]
p = dgs.variables['plev'][:]

n_lon = len(lon)
n_lat = len(lat)
layers= len(p)

var = dgs.variables[name][:]
if (len(sys.argv) > 2):
    var2 = dgs2.variables[name][:]
    var = var - var2
    
try:
    width     = dgs.variables['bandwidth'][:]
    wl_short  = dgs.variables['wl_short'][:]
    wl_long   = dgs.variables['wl_long'][:]
    n_channel = len(width)    
except:
    vmean = np.zeros(layers)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[i, :, :]) / (n_lon * n_lat)
    ax1 = plt.figure().add_subplot(111)
    ax1.plot(vmean, -np.log(p/max(p))*CONST)
    ax1.set_title('Average profile')
    ax1.set_xlabel(name)
    ax1.set_ylabel('Approx height (km)')
else:
    fig=plt.figure()
    vmean = np.zeros(layers)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[:, i, :, :]) / (n_lon * n_lat)
    ax1 = fig.add_subplot(121)
    ax1.plot(vmean, -np.log(p/max(p))*CONST)
    ax1.set_title('Average profile')
    ax1.set_xlabel(name)
    ax1.set_ylabel('Approx height (km)')
    wn = 0.5e-2/wl_short + 0.5e-2/wl_long
    wl = 0.5e6*(wl_short + wl_long)
    toa_spec = np.zeros(n_channel)
    surf_spec = np.zeros(n_channel)
    for ch in range(0, n_channel):
        toa_spec[ch]  = np.sum(var[ch,0,       :,:])/(width[ch]*n_lon*n_lat)
        surf_spec[ch] = np.sum(var[ch,layers-1,:,:])/(width[ch]*n_lon*n_lat)
    ax2 = fig.add_subplot(122)
    ax2.plot(wl, toa_spec, color='blue', label='TOA')
    ax2.plot(wl, surf_spec, color='green', label='Surface')
    ax2.set_xscale('log')
    ax2.set_title('TOA & surface spectrum')
    ax2.set_xlabel('Wavelength (micron)')
    ax2.set_ylabel('Flux (Wm-2m-1)')
    plt.legend()

plt.tight_layout()
plt.show()
