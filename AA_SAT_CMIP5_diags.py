#AA_SAT_CMIP5_diags
#This module reads in monthly 2D CMIP5 data, calculates seasonal averages and interpolates onto 2x2 lat-lon grid
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline
#from matplotlib import pyplot as plt
#%matplotlib inline

workdir = '/home/ksmith/Insight'
datdir  = 'XDIRX'
#datdir = '/Users/karensmith/Desktop'
varname = 'XVARX'
#varname = 'psl'
memname = 'XMEMX'
#memname = 'r1i1p1'
modname = 'XMODELX'
#modname = 'CCSM4'
fileend = '.allyrs.nc'
month = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

#open netcdf files
for i in range(12):
    fname   = datdir + '/' + modname + '/' + memname + '/' + varname + '.' + month[i] + fileend
    #fname = datdir + '/' + varname + '.' + month[i] + fileend
    nc      = Dataset(fname)
    var     = nc.variables[varname][:,0]
    lat_mod = nc.variables['lat'][:]
    lon_mod = nc.variables['lon'][:]
    if i == 0:
        var_month = var
    elif i == 1:
        var_month = np.stack((var_month,var),axis=-1)
    else:
        var_month = np.concatenate((var_month,var[:,:,:,np.newaxis]),axis=-1)

#create seasonal averages
var_djf = (var_month[0:var_month.shape[0]-1,:,:,11] + var_month[1:var_month.shape[0],:,:,0] + var_month[0:var_month.shape[0]-1,:,:,1])/3
var_mam = (var_month[:,:,:,2] + var_month[:,:,:,3] + var_month[:,:,:,4])/3
var_jja = (var_month[:,:,:,5] + var_month[:,:,:,6] + var_month[:,:,:,7])/3
var_son = (var_month[:,:,:,8] + var_month[:,:,:,9] + var_month[:,:,:,10])/3
var_ann = np.mean(var_month, axis=3)

#load file with coarser 2D-grid
gridname = workdir + '/' + 'latlon_interp.npz'
g = np.load(gridname)
lat = g['lat']
lon = g['lon']

#interpolate onto coarser 2D-grid (DJF)
var_djf_interp = np.zeros([var_djf.shape[0],lat.size,lon.size])
for i in range(var_djf.shape[0]):
    ip = RectBivariateSpline(lat_mod, lon_mod, np.squeeze(var_djf[i,:,:])) 
    var_djf_interp[i,:,:] = ip(lat, lon)

#interpolate onto coarser 2D-grid (MAM, JJA, SON)
var_mam_interp = np.zeros([var_month.shape[0],lat.size,lon.size])
var_jja_interp = np.zeros([var_month.shape[0],lat.size,lon.size])
var_son_interp = np.zeros([var_month.shape[0],lat.size,lon.size])
var_ann_interp = np.zeros([var_month.shape[0],lat.size,lon.size])
for i in range(var_month.shape[0]):
    ip = RectBivariateSpline(lat_mod, lon_mod, np.squeeze(var_mam[i,:,:])) 
    var_mam_interp[i,:,:] = ip(lat, lon)
    ip = RectBivariateSpline(lat_mod, lon_mod, np.squeeze(var_jja[i,:,:])) 
    var_jja_interp[i,:,:] = ip(lat, lon)
    ip = RectBivariateSpline(lat_mod, lon_mod, np.squeeze(var_son[i,:,:])) 
    var_son_interp[i,:,:] = ip(lat, lon)
    ip = RectBivariateSpline(lat_mod, lon_mod, np.squeeze(var_ann[i,:,:])) 
    var_ann_interp[i,:,:] = ip(lat, lon)

#save interpolated seasonal averages in .npz file
from tempfile import TemporaryFile
npzname = workdir + '/' + varname + '_' + modname + '_' + memname
dictname = {}
dictname[npzname] = TemporaryFile()
np.savez(npzname, lon=lon, lat=lat, var_djf_interp=var_djf_interp, var_mam_interp=var_mam_interp, var_jja_interp=var_jja_interp, var_son_interp=var_son_interp, var_ann_interp=var_ann_interp)
dictname[npzname].seek(0)
