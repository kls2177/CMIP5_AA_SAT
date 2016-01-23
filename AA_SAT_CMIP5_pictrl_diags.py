#AA_SAT_CMIP5_pictrl_diags
#This module reads in monthly 2D piCTRL CMIP5 data, calculates seasonal averages, interpolates onto 2x2 lat-lon grid, area averages and calculates linear trends
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline

workdir = '/home/ksmith/Insight'
datdir  = 'XDIRX'
#datdir = '/Users/karensmith/Desktop'
varname = 'XVARX'
#varname = 'tas'
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
    var     = nc.variables[varname][:]
    lat_mod = nc.variables['lat'][:]
    lon_mod = nc.variables['lon'][:]
    if i == 0:
        var_month = np.squeeze(var)
    elif i == 1:
        var_month = np.stack((var_month,np.squeeze(var)),axis=-1)
    else:
        var = np.squeeze(var)
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

#do area averages
g = np.load(workdir + '/AA_masks.npz')
AA_mask = g['AA_mask'].T
E_AA_mask = g['E_AA_mask'].T
W_AA_mask = g['W_AA_mask'].T

from area_diags import area_avg
p = area_avg.averages()
AA_avg_ann = p.AAavg(var=var_ann_interp, lat=lat, lon=lon, mask=AA_mask)
AA_avg_djf = p.AAavg(var=var_djf_interp, lat=lat, lon=lon, mask=AA_mask)
AA_avg_mam = p.AAavg(var=var_mam_interp, lat=lat, lon=lon, mask=AA_mask)
AA_avg_jja = p.AAavg(var=var_jja_interp, lat=lat, lon=lon, mask=AA_mask)
AA_avg_son = p.AAavg(var=var_son_interp, lat=lat, lon=lon, mask=AA_mask)

E_AA_avg_ann = p.AAavg(var=var_ann_interp, lat=lat, lon=lon, mask=E_AA_mask)
E_AA_avg_djf = p.AAavg(var=var_djf_interp, lat=lat, lon=lon, mask=E_AA_mask)
E_AA_avg_mam = p.AAavg(var=var_mam_interp, lat=lat, lon=lon, mask=E_AA_mask)
E_AA_avg_jja = p.AAavg(var=var_jja_interp, lat=lat, lon=lon, mask=E_AA_mask)
E_AA_avg_son = p.AAavg(var=var_son_interp, lat=lat, lon=lon, mask=E_AA_mask)

W_AA_avg_ann = p.AAavg(var=var_ann_interp, lat=lat, lon=lon, mask=W_AA_mask)
W_AA_avg_djf = p.AAavg(var=var_djf_interp, lat=lat, lon=lon, mask=W_AA_mask)
W_AA_avg_mam = p.AAavg(var=var_mam_interp, lat=lat, lon=lon, mask=W_AA_mask)
W_AA_avg_jja = p.AAavg(var=var_jja_interp, lat=lat, lon=lon, mask=W_AA_mask)
W_AA_avg_son = p.AAavg(var=var_son_interp, lat=lat, lon=lon, mask=W_AA_mask)

#calculate linear trends for 49-year and 27-year chunks
year = np.arange(1,AA_avg_ann.size+1)
year_djf = year[0:year.size-2]
n = 27;
AA_trend_ann = np.zeros([AA_avg_djf.size - n])
AA_trend_jja = np.zeros([AA_avg_djf.size - n])
AA_trend_son = np.zeros([AA_avg_djf.size - n])
AA_trend_mam = np.zeros([AA_avg_djf.size - n])
AA_trend_djf = np.zeros([AA_avg_djf.size - n])
E_AA_trend_ann = np.zeros([AA_avg_djf.size - n])
E_AA_trend_jja = np.zeros([AA_avg_djf.size - n])
E_AA_trend_son = np.zeros([AA_avg_djf.size - n])
E_AA_trend_mam = np.zeros([AA_avg_djf.size - n])
E_AA_trend_djf = np.zeros([AA_avg_djf.size - n])
W_AA_trend_ann = np.zeros([AA_avg_djf.size - n])
W_AA_trend_jja = np.zeros([AA_avg_djf.size - n])
W_AA_trend_son = np.zeros([AA_avg_djf.size - n])
W_AA_trend_mam = np.zeros([AA_avg_djf.size - n])
W_AA_trend_djf = np.zeros([AA_avg_djf.size - n])

for r in range(AA_avg_djf.size - n):
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(AA_avg_ann[r:r+(n-1)]),1)
    AA_trend_ann[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(AA_avg_jja[r:r+(n-1)]),1)
    AA_trend_jja[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(AA_avg_son[r:r+(n-1)]),1)
    AA_trend_son[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(AA_avg_mam[r:r+(n-1)]),1)
    AA_trend_mam[r] = p1[0]*10

    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(E_AA_avg_ann[r:r+(n-1)]),1)
    E_AA_trend_ann[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(E_AA_avg_jja[r:r+(n-1)]),1)
    E_AA_trend_jja[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(E_AA_avg_son[r:r+(n-1)]),1)
    E_AA_trend_son[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(E_AA_avg_mam[r:r+(n-1)]),1)
    E_AA_trend_mam[r] = p1[0]*10

    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(W_AA_avg_ann[r:r+(n-1)]),1)
    W_AA_trend_ann[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(W_AA_avg_jja[r:r+(n-1)]),1)
    W_AA_trend_jja[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(W_AA_avg_son[r:r+(n-1)]),1)
    W_AA_trend_son[r] = p1[0]*10
    p1 = np.polyfit(year[r:r+(n-1)],np.squeeze(W_AA_avg_mam[r:r+(n-1)]),1)
    W_AA_trend_mam[r] = p1[0]*10

r = 1
for r in range(AA_avg_djf.size - n):
    p1 = np.polyfit(year_djf[r:r+(n-1)],np.squeeze(AA_avg_djf[r:r+(n-1)]),1)
    AA_trend_djf[r] = p1[0]*10
    p1 = np.polyfit(year_djf[r:r+(n-1)],np.squeeze(AA_avg_djf[r:r+(n-1)]),1)
    E_AA_trend_djf[r] = p1[0]*10
    p1 = np.polyfit(year_djf[r:r+(n-1)],np.squeeze(AA_avg_djf[r:r+(n-1)]),1)
    W_AA_trend_djf[r] = p1[0]*10

#concatenate seasons together for each region
AA_trends = np.stack((AA_trend_ann,AA_trend_mam,AA_trend_jja,AA_trend_son,AA_trend_djf),axis=-1)
E_AA_trends = np.stack((E_AA_trend_ann,E_AA_trend_mam,E_AA_trend_jja,E_AA_trend_son,E_AA_trend_djf),axis=-1)
W_AA_trends = np.stack((W_AA_trend_ann,W_AA_trend_mam,W_AA_trend_jja,W_AA_trend_son,W_AA_trend_djf),axis=-1)

#save interpolated seasonal averages in .npz file
from tempfile import TemporaryFile
npzname = workdir + '/' + varname + '_' + modname + '_' + memname + '_pictrl_27yrs'
dictname = {}
dictname[npzname] = TemporaryFile()
np.savez(npzname, lon=lon, lat=lat, AA_trends=AA_trends, E_AA_trends=E_AA_trends, W_AA_trends=W_AA_trends)
dictname[npzname].seek(0)
