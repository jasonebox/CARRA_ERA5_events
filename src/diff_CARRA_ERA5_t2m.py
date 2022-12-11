#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 09:22:24 2022

@author: jason
"""

import numpy as np
import os
# from glob import glob
# from netCDF4 import Dataset
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl

# from datetime import datetime 
# import calendar
# import netCDF4
# from cftime import num2date, date2num
import matplotlib.pyplot as plt
import xarray as xr

ly='x'

res='l'
if ly=='p':res='h'


if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/CARRA/CARRA_ERA5_events'

th=2 # line thickness
formatx='{x:,.3f}' ; fs=16
plt.rcParams["font.size"] = fs
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.8
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'
plt.rcParams['grid.linewidth'] = th/2
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['figure.figsize'] = 12, 20

#%%
# CARRA West grid dims
ni = 1269 ; nj = 1069

map_version=1

if map_version:
    fn='./meta/CARRA/2.5km_CARRA_west_lat_1269x1069.npy'
    lat=np.fromfile(fn, dtype=np.float32)
    lat=lat.reshape(ni, nj)

    fn='./meta/CARRA/2.5km_CARRA_west_lon_1269x1069.npy'
    lon=np.fromfile(fn, dtype=np.float32)
    lon=lon.reshape(ni, nj)


    # latx=np.rot90(lat.T)
    # lonx=np.rot90(lon.T)
    offset=0
    lon=lon[offset:ni-offset,offset:nj-offset]
    lat=lat[offset:ni-offset,offset:nj-offset]
    ni-=offset*2
    nj-=offset*2
    # print(ni,nj)
    LLlat=lat[0,0]
    LLlon=lon[0,0]-360
    # print("LL",LLlat,LLlon)
    # print("UL",lat[ni-1,0],lon[ni-1,0]-360)
    lon0=lon[int(round(ni/2)),int(round(nj/2))]-360
    lat0=lat[int(round(ni/2)),int(round(nj/2))]
    # print("mid lat lon",lat0,lon0)
    
    URlat=lat[ni-1,nj-1]
    URlon=lon[ni-1,nj-1]
    # print("LR",lat[0,nj-1],lon[0,nj-1]-360)
    # print("UR",URlat,URlon)
    
    # m = Basemap(llcrnrlon=LLlon, llcrnrlat=LLlat, urcrnrlon=URlon, urcrnrlat=URlat, lat_0=72, lon_0=-36, resolution='l', projection='lcc')
    m = Basemap(llcrnrlon=LLlon, llcrnrlat=LLlat, urcrnrlon=URlon, urcrnrlat=URlat, lat_0=lat0, lon_0=lon0, resolution=res, projection='lcc')
    
    x, y = m(lat, lon)
    lon-=360

#%%

datex='2021-08-14T00'
ERA5_data=np.fromfile('/Users/jason/0_dat/ERA5/events/resampled/t2m/'+datex+'_1269x1069.npy', dtype=np.float16)
ERA5_data=ERA5_data.reshape(ni, nj)

CARRA_data=np.fromfile('/Users/jason/0_dat/CARRA/output/t2m_2021-08-14-00_1269x1069.npy', dtype=np.float16)
CARRA_data=CARRA_data.reshape(ni, nj)

dx=ERA5_data-CARRA_data
lo=-5 

map_version=0
ly='x'

if map_version==0:
    fig, ax = plt.subplots(figsize=(10,10))
else:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14,10),gridspec_kw={'width_ratios': [1, 0.5]})

plt.imshow(dx, cmap=mpl.colormaps['bwr'],vmin=lo,vmax=-lo)
plt.axis('off')
plt.title('ERA5 minus CARRA t2m '+datex)

# plt.imshow(CARRA_data)
# plt.imshow(ERA5_data)
clb=plt.colorbar()
clb.ax.tick_params(labelsize=fs) 
clb.ax.set_title('°C',fontsize=fs)