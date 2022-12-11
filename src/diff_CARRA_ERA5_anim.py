#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 09:22:24 2022

@author: jason
"""

import numpy as np
import os
from glob import glob
# from netCDF4 import Dataset
# import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# from datetime import datetime 
# import calendar
# import netCDF4
# from cftime import num2date, date2num
import matplotlib.pyplot as plt
# import xarray as xr

ly='x'

if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/'

os.chdir(base_path)

th=2 # line thickness
formatx='{x:,.3f}' ; font_size=16
plt.rcParams["font.size"] = font_size
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

# res='l'
# if ly=='p':res='h'
res='h'

if map_version:
    fn='./meta/CARRA/2.5km_CARRA_west_lat_1269x1069.npy'
    lat=np.fromfile(fn, dtype=np.float32)
    lat=lat.reshape(ni, nj)

    fn='./meta/CARRA/2.5km_CARRA_west_lon_1269x1069.npy'
    lon=np.fromfile(fn, dtype=np.float32)
    lon=lon.reshape(ni, nj)


    # lat=np.rot90(lat.T)
    # lon=np.rot90(lon.T)
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
    
    m = Basemap(llcrnrlon=LLlon, llcrnrlat=LLlat, urcrnrlon=URlon, urcrnrlat=URlat, lat_0=lat0, lon_0=lon0, resolution=res, projection='lcc')
    
    x, y = m(lat, lon)
    lon-=360

#%%

def imshow_var(var,lo,hi,cmap,source,datex):
    plt.close()
    fig, ax = plt.subplots(figsize=(10,10))
    font_size=14
    plt.imshow(var, cmap=mpl.colormaps[cmap],vmin=lo,vmax=hi)
    plt.axis('off')
    plt.title(source+' '+datex)
    
    # plt.imshow(CARRA_data)
    # plt.imshow(ERA5_data)
    clb=plt.colorbar()
    clb.ax.tick_params(labelsize=font_size) 
    clb.ax.set_title(units,fontsize=font_size)
    DPI=200
    ly='p'
            
    if ly == 'p':
        figpath='/Users/jason/0_dat/ERA5_vs_CARRA/Figs/'+choice+'/'
        # figpath='/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/Figs/ERA5/'
        os.system('mkdir -p '+figpath)
        plt.savefig(figpath+source+'_'+str(datex)+'.png', bbox_inches='tight', pad_inches=0.04, dpi=DPI, facecolor='w', edgecolor='k')
        # plt.savefig(figpath+select_period+'JJA_'+hgt+'z_anom.eps', bbox_inches='tight')

def map_var(var,tit,lo,hi,cmap,source,datex):
    plt.close()
    fig, ax = plt.subplots(figsize=(10,10))
    # font_size=16
    pp=m.imshow(var, cmap = cmap,vmin=lo,vmax=hi) 
    # m.axis('off')
            # mine=m.contourf(x,y,plotvar,clevs,cmap=mcmap,linewidths=th/1.)
            # clevs=np.arange(elev0,elev1,delev)
            # plt.clim(lo,hi
            # du=1
            # if du:
            #     clb = plt.colorbar(fraction=0.046/2., pad=0.04)
            #     clb.ax.set_title(units,fontsize=font_size,c='k')
        
    m.drawcoastlines(color='k',linewidth=0.5)
    # clevs=np.arange(1,10)

    m.drawparallels([66.6],color='gray')
    m.drawparallels([60,70,80],dashes=[2,4],color='k')
    m.drawmeridians(np.arange(0.,420.,10.))

    # ax = plt.gca()     
    # lons, lats = m(lon, lat)
    # m.contour(lons, lats,dx,[0],linewidths=1,colors='grey')#,linestyle='dashed')
    plt.title(tit)

    # m.scatter(lons[plotvar_non_fuzzy==maxval2],lats[plotvar_non_fuzzy==maxval2], s=380, facecolors='none', edgecolors='m')
    # m.scatter(lons[plotvar_non_fuzzy==maxval],lats[plotvar_non_fuzzy==maxval], s=780, facecolors='none', edgecolors='k')
    # print(alllat,alllon)# m.scatter(lons[plotvar_non_fuzzy==minval],lats[plotvar_non_fuzzy==minval], s=380, facecolors='none', edgecolors='m')
 
    clb=plt.colorbar()
    clb.ax.tick_params(labelsize=font_size) 
    clb.ax.set_title(units,fontsize=font_size)
    DPI=200
    ly='p'
            
    if ly == 'p':
        figpath='/Users/jason/0_dat/ERA5_vs_CARRA/Figs/'+choice+'/'
        # figpath='/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/Figs/ERA5/'
        os.system('mkdir -p '+figpath)
        plt.savefig(figpath+source+'_'+str(search_date)+'.png', bbox_inches='tight', pad_inches=0.04, dpi=DPI, facecolor='w', edgecolor='k')
        # plt.savefig(figpath+select_period+'JJA_'+hgt+'z_anom.eps', bbox_inches='tight')


choices=['tp']
cmap='bwr'
choices=['tp']

for choice in choices:
    if choice=='tp':
        cmap='BrBG'
        lo=-20
        units='mm'
    if choice=='t2m':
        cmap='BrBG'
        lo=-5
        units='°C'
    CARRA_files=sorted(glob('/Users/jason/0_dat/CARRA/output/'+choice+'/*'))
    CARRA_files=np.array(CARRA_files)
    CARRA_date=[]
    for fn in CARRA_files:
        # print(fn.split('/')[-1][3:16])
        CARRA_date.append(fn.split('/')[-1][0:13])
    
    CARRA_date=np.array(CARRA_date)
    print(CARRA_date)
    
    dirs = sorted(glob('/Users/jason/0_dat/ERA5/events/resampled/'+choice+'/'+'*'))

    # for fn_index,fn in enumerate(dirs[10:16]):
    for fn_index,fn in enumerate(dirs):
        if fn_index<=1000:
            search_date=fn.split('/')[-1][0:13]
            v=np.where(CARRA_date==search_date)
            
            print(search_date)
            
            v=v[0]
            if np.size(v):
                v=v[0]
                print('yes',search_date)
        
                ERA5_data=np.fromfile(fn, dtype=np.float16)
                ERA5_data=ERA5_data.reshape(ni, nj)*3
                ERA5_data=np.flipud(ERA5_data)

                CARRA_data=np.fromfile(CARRA_files[v], dtype=np.float16)
                CARRA_data=CARRA_data.reshape(ni, nj)
                CARRA_data=np.flipud(CARRA_data)
                
                dx=ERA5_data-CARRA_data
                # dx=np.rot90(dx).T

                # imshow_var(CARRA_data,0,6,'viridis','CARRA',search_date)
                # imshow_var(ERA5_data,0,6,'viridis','ERA5',search_date)
                
                map_version=1
                ly='x'

                if map_version:
                    # map_var(dx,lo,-lo,cmap,'ERA5 minus CARRA '+choice,search_date)
                    # map_var(CARRA_data,'CARRA '+choice+' '+search_date,0,7,'Blues','CARRA',search_date)
                    # map_var(ERA5_data,'ERA5 '+choice+' '+search_date,0,7,'Blues','ERA5',search_date)
                    map_var(dx,'ERA5 minus CARRA, '+choice+' '+search_date,-7,7,'BrBG','diff',search_date)
                    
                # else:
                #     plt.imshow(dx, cmap=mpl.colormaps[cmap],vmin=lo,vmax=-lo)
                #     plt.axis('off')
                #     plt.title('ERA5 minus CARRA '+choice+' '+search_date)
                    
                    # plt.imshow(CARRA_data)
                    # plt.imshow(ERA5_data)
                    # clb=plt.colorbar()
                    # clb.ax.tick_params(labelsize=font_size) 
                    # clb.ax.set_title(units,fontsize=font_size)
    
                    DPI=200
                    ly='p'
                            
                    if ly == 'p':
                        figpath='/Users/jason/0_dat/ERA5_vs_CARRA/Figs/'+choice+'/'
                        # figpath='/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/Figs/ERA5/'
                        os.system('mkdir -p '+figpath)
                        plt.savefig(figpath+'diff_'+str(search_date)+'.png', bbox_inches='tight', pad_inches=0.04, dpi=DPI, facecolor='w', edgecolor='k')
                        # plt.savefig(figpath+select_period+'JJA_'+hgt+'z_anom.eps', bbox_inches='tight')
            
