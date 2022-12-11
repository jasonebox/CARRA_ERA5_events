#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gather CARRA or ERA5 grib data from CDS

followed by:
    /Users/jason/Dropbox/CARRA/CARRA_ERA5_events/src/ERA5/Resampling_ERA5_to_CARRA_code_from_extremes_list.py
    /Users/jason/Dropbox/CARRA/CARRA_ERA5_events/src/CARRA/extract_CARRA_rf_from_tpt2m_from_extremes_list.py

updated Dec 2022
@author: Jason Box, GEUS, jeb@geus.dk
"""
import cdsapi
import os
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import calendar
from datetime import timedelta

if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/'

os.chdir(base_path)


fn='/Users/jason/Dropbox/CARRA/CARRA_rainfall_study/stats/rf_extremes.csv'
os.system('ls -lF '+fn)
# os.system('head '+fn)
# os.system('tail '+fn)
df=pd.read_csv(fn)

# print(df.columns)
##%%
# minval=300 ;scaling_factor=3000
# v=np.where(df.maxlocalrate>minval)

v=np.where(df.Gt_overall>2.2)
# v=np.where(df.Gt_overall>2.7)
# v=v[0]
# print(v[0])

# print()
minval=np.min(df.maxlocalrate[v[0]]) ; scaling_factor=1500
# listx=[v[0][3]] # 2021
listx=[v[0][12]] # 2017
# listx=v[0]
# listx=[v[0][0]] # 2012a


model='CARRA'
# model='ERA5'
models=['CARRA','ERA5']

for model in models:
    if model=='CARRA':
        choices=['tpt2m']
    else:
        # choices=['tpsf','tzuv','tcwv']
        choices=['tpsf']
    
    for ii,i in enumerate(listx):
        print(model,df.date[i],df.Gt_overall[i],ii)
        ##%%
        timex=pd.to_datetime(df.date[i]) 
        # timex=pd.to_datetime('17/06/2022') # force this date
        # print(df.date[i])
        # timex=[]
    
        ddx=int(timex.strftime('%j'))
        ddx0=ddx-1 ; ddx1=ddx+1
        ymd=timex.strftime('%Y-%m-%d')
        year=timex.strftime('%Y')
        month=timex.strftime('%m')
        day=timex.strftime('%d')
    
        timex0=pd.to_datetime(df.date[i])- timedelta(days=1)
        ymd0=timex0.strftime('%Y-%m-%d')
        day_before=timex0.strftime('%d')
        
        timex1=pd.to_datetime(df.date[i])+ timedelta(days=1)
        ymd1=timex1.strftime('%Y-%m-%d')
        
        print(i,ymd0,ymd,ymd1,ddx,df.lon[i],df.lat[i],df.Gt_overall[i],df.maxlocalrate[i])
    
        for choice_index,choice in enumerate(choices):
        
            # path='/Users/jason/Dropbox/CARRA/CARRA_rain/'
            
            # os.chdir(path)
            # ofile="./ancil/events/"+varnams[j]+"_events.csv"
            # # read in events df
            # df=pd.read_csv(ofile)
            # print(df.columns)
            
            # path='./Figs/event/'+varnams[j]+'/'
            
            
            # files = sorted(glob(path+"*.png"), reverse=True)
            
            # numpy.savetxt("foo.csv", a, delimiter=",")
            
            # os.system('mkdir -p /Users/jason/0_dat/ERA5/events/')
            # os.system('mkdir -p /Users/jason/0_dat/ERA5/events//')
            
            # print(files)
                
            yearx=[]
            monx=[]
            dayx=[]
            Gt=[]
        
            day_before=str(int(day)-1).zfill(2)
            day_three=str(int(day)+1).zfill(2)
            # day_eight=str(int(day)+6).zfill(2)
            print(choice,year,month,day_before,day,day_three)
            
            opath='/Users/jason/0_dat/'+model+'/events/'+choice+'/'
            # opath='/Users/jason/Dropbox/CARRA/CARRA_ERA5_events/data_raw/ERA5/'+choice+'/'
            os.system('mkdir -p '+opath)
            
            ofile=opath+str(year)+str(month).zfill(2)+str(day_before).zfill(2)+'-'+day_three+'_3hourly_'+choice+'.grib'
        
            c = cdsapi.Client()
        
            if choice=='tcwv':
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'grib',
                        'variable': 'total_column_water_vapour',
                        'year': year,
                        'month': str(month).zfill(2),
                        'day': [
                            str(int(day_before)).zfill(2),
                            str(day).zfill(2),
                            str(int(day_three)).zfill(2),
                            ],
                        'time': [
                            '00:00', '03:00', '06:00',
                            '09:00', '12:00', '15:00',
                            '18:00', '21:00',
                            ],
                    },
                    ofile)
        
            if choice=='tzuv':        
                c.retrieve(
                    'reanalysis-era5-pressure-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'grib',
                        'variable': [
                            'temperature', 'u_component_of_wind', 'v_component_of_wind',
                        ],
                        'pressure_level': '850',
                        'year': year,
                        'month': str(month).zfill(2),
                        'day': [
                            str(int(day_before)).zfill(2),
                            str(day).zfill(2),
                            str(int(day_three)).zfill(2),
                            ],
                        'time': [
                            '00:00', '03:00', '06:00',
                            '09:00', '12:00', '15:00',
                            '18:00', '21:00',
                            ],
                    },
                    ofile)
        
            if choice=='tpsf':        
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'grib',
                        'variable': ['snowfall', 'total_precipitation',],
                        'year': year,
                        'month': str(month).zfill(2),
                        'day': [
                            str(int(day_before)).zfill(2),
                            str(day).zfill(2),
                            str(int(day_three)).zfill(2),
                            ],
                        'time': [
                            '00:00', '03:00', '06:00',
                            '09:00', '12:00', '15:00',
                            '18:00', '21:00',
                        ],
                    },
                ofile)
                
            if choice=='tpt2m':
                c.retrieve(
                    'reanalysis-carra-single-levels',
                    {
                        'format': 'grib',
                        'domain': 'west_domain',
                        'level_type': 'surface_or_atmosphere',
                        'variable': [
                            '2m_temperature', 'total_precipitation',
                        ],
                        'product_type': 'forecast',
                        'time': [
                            '00:00', '03:00', '06:00',
                            '09:00', '12:00', '15:00',
                            '18:00', '21:00',
                        ],
                        'year': year,
                        'month': str(month).zfill(2),
                        'day': [
                            str(int(day_before)).zfill(2),
                            str(day).zfill(2),
                            str(int(day_three)).zfill(2),
                            ],
                        'leadtime_hour': '3',
                    },
                    ofile)


