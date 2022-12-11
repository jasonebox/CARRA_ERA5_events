#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:12:32 2022

@author: jason
"""

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-carra-pressure-levels',
    {
        'format': 'grib',
        'domain': 'west_domain',
        'variable': [
            'geometric_vertical_velocity', 'temperature', 'u_component_of_wind',
            'v_component_of_wind',
        ],
        'pressure_level': '850',
        'product_type': 'analysis',
        'time': [
            '00:00', '03:00', '06:00',
            '09:00', '12:00', '15:00',
            '18:00', '21:00',
        ],
        'year': '2022',
        'month': '06',
        'day': '17',
    },
    '/Users/jason/0_dat/CARRA_2022_June/20220617_3h_CARRA_UVWT.grib')