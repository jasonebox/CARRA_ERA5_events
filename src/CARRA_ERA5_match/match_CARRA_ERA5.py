#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Adrien Wehrlé, University of Zurich, Switzerland

"""

import pygrib
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

sys.path.append("/home/adrien/EO-IO/geomatcher")
import geomatcher.geomatcher as gm


def lon360_to_lon180(lon360):

    # reduce the angle
    lon180 = lon360 % 360

    # force it to be the positive remainder, so that 0 <= angle < 360
    lon180 = (lon180 + 360) % 360

    # force into the minimum absolute value residue class, so that -180 < angle <= 180
    lon180[lon180 > 180] -= 360

    return lon180


# %% read ERA5 example

grbs = pygrib.open(
    "/home/adrien/EO-IO/CARRA_ERA5_events/data/ERA5/tcwv/202206_3hourly_tcwv.grib"
)

grbs.seek(0)
for grb in grbs:
    print(grb)

selected_grb = grbs.select(name="Total column water vapour")[0]

era5_ex, era5_lats, era5_lons_raw = selected_grb.data()
era5_lons = lon360_to_lon180(era5_lons_raw)
# era5_lons = era5_lons_raw

era5_grid = np.dstack([era5_ex, era5_lons, era5_lats])

# %% read CARRA grid

ni = 1269
nj = 1069
fn = "/home/adrien/EO-IO/CARRA_rain/ancil/2.5km_CARRA_west_lat_1269x1069.npy"
lat = np.fromfile(fn, dtype=np.float32)
lat = lat.reshape(ni, nj)

fn = "./ancil/2.5km_CARRA_west_lon_1269x1069.npy"
lon = np.fromfile(fn, dtype=np.float32)
lon = lon.reshape(ni, nj)


fn = "/home/adrien/EO-IO/CARRA_rain/ancil/CARRA_W_elev_lat_lon.nc"

carra_ds = xr.open_dataset(fn)

carra_elev = np.array(carra_ds["z"])[::-1]
carra_lat = np.array(carra_ds["latitude"])[::-1]

carra_lon_raw = np.array(carra_ds["longitude"])[::-1]
# carra_lon = lon360_to_lon180(carra_lon_raw)
# carra_lon = carra_lon_raw

carra_lon_raw2 = np.fromfile(
    "/home/adrien/EO-IO/CARRA_ERA5_events/meta/CARRA/2.5km_CARRA_west_lon_1269x1069.npy"
)

carra_grid = np.dstack([carra_elev, carra_lon, carra_lat])

# %%
plt.figure()
ax1 = plt.subplot(211)
ax1.plot(era5_lons, era5_lats, alpha=0.1, color="gray")
ax1.plot(carra_grid[:, :, 0], carra_grid[:, :, 1], alpha=0.1, color="red")
ax2 = plt.subplot(212)
ax2.imshow(era5_ex, aspect="auto")

# %%

plt.figure()
ax1 = plt.subplot(211)
ax1.imshow(era5_lons, aspect="auto")
ax2 = plt.subplot(212)
ax2.imshow(carra_lon, aspect="auto")

# %% geomatch

indexes = gm.match_m2m(carra_grid, era5_grid, only_indexes=True)

# %% check results

plt.figure()
plt.imshow(era5_ex.flatten()[indexes], vmin=0, vmax=5, aspect="auto")
