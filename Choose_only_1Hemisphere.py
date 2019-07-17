#!/usr/bin/python
# Choose only 1 hemisphere in case of global data
# Arun Rana
#
import netCDF4
import numpy as np
import sys
import scipy.stats
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
from make_cellarea import make_cellarea
from netCDF4 import Dataset

flat = Dataset('/home/arunr/UCL_BE/Diagnostics/tools/NOAA_OI_lat.nc')
lat = flat["lat"][:]
flon = Dataset('/home/arunr/UCL_BE/Diagnostics/tools/NOAA_OI_lon.nc')
lon = flon["lon"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] > 0:
      a[:][i] = 1
   else:
      a[:][i] = 0

