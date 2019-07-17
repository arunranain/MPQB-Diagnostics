#!/usr/bin/python
# Common sea ice diags
# Arun Rana
#https://matplotlib.org/basemap/api/basemap_api.html

import netCDF4
import numpy as np
import sys
import scipy.stats
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
from mpl_toolkits.basemap import Basemap
from seaice_commondiags import compute_extent
from seaice_commondiags import compute_area
from netCDF4 import Dataset

#Get Lat/Lon to plot since this is regridded 1*1 (you can use any interp dataset for that) data and 1992-2007 data as ESACCI_SSMI reference period
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
lon, lat = np.meshgrid(lon, lat)
del data

#NSIDC0051 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc') #lat 448 lon 304 (25*25 resolution)
conc = data.variables['siconc'][156:360]
del data
std_conc_NSIDC0051_NH = np.std(conc, axis=0)

fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,1)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image=m.pcolormesh(lon, lat, std_conc_NSIDC0051_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("NSIDC0051")

#OSISAF409a - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][156:360]
del data
std_conc_OSISAF409a_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,2)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_OSISAF409a_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("OSISAF409a")

#OSISAF450 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][24:228]
del data
std_conc_OSISAF450_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,3)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_OSISAF450_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("OSISAF450")

#NSIDC_G02202_v3 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc')
conc = data.variables['siconc'][36:240]
del data
std_conc_NSIDC_G02202_v3_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,4)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_NSIDC_G02202_v3_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("NSIDC_G02202_v3")

#ICDC_ASI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc')
conc = data.variables['siconc'][0:204]
del data
std_conc_ICDC_ASI_SSMI_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,5)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_ICDC_ASI_SSMI_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("ICDC_ASI_SSMI")

#Hadisst1 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc = data.variables['siconc'][1464:1668]
del data
std_conc_Hadisst1_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,6)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_Hadisst1_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("Hadisst1")

#Hadisst2 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/interp/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc = data.variables['siconc'][1704:1908]
del data
std_conc_Hadisst2_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,7)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_Hadisst2_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("Hadisst2")


#ESACCI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/interp/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc') 
conc = data.variables['siconc'][0:204]
del data
std_conc_ESACCI_SSMI_NH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,8)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_ESACCI_SSMI_NH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("ESACCI_SSMI")

#colorbars and title
cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image, cax, orientation='vertical')
plt.suptitle("Sea-Ice Concentration - Standard Error (1992-2007) - Northern Hemisphere")
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Std_nh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Std_nh.pdf', format='pdf')
plt.close(fig)

######################################################################################################################################################################################################################################################################
#NSIDC0051 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc') #lat 448 lon 304 (25*25 resolution)
conc = data.variables['siconc'][156:360]
del data
std_conc_NSIDC0051_SH = np.std(conc, axis=0)

fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,1)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image=m.pcolormesh(lon, lat, std_conc_NSIDC0051_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("NSIDC0051")

#OSISAF409a - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][156:360]
del data
std_conc_OSISAF409a_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,2)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for Southh
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_OSISAF409a_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("OSISAF409a")

#OSISAF450 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][24:228]
del data
std_conc_OSISAF450_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,3)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_OSISAF450_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("OSISAF450")

#NSIDC_G02202_v3 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc')
conc = data.variables['siconc'][36:240]
del data
std_conc_NSIDC_G02202_v3_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,4)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_NSIDC_G02202_v3_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("NSIDC_G02202_v3")

#ICDC_ASI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc')
conc = data.variables['siconc'][0:204]
del data
std_conc_ICDC_ASI_SSMI_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,5)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_ICDC_ASI_SSMI_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("ICDC_ASI_SSMI")

#Hadisst1 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc = data.variables['siconc'][1464:1668]
del data
std_conc_Hadisst1_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,6)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_Hadisst1_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("Hadisst1")

#Hadisst2 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/interp/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc = data.variables['siconc'][1704:1908]
del data
std_conc_Hadisst2_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,7)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_Hadisst2_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("Hadisst2")

#ESACCI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/interp/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc') 
conc = data.variables['siconc'][0:204]
del data
std_conc_ESACCI_SSMI_SH = np.std(conc, axis=0)

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,8)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, std_conc_ESACCI_SSMI_SH,
             latlon=True, cmap='Reds')
plt.clim(0, 50)
plt.title("ESACCI_SSMI")

#colorbars and title
cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image, cax, orientation='vertical')
plt.suptitle("Sea-Ice Concentration - Standard Error (1992-2007) - Southern Hemisphere")
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Std_sh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Std_sh.pdf', format='pdf')
plt.close(fig)
