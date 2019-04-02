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
from scipy.stats import pearsonr
from netCDF4 import Dataset

#Get Lat/Lon to plot since this is regridded 1*1 (you can use any interp dataset for that) data else you define for each of the obs dataset as is done with plotting average concentration without regridding i.e on native grid 
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
lon, lat = np.meshgrid(lon, lat)
del data

#Dataset to compare/ reference data (1992-2008) and that is years we choose in all datasets
#ESACCI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/interp/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc') 
conc = data.variables['siconc'][:]
del data
average_conc_ESACCI_SSMI_NH = np.mean(conc, axis=0)
corr_coef_ESACCI_SSMI_NH = conc

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,9)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, average_conc_ESACCI_SSMI_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 100)
plt.title("ESACCI_SSMI")

#NSIDC0051 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc') #lat 448 lon 304 (25*25 resolution)
conc = data.variables['siconc'][156:360]
del data
corr_coef_NSIDC0051_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_NSIDC0051_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_NSIDC0051_NH[corr_coef_NSIDC0051_NH==0.] = np.nan # when your variable = 0, you set it to NaN
corr_coef_NSIDC0051_NH = np.ma.array(corr_coef_NSIDC0051_NH,mask=np.isnan(corr_coef_NSIDC0051_NH)) # you create a mask where all NaN values are taken apart from the rest (allow to make them white in a map)

fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,1)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image=m.pcolormesh(lon, lat, corr_coef_NSIDC0051_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("NSIDC0051")

#OSISAF409a - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][156:360]
del data
corr_coef_OSISAF409a_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_OSISAF409a_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_OSISAF409a_NH[corr_coef_OSISAF409a_NH==0.] = np.nan
corr_coef_OSISAF409a_NH = np.ma.array(corr_coef_OSISAF409a_NH,mask=np.isnan(corr_coef_OSISAF409a_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,2)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_OSISAF409a_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("OSISAF409a")

#OSISAF450 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][24:228]
del data
corr_coef_OSISAF450_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_OSISAF450_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_OSISAF450_NH[corr_coef_OSISAF450_NH==0.] = np.nan
corr_coef_OSISAF450_NH = np.ma.array(corr_coef_OSISAF450_NH,mask=np.isnan(corr_coef_OSISAF450_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,3)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_OSISAF450_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("OSISAF450")

#NSIDC_G02202_v3 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc')
conc = data.variables['siconc'][36:240]
del data
corr_coef_NSIDC_G02202_v3_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_NSIDC_G02202_v3_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_NSIDC_G02202_v3_NH[corr_coef_NSIDC_G02202_v3_NH==0.] = np.nan
corr_coef_NSIDC_G02202_v3_NH = np.ma.array(corr_coef_NSIDC_G02202_v3_NH,mask=np.isnan(corr_coef_NSIDC_G02202_v3_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,4)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_NSIDC_G02202_v3_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("NSIDC_G02202_v3")

#ICDC_ASI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc')
conc = data.variables['siconc'][0:204]
del data
corr_coef_ICDC_ASI_SSMI_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_ICDC_ASI_SSMI_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_ICDC_ASI_SSMI_NH[corr_coef_ICDC_ASI_SSMI_NH==0.] = np.nan
corr_coef_ICDC_ASI_SSMI_NH = np.ma.array(corr_coef_ICDC_ASI_SSMI_NH,mask=np.isnan(corr_coef_ICDC_ASI_SSMI_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,5)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_ICDC_ASI_SSMI_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("ICDC_ASI_SSMI")

#Hadisst1 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc = data.variables['siconc'][1464:1668]
del data
corr_coef_Hadisst1_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_Hadisst1_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_Hadisst1_NH[corr_coef_Hadisst1_NH==0.] = np.nan
corr_coef_Hadisst1_NH = np.ma.array(corr_coef_Hadisst1_NH,mask=np.isnan(corr_coef_Hadisst1_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,6)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_Hadisst1_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("Hadisst1")

#Hadisst2 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/interp/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc = data.variables['siconc'][1704:1908]
del data
corr_coef_Hadisst2_NH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_Hadisst2_NH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_NH[:,n,m])[0]
corr_coef_Hadisst2_NH[corr_coef_Hadisst2_NH==0.] = np.nan
corr_coef_Hadisst2_NH = np.ma.array(corr_coef_Hadisst2_NH,mask=np.isnan(corr_coef_Hadisst2_NH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,7)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_Hadisst2_NH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("Hadisst2")

#Plotting ESACCI_SSMI in the end for absolute values
plt.subplot(3,3,9)
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
#m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image2=m.pcolormesh(lon, lat, average_conc_ESACCI_SSMI_NH,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("ESACCI_SSMI")
cax2 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image2, cax2, orientation='horizontal')

#colorbars and title
cax = fig.add_axes([0.93, 0.2, 0.02, 0.6])
plt.colorbar(image, cax, orientation='vertical')
plt.suptitle("Sea-Ice Concentration - Correlation Coefficient with ESACCI_SSMI (1992-2007) - Northern Hemisphere")
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_CorrCoef_nh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_CorrCoef_nh.pdf', format='pdf')
plt.close(fig)

######################################################################################################################################################################################################################################################################
#Dataset to compare/ reference data (1992-2008)
#ESACCI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/interp/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc') 
conc = data.variables['siconc'][:]
del data
average_conc_ESACCI_SSMI_SH = np.mean(conc, axis=0)
corr_coef_ESACCI_SSMI_SH = conc

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,9)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, average_conc_ESACCI_SSMI_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 100)
plt.title("ESACCI_SSMI")

#NSIDC0051 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc') #lat 448 lon 304 (25*25 resolution)
conc = data.variables['siconc'][156:360]
del data
corr_coef_NSIDC0051_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_NSIDC0051_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_NSIDC0051_SH[corr_coef_NSIDC0051_SH==0.] = np.nan
corr_coef_NSIDC0051_SH = np.ma.array(corr_coef_NSIDC0051_SH,mask=np.isnan(corr_coef_NSIDC0051_SH))

fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,1)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image=m.pcolormesh(lon, lat, corr_coef_NSIDC0051_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("NSIDC0051")

#OSISAF409a - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][156:360]
del data
corr_coef_OSISAF409a_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_OSISAF409a_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_OSISAF409a_SH[corr_coef_OSISAF409a_SH==0.] = np.nan
corr_coef_OSISAF409a_SH = np.ma.array(corr_coef_OSISAF409a_SH,mask=np.isnan(corr_coef_OSISAF409a_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,2)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for Southh
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_OSISAF409a_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("OSISAF409a")

#OSISAF450 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc = data.variables['siconc'][24:228]
del data
corr_coef_OSISAF450_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_OSISAF450_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_OSISAF450_SH[corr_coef_OSISAF450_SH==0.] = np.nan
corr_coef_OSISAF450_SH = np.ma.array(corr_coef_OSISAF450_SH,mask=np.isnan(corr_coef_OSISAF450_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,3)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_OSISAF450_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("OSISAF450")

#NSIDC_G02202_v3 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc')
conc = data.variables['siconc'][36:240]
del data
corr_coef_NSIDC_G02202_v3_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_NSIDC_G02202_v3_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_NSIDC_G02202_v3_SH[corr_coef_NSIDC_G02202_v3_SH==0.] = np.nan
corr_coef_NSIDC_G02202_v3_SH = np.ma.array(corr_coef_NSIDC_G02202_v3_SH,mask=np.isnan(corr_coef_NSIDC_G02202_v3_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,4)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_NSIDC_G02202_v3_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("NSIDC_G02202_v3")

#ICDC_ASI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc')
conc = data.variables['siconc'][0:204]
del data
corr_coef_ICDC_ASI_SSMI_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_ICDC_ASI_SSMI_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_ICDC_ASI_SSMI_SH[corr_coef_ICDC_ASI_SSMI_SH==0.] = np.nan
corr_coef_ICDC_ASI_SSMI_SH = np.ma.array(corr_coef_ICDC_ASI_SSMI_SH,mask=np.isnan(corr_coef_ICDC_ASI_SSMI_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,5)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_ICDC_ASI_SSMI_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("ICDC_ASI_SSMI")

#Hadisst1 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/interp/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc = data.variables['siconc'][1464:1668]
del data
corr_coef_Hadisst1_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_Hadisst1_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_Hadisst1_SH[corr_coef_Hadisst1_SH==0.] = np.nan
corr_coef_Hadisst1_SH = np.ma.array(corr_coef_Hadisst1_SH,mask=np.isnan(corr_coef_Hadisst1_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,6)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_Hadisst1_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("Hadisst1")

#Hadisst2 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/interp/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc = data.variables['siconc'][1704:1908]
del data
corr_coef_Hadisst2_SH = np.zeros((conc.shape[1], conc.shape[2]))
for n in range(conc.shape[1]):
   for m in range(conc.shape[2]):
      corr_coef_Hadisst2_SH[n, m] = pearsonr(conc[:,n,m], corr_coef_ESACCI_SSMI_SH[:,n,m])[0]
corr_coef_Hadisst2_SH[corr_coef_Hadisst2_SH==0.] = np.nan
corr_coef_Hadisst2_SH = np.ma.array(corr_coef_Hadisst2_SH,mask=np.isnan(corr_coef_Hadisst2_SH))

#fig = plt.figure(figsize=(10, 8))
plt.subplot(3,3,7)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
m.pcolormesh(lon, lat, corr_coef_Hadisst2_SH,
             latlon=True, cmap='RdBu_r')
plt.clim(0, 1)
plt.title("Hadisst2")

#Plotting ESACCI_SSMI in the end for absolute values
plt.subplot(3,3,9)
#m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # for North
m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l') # for South
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0] )
m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,1,0,0] )
image2=m.pcolormesh(lon, lat, average_conc_ESACCI_SSMI_SH,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("ESACCI_SSMI")
cax2 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image2, cax2, orientation='horizontal')

#colorbars and title
cax = fig.add_axes([0.93, 0.2, 0.02, 0.6])
plt.colorbar(image, cax, orientation='vertical')
plt.suptitle("Sea-Ice Concentration - Correlation Coefficient with ESACCI_SSMI (1992-2007) - Southern Hemisphere")
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_CorrCoef_sh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_CorrCoef_sh.pdf', format='pdf')
plt.close(fig)
