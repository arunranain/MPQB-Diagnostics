#!/usr/bin/python
# Common sea ice diags
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
from seaice_commondiags import compute_extent
from seaice_commondiags import compute_area
from netCDF4 import Dataset

#NSIDC0051 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc') #lat 448 lon 304 (25*25 resolution)
conc_NSIDC0051_NH = data['siconc'][:]
area_NSIDC0051_NH = data['areacello'][:]/100
mask_NSIDC0051_NH = data['sftof'][:]
ext_NSIDC0051_NH = compute_extent(conc_NSIDC0051_NH, area_NSIDC0051_NH, threshold = 15.0, mask = mask_NSIDC0051_NH)
are_NSIDC0051_NH = compute_area(conc_NSIDC0051_NH, area_NSIDC0051_NH, mask = mask_NSIDC0051_NH)
sea_NSIDC0051_NH = compute_extent(conc_NSIDC0051_NH, area_NSIDC0051_NH, threshold = -1.0, mask = mask_NSIDC0051_NH)

#OSISAF409a - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF409a_NH = data['siconc'][:]
area_OSISAF409a_NH = data['areacello'][0,:,:]
mask_OSISAF409a_NH = data['sftof'][0,:,:]
ext_OSISAF409a_NH = compute_extent(conc_OSISAF409a_NH, area_OSISAF409a_NH, threshold = 15.0, mask = mask_OSISAF409a_NH) 
are_OSISAF409a_NH = compute_area(conc_OSISAF409a_NH, area_OSISAF409a_NH, mask = mask_OSISAF409a_NH)
sea_OSISAF409a_NH = compute_extent(conc_OSISAF409a_NH, area_OSISAF409a_NH, threshold = -1.0, mask = mask_OSISAF409a_NH)

#OSISAF450 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF450_NH = data['siconc'][:]
area_OSISAF450_NH = data['areacello'][0,:,:]
mask_OSISAF450_NH = data['sftof'][0,:,:]
ext_OSISAF450_NH = compute_extent(conc_OSISAF450_NH, area_OSISAF450_NH, threshold = 15.0, mask = mask_OSISAF450_NH) 
are_OSISAF450_NH = compute_area(conc_OSISAF450_NH, area_OSISAF450_NH, mask = mask_OSISAF450_NH)
sea_OSISAF450_NH = compute_extent(conc_OSISAF450_NH, area_OSISAF450_NH, threshold = -1.0, mask = mask_OSISAF450_NH)

#NSIDC_G02202_v3 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc')
conc_NSIDC_G02202_v3_NH = data['siconc'][:]
area_NSIDC_G02202_v3_NH = data['areacello'][0,:,:]
mask_NSIDC_G02202_v3_NH = data['sftof'][0,:,:]
ext_NSIDC_G02202_v3_NH = compute_extent(conc_NSIDC_G02202_v3_NH, area_NSIDC_G02202_v3_NH, threshold = 15.0, mask = mask_NSIDC_G02202_v3_NH) 
are_NSIDC_G02202_v3_NH = compute_area(conc_NSIDC_G02202_v3_NH, area_NSIDC_G02202_v3_NH, mask = mask_NSIDC_G02202_v3_NH)
sea_NSIDC_G02202_v3_NH = compute_extent(conc_NSIDC_G02202_v3_NH, area_NSIDC_G02202_v3_NH, threshold = -1.0, mask = mask_NSIDC_G02202_v3_NH)

#ESACCI_AMSR - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/AMSR/processed/native/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_NH.nc')
conc_ESACCI_AMSR_NH = data['siconc'][:]
area_ESACCI_AMSR_NH = data['areacello'][0,:,:]
mask_ESACCI_AMSR_NH = data['sftof'][0,:,:]
ext_ESACCI_AMSR_NH = compute_extent(conc_ESACCI_AMSR_NH, area_ESACCI_AMSR_NH, threshold = 15.0, mask = mask_ESACCI_AMSR_NH) 
are_ESACCI_AMSR_NH = compute_area(conc_ESACCI_AMSR_NH, area_ESACCI_AMSR_NH, mask = mask_ESACCI_AMSR_NH)
sea_ESACCI_AMSR_NH = compute_extent(conc_ESACCI_AMSR_NH, area_ESACCI_AMSR_NH, threshold = -1.0, mask = mask_ESACCI_AMSR_NH)

#ICDC_ASI_AMSRE - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-AMSRE/processed/native/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_n.nc')
conc_ICDC_ASI_AMSRE_NH = data['siconc'][:]
area_ICDC_ASI_AMSRE_NH = data['areacello'][0,:,:]
mask_ICDC_ASI_AMSRE_NH = data['sftof'][0,:,:]
ext_ICDC_ASI_AMSRE_NH = compute_extent(conc_ICDC_ASI_AMSRE_NH, area_ICDC_ASI_AMSRE_NH, threshold = 15.0, mask = mask_ICDC_ASI_AMSRE_NH) 
are_ICDC_ASI_AMSRE_NH = compute_area(conc_ICDC_ASI_AMSRE_NH, area_ICDC_ASI_AMSRE_NH, mask = mask_ICDC_ASI_AMSRE_NH)
sea_ICDC_ASI_AMSRE_NH = compute_extent(conc_ICDC_ASI_AMSRE_NH, area_ICDC_ASI_AMSRE_NH, threshold = -1.0, mask = mask_ICDC_ASI_AMSRE_NH)

#ICDC_ASI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc')
conc_ICDC_ASI_SSMI_NH = data['siconc'][:]
area_ICDC_ASI_SSMI_NH = data['areacello'][0,:,:]
mask_ICDC_ASI_SSMI_NH = data['sftof'][0,:,:]
ext_ICDC_ASI_SSMI_NH = compute_extent(conc_ICDC_ASI_SSMI_NH, area_ICDC_ASI_SSMI_NH, threshold = 15.0, mask = mask_ICDC_ASI_SSMI_NH) 
are_ICDC_ASI_SSMI_NH = compute_area(conc_ICDC_ASI_SSMI_NH, area_ICDC_ASI_SSMI_NH, mask = mask_ICDC_ASI_SSMI_NH)
sea_ICDC_ASI_SSMI_NH = compute_extent(conc_ICDC_ASI_SSMI_NH, area_ICDC_ASI_SSMI_NH, threshold = -1.0, mask = mask_ICDC_ASI_SSMI_NH)

#Hadisst1 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc_Hadisst1_NH = data['siconc'][:]
area_Hadisst1_NH = data['areacello'][:]
mask_Hadisst1_NH = data['sftof'][0,:,:]
#### selecting just northern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/HadISST_lat.nc')
lat = flat["latitude"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/HadISST_lon.nc')
lon = flon["longitude"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] > 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_Hadisst1_NH = mask_Hadisst1_NH*a
ext_Hadisst1_NH = compute_extent(conc_Hadisst1_NH, area_Hadisst1_NH, threshold = 15.0, mask = mask_Hadisst1_NH) 
are_Hadisst1_NH = compute_area(conc_Hadisst1_NH, area_Hadisst1_NH, mask = mask_Hadisst1_NH)
sea_Hadisst1_NH = compute_extent(conc_Hadisst1_NH, area_Hadisst1_NH, threshold = -1.0, mask = mask_Hadisst1_NH)

#Hadisst2 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/native/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc_Hadisst2_NH = data['siconc'][:]
area_Hadisst2_NH = data['areacello'][:]
mask_Hadisst2_NH = data['sftof'][0,:,:]
#### selecting just northern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/HadISST2_lat.nc')
lat = flat["latitude"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/HadISST2_lon.nc')
lon = flon["longitude"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] > 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_Hadisst2_NH = mask_Hadisst2_NH*a
ext_Hadisst2_NH = compute_extent(conc_Hadisst2_NH, area_Hadisst2_NH, threshold = 15.0, mask = mask_Hadisst2_NH) 
are_Hadisst2_NH = compute_area(conc_Hadisst2_NH, area_Hadisst2_NH, mask = mask_Hadisst2_NH)
sea_Hadisst2_NH = compute_extent(conc_Hadisst2_NH, area_Hadisst2_NH, threshold = -1.0, mask = mask_Hadisst2_NH)

#NOAA_OI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NOAA/NOAA-OI/processed/native/siconc_SImon_NOAA_OI_r1i1p1_198912-201809_reg1.0.nc')
conc_NOAA_OI_NH = data['siconc'][:]
area_NOAA_OI_NH = data['areacello'][:]
mask_NOAA_OI_NH = data['sftof'][0,:,:]
#### selecting just northern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/NOAA_OI_lat.nc')
lat = flat["lat"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/NOAA_OI_lon.nc')
lon = flon["lon"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] > 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_NOAA_OI_NH = mask_NOAA_OI_NH*a
ext_NOAA_OI_NH = compute_extent(conc_NOAA_OI_NH, area_NOAA_OI_NH, threshold = 15.0, mask = mask_NOAA_OI_NH) 
are_NOAA_OI_NH = compute_area(conc_NOAA_OI_NH, area_NOAA_OI_NH, mask = mask_NOAA_OI_NH)
sea_NOAA_OI_NH = compute_extent(conc_NOAA_OI_NH, area_NOAA_OI_NH, threshold = -1.0, mask = mask_NOAA_OI_NH)

#ESACCI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc') 
conc_ESACCI_SSMI_NH = data['siconc'][:]
area_ESACCI_SSMI_NH = data['areacello'][0,:,:]
mask_ESACCI_SSMI_NH = data['sftof'][0,:,:]
ext_ESACCI_SSMI_NH = compute_extent(conc_ESACCI_SSMI_NH, area_ESACCI_SSMI_NH, threshold = 15.0, mask = mask_ESACCI_SSMI_NH) 
are_ESACCI_SSMI_NH = compute_area(conc_ESACCI_SSMI_NH, area_ESACCI_SSMI_NH, mask = mask_ESACCI_SSMI_NH)
sea_ESACCI_SSMI_NH = compute_extent(conc_ESACCI_SSMI_NH, area_ESACCI_SSMI_NH, threshold = -1.0, mask = mask_ESACCI_SSMI_NH)

# create time data to plot
dates = []
for years in range(1977, 2017):
    for months in range(1, 13):
        dates.append(dt.datetime(year=years, month=months, day=1))

fig = plt.figure()
ax = plt.axes()
plt.plot(dates[0:468], ext_NSIDC0051_NH, label="NSIDC0051")   #start year 197901 end year 201712
plt.plot(dates[0:444], ext_OSISAF409a_NH, label="OSISAF409a")  #start year 197901 end year 201512
plt.plot(dates[132:444], ext_OSISAF450_NH, label="OSISAF450")  #start year 199001 end year 201512
plt.plot(dates[120:468], ext_NSIDC_G02202_v3_NH, label="NSIDC_G02202_v3")  #start year 198901 end year 201712
plt.plot(dates[288:384], ext_ESACCI_AMSR_NH, label="ESACCI_AMSR")  #start year 200301 end year 201012
plt.plot(dates[156:360], ext_ESACCI_SSMI_NH, label="ESACCI_SSMI")  #start year 199201 end year 200812
plt.plot(dates[288:396], ext_ICDC_ASI_AMSRE_NH, label="ICDC_ASI_AMSRE")  #start year 200301 end year 201112
plt.plot(dates[156:480], ext_ICDC_ASI_SSMI_NH, label="ICDC_ASI_SSMI")  #start year 199201 end year 201812
#plt.plot(dates[120:466], ext_NOAA_OI_NH, label="NOAA_OI")  #start year 198901 end year 201809 -- Here monthly data is averaged from weekly values and thus the time between 1989 to 2018 is 346 steps and not 357
plt.plot(dates[0:468], ext_Hadisst1_NH[1308:1776], label="Hadisst1")  #start year 187001 end year 201805 but we plot from 1979 onwards till 201712
plt.plot(dates[0:468], ext_Hadisst2_NH[1548:2016], label="Hadisst2")  #start year 185001 end year 201803 but we plot from 1979 onwards till 201712

# multiple line plot in various subplots

#plt.figure(1)
#plt.subplot(211)
#plt.plot(t, s1)
#plt.subplot(212)
#plt.plot(t, 2*s1)

#or in same plot window

#plt.plot( 'x', 'y1', data=df, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
#plt.plot( 'x', 'y2', data=df, marker='', color='olive', linewidth=2)
#plt.plot( 'x', 'y3', data=df, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")

# format the ticks
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Northern Hemisphere")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh.pdf', format='pdf')
plt.close(fig)

# create month data to plot - September
fig = plt.figure()
ax = plt.axes()

ext_NSIDC0051_NH_SEP = ext_NSIDC0051_NH [8::12]
ext_OSISAF409a_NH_SEP = ext_OSISAF409a_NH [8::12]
ext_OSISAF450_NH_SEP = ext_OSISAF450_NH [8::12]
ext_NSIDC_G02202_v3_NH_SEP = ext_NSIDC_G02202_v3_NH [8::12]
ext_ESACCI_AMSR_NH_SEP = ext_ESACCI_AMSR_NH [8::12]
ext_ESACCI_SSMI_NH_SEP = ext_ESACCI_SSMI_NH [8::12]
ext_ICDC_ASI_AMSRE_NH_SEP = ext_ICDC_ASI_AMSRE_NH [8::12]
ext_ICDC_ASI_SSMI_NH_SEP = ext_ICDC_ASI_SSMI_NH [8::12]
ext_NOAA_OI_NH_SEP = ext_NOAA_OI_NH [8::12]
ext_Hadisst1_NH_SEP = ext_Hadisst1_NH [8::12]
ext_Hadisst2_NH_SEP = ext_Hadisst2_NH [8::12]
dates_SEP = dates [8::12]

plt.plot(dates_SEP[0:39], ext_NSIDC0051_NH_SEP, label="NSIDC0051")   
plt.plot(dates_SEP[0:37], ext_OSISAF409a_NH_SEP, label="OSISAF409a") 
plt.plot(dates_SEP[11:37], ext_OSISAF450_NH_SEP, label="OSISAF450") 
plt.plot(dates_SEP[10:39], ext_NSIDC_G02202_v3_NH_SEP, label="NSIDC_G02202_v3")  
plt.plot(dates_SEP[24:32], ext_ESACCI_AMSR_NH_SEP, label="ESACCI_AMSR") 
plt.plot(dates_SEP[13:30], ext_ESACCI_SSMI_NH_SEP, label="ESACCI_SSMI")  
plt.plot(dates_SEP[24:33], ext_ICDC_ASI_AMSRE_NH_SEP, label="ICDC_ASI_AMSRE")  
plt.plot(dates_SEP[13:40], ext_ICDC_ASI_SSMI_NH_SEP, label="ICDC_ASI_SSMI")  
#plt.plot(dates_SEP[10:39], ext_NOAA_OI_NH_SEP, label="NOAA_OI")  
plt.plot(dates_SEP[0:39], ext_Hadisst1_NH_SEP[109:148], label="Hadisst1") 
plt.plot(dates_SEP[0:39], ext_Hadisst2_NH_SEP[129:168], label="Hadisst2") 

# format the ticks
years = mdates.YearLocator()   # every year
#months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Northern Hemisphere -September")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh_SEP.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh_SEP.pdf', format='pdf')
plt.close(fig)

# create month data to plot - March
fig = plt.figure()
ax = plt.axes()

ext_NSIDC0051_NH_MAR = ext_NSIDC0051_NH [2::12]
ext_OSISAF409a_NH_MAR = ext_OSISAF409a_NH [2::12]
ext_OSISAF450_NH_MAR = ext_OSISAF450_NH [2::12]
ext_NSIDC_G02202_v3_NH_MAR = ext_NSIDC_G02202_v3_NH [2::12]
ext_ESACCI_AMSR_NH_MAR = ext_ESACCI_AMSR_NH [2::12]
ext_ESACCI_SSMI_NH_MAR = ext_ESACCI_SSMI_NH [2::12]
ext_ICDC_ASI_AMSRE_NH_MAR = ext_ICDC_ASI_AMSRE_NH [2::12]
ext_ICDC_ASI_SSMI_NH_MAR = ext_ICDC_ASI_SSMI_NH [2::12]
ext_NOAA_OI_NH_MAR = ext_NOAA_OI_NH [2::12]
ext_Hadisst1_NH_MAR = ext_Hadisst1_NH [2::12]
ext_Hadisst2_NH_MAR = ext_Hadisst2_NH [2::12]
dates_MAR = dates [2::12]

plt.plot(dates_MAR[0:39], ext_NSIDC0051_NH_MAR, label="NSIDC0051")   
plt.plot(dates_MAR[0:37], ext_OSISAF409a_NH_MAR, label="OSISAF409a") 
plt.plot(dates_MAR[11:37], ext_OSISAF450_NH_MAR, label="OSISAF450") 
plt.plot(dates_MAR[10:39], ext_NSIDC_G02202_v3_NH_MAR, label="NSIDC_G02202_v3")  
plt.plot(dates_MAR[24:32], ext_ESACCI_AMSR_NH_MAR, label="ESACCI_AMSR") 
plt.plot(dates_MAR[13:30], ext_ESACCI_SSMI_NH_MAR, label="ESACCI_SSMI")  
plt.plot(dates_MAR[24:33], ext_ICDC_ASI_AMSRE_NH_MAR, label="ICDC_ASI_AMSRE")  
plt.plot(dates_MAR[13:40], ext_ICDC_ASI_SSMI_NH_MAR, label="ICDC_ASI_SSMI")  
#plt.plot(dates_MAR[10:39], ext_NOAA_OI_NH_MAR, label="NOAA_OI")  
plt.plot(dates_MAR[0:39], ext_Hadisst1_NH_MAR[109:148], label="Hadisst1") 
plt.plot(dates_MAR[0:39], ext_Hadisst2_NH_MAR[129:168], label="Hadisst2") 

# format the ticks
years = mdates.YearLocator()   # every year
#months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Northern Hemisphere - March")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh_MAR.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_nh_MAR.pdf', format='pdf')
plt.close(fig)

######################################################################################################################################################################################################################################################################
#NSIDC0051 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc') #lat 448 lon 304 (25*25 resolution)
conc_NSIDC0051_SH = data['siconc'][:]
area_NSIDC0051_SH = data['areacello'][:]/100
mask_NSIDC0051_SH = data['sftof'][:]
ext_NSIDC0051_SH = compute_extent(conc_NSIDC0051_SH, area_NSIDC0051_SH, threshold = 15.0, mask = mask_NSIDC0051_SH)
are_NSIDC0051_SH = compute_area(conc_NSIDC0051_SH, area_NSIDC0051_SH, mask = mask_NSIDC0051_SH)
sea_NSIDC0051_SH = compute_extent(conc_NSIDC0051_SH, area_NSIDC0051_SH, threshold = -1.0, mask = mask_NSIDC0051_SH)

#OSISAF409a - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF409a_SH = data['siconc'][:]
area_OSISAF409a_SH = data['areacello'][0,:,:]
mask_OSISAF409a_SH = data['sftof'][0,:,:]
ext_OSISAF409a_SH = compute_extent(conc_OSISAF409a_SH, area_OSISAF409a_SH, threshold = 15.0, mask = mask_OSISAF409a_SH) 
are_OSISAF409a_SH = compute_area(conc_OSISAF409a_SH, area_OSISAF409a_SH, mask = mask_OSISAF409a_SH)
sea_OSISAF409a_SH = compute_extent(conc_OSISAF409a_SH, area_OSISAF409a_SH, threshold = -1.0, mask = mask_OSISAF409a_SH)

#OSISAF450 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF450_SH = data['siconc'][:]
area_OSISAF450_SH = data['areacello'][0,:,:]
mask_OSISAF450_SH = data['sftof'][0,:,:]
ext_OSISAF450_SH = compute_extent(conc_OSISAF450_SH, area_OSISAF450_SH, threshold = 15.0, mask = mask_OSISAF450_SH) 
are_OSISAF450_SH = compute_area(conc_OSISAF450_SH, area_OSISAF450_SH, mask = mask_OSISAF450_SH)
sea_OSISAF450_SH = compute_extent(conc_OSISAF450_SH, area_OSISAF450_SH, threshold = -1.0, mask = mask_OSISAF450_SH)

#NSIDC_G02202_v3 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc')
conc_NSIDC_G02202_v3_SH = data['siconc'][:]
area_NSIDC_G02202_v3_SH = data['areacello'][0,:,:]
mask_NSIDC_G02202_v3_SH = data['sftof'][0,:,:]
ext_NSIDC_G02202_v3_SH = compute_extent(conc_NSIDC_G02202_v3_SH, area_NSIDC_G02202_v3_SH, threshold = 15.0, mask = mask_NSIDC_G02202_v3_SH) 
are_NSIDC_G02202_v3_SH = compute_area(conc_NSIDC_G02202_v3_SH, area_NSIDC_G02202_v3_SH, mask = mask_NSIDC_G02202_v3_SH)
sea_NSIDC_G02202_v3_SH = compute_extent(conc_NSIDC_G02202_v3_SH, area_NSIDC_G02202_v3_SH, threshold = -1.0, mask = mask_NSIDC_G02202_v3_SH)

#ESACCI_AMSR - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/AMSR/processed/native/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_SH.nc')
conc_ESACCI_AMSR_SH = data['siconc'][:]
area_ESACCI_AMSR_SH = data['areacello'][0,:,:]
mask_ESACCI_AMSR_SH = data['sftof'][0,:,:]
ext_ESACCI_AMSR_SH = compute_extent(conc_ESACCI_AMSR_SH, area_ESACCI_AMSR_SH, threshold = 15.0, mask = mask_ESACCI_AMSR_SH) 
are_ESACCI_AMSR_SH = compute_area(conc_ESACCI_AMSR_SH, area_ESACCI_AMSR_SH, mask = mask_ESACCI_AMSR_SH)
sea_ESACCI_AMSR_SH = compute_extent(conc_ESACCI_AMSR_SH, area_ESACCI_AMSR_SH, threshold = -1.0, mask = mask_ESACCI_AMSR_SH) 

#ICDC_ASI_AMSRE - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-AMSRE/processed/native/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_s.nc')
conc_ICDC_ASI_AMSRE_SH = data['siconc'][:]
area_ICDC_ASI_AMSRE_SH = data['areacello'][0,:,:]
mask_ICDC_ASI_AMSRE_SH = data['sftof'][0,:,:]
ext_ICDC_ASI_AMSRE_SH = compute_extent(conc_ICDC_ASI_AMSRE_SH, area_ICDC_ASI_AMSRE_SH, threshold = 15.0, mask = mask_ICDC_ASI_AMSRE_SH) 
are_ICDC_ASI_AMSRE_SH = compute_area(conc_ICDC_ASI_AMSRE_SH, area_ICDC_ASI_AMSRE_SH, mask = mask_ICDC_ASI_AMSRE_SH)
sea_ICDC_ASI_AMSRE_SH = compute_extent(conc_ICDC_ASI_AMSRE_SH, area_ICDC_ASI_AMSRE_SH, threshold = -1.0, mask = mask_ICDC_ASI_AMSRE_SH) 

#ICDC_ASI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc')
conc_ICDC_ASI_SSMI_SH = data['siconc'][:]
area_ICDC_ASI_SSMI_SH = data['areacello'][0,:,:]
mask_ICDC_ASI_SSMI_SH = data['sftof'][0,:,:]
ext_ICDC_ASI_SSMI_SH = compute_extent(conc_ICDC_ASI_SSMI_SH, area_ICDC_ASI_SSMI_SH, threshold = 15.0, mask = mask_ICDC_ASI_SSMI_SH) 
are_ICDC_ASI_SSMI_SH = compute_area(conc_ICDC_ASI_SSMI_SH, area_ICDC_ASI_SSMI_SH, mask = mask_ICDC_ASI_SSMI_SH)
sea_ICDC_ASI_SSMI_SH = compute_extent(conc_ICDC_ASI_SSMI_SH, area_ICDC_ASI_SSMI_SH, threshold = -1.0, mask = mask_ICDC_ASI_SSMI_SH) 

#Hadisst1 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc_Hadisst1_SH = data['siconc'][:]
area_Hadisst1_SH = data['areacello'][:]
mask_Hadisst1_SH = data['sftof'][0,:,:]
#### selecting just southern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/HadISST_lat.nc')
lat = flat["latitude"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/HadISST_lon.nc')
lon = flon["longitude"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] < 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_Hadisst1_SH = mask_Hadisst1_SH*a
ext_Hadisst1_SH = compute_extent(conc_Hadisst1_SH, area_Hadisst1_SH, threshold = 15.0, mask = mask_Hadisst1_SH) 
are_Hadisst1_SH = compute_area(conc_Hadisst1_SH, area_Hadisst1_SH, mask = mask_Hadisst1_SH)
sea_Hadisst1_SH = compute_extent(conc_Hadisst1_SH, area_Hadisst1_SH, threshold = -1.0, mask = mask_Hadisst1_SH)

#Hadisst2 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/native/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc_Hadisst2_SH = data['siconc'][:]
area_Hadisst2_SH = data['areacello'][:]
mask_Hadisst2_SH = data['sftof'][0,:,:]
#### selecting just southern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/HadISST2_lat.nc')
lat = flat["latitude"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/HadISST2_lon.nc')
lon = flon["longitude"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] < 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_Hadisst2_SH = mask_Hadisst2_SH*a
ext_Hadisst2_SH = compute_extent(conc_Hadisst2_SH, area_Hadisst2_SH, threshold = 15.0, mask = mask_Hadisst2_SH) 
are_Hadisst2_SH = compute_area(conc_Hadisst2_SH, area_Hadisst2_SH, mask = mask_Hadisst2_SH)
sea_Hadisst2_SH = compute_extent(conc_Hadisst2_SH, area_Hadisst2_SH, threshold = -1.0, mask = mask_Hadisst2_SH) 

#NOAA_OI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NOAA/NOAA-OI/processed/native/siconc_SImon_NOAA_OI_r1i1p1_198912-201809_reg1.0.nc')
conc_NOAA_OI_SH = data['siconc'][:]
area_NOAA_OI_SH = data['areacello'][:]
mask_NOAA_OI_SH = data['sftof'][0,:,:]
#### selecting just southern hemisphere
flat = Dataset('/storepelican/CLIMDATA/tools/NOAA_OI_lat.nc')
lat = flat["lat"][:]
flon = Dataset('/storepelican/CLIMDATA/tools/NOAA_OI_lon.nc')
lon = flon["lon"][:]
x, y = np.meshgrid(lon, lat, indexing='ij')
coordinate_grid = np.array([x, y])
a=np.zeros(shape=(len(lat),len(lon)))
for i in range(len(coordinate_grid[1][1][:])):
   if coordinate_grid[1][1][i] < 0:
      a[:][i] = 1
   else:
      a[:][i] = 0
mask_NOAA_OI_SH = mask_NOAA_OI_SH*a
ext_NOAA_OI_SH = compute_extent(conc_NOAA_OI_SH, area_NOAA_OI_SH, threshold = 15.0, mask = mask_NOAA_OI_SH) 
are_NOAA_OI_SH = compute_area(conc_NOAA_OI_SH, area_NOAA_OI_SH, mask = mask_NOAA_OI_SH)
sea_NOAA_OI_SH = compute_extent(conc_NOAA_OI_SH, area_NOAA_OI_SH, threshold = -1.0, mask = mask_NOAA_OI_SH) 

#ESACCI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc') 
conc_ESACCI_SSMI_SH = data['siconc'][:]
area_ESACCI_SSMI_SH = data['areacello'][0,:,:]
mask_ESACCI_SSMI_SH = data['sftof'][0,:,:]
ext_ESACCI_SSMI_SH = compute_extent(conc_ESACCI_SSMI_SH, area_ESACCI_SSMI_SH, threshold = 15.0, mask = mask_ESACCI_SSMI_SH) 
are_ESACCI_SSMI_SH = compute_area(conc_ESACCI_SSMI_SH, area_ESACCI_SSMI_SH, mask = mask_ESACCI_SSMI_SH)
sea_ESACCI_SSMI_SH = compute_extent(conc_ESACCI_SSMI_SH, area_ESACCI_SSMI_SH, threshold = -1.0, mask = mask_ESACCI_SSMI_SH)

# create time data to plot
dates = []
for years in range(1977, 2017):
    for months in range(1, 13):
        dates.append(dt.datetime(year=years, month=months, day=1))

fig = plt.figure()
ax = plt.axes()
plt.plot(dates[0:468], ext_NSIDC0051_SH, label="NSIDC0051")   #start year 197901 end year 201712
plt.plot(dates[0:444], ext_OSISAF409a_SH, label="OSISAF409a")  #start year 197901 end year 201512
plt.plot(dates[132:444], ext_OSISAF450_SH, label="OSISAF450")  #start year 197901 end year 201512
plt.plot(dates[120:468], ext_NSIDC_G02202_v3_SH, label="NSIDC_G02202_v3")  #start year 198901 end year 201712
plt.plot(dates[288:384], ext_ESACCI_AMSR_SH, label="ESACCI_AMSR")  #start year 200301 end year 201012
plt.plot(dates[156:360], ext_ESACCI_SSMI_SH, label="ESACCI_SSMI")  #start year 199201 end year 200812
plt.plot(dates[288:396], ext_ICDC_ASI_AMSRE_SH, label="ICDC_ASI_AMSRE")  #start year 200301 end year 201112
plt.plot(dates[156:480], ext_ICDC_ASI_SSMI_SH, label="ICDC_ASI_SSMI")  #start year 199201 end year 201812
#plt.plot(dates[120:466], ext_NOAA_OI_NH, label="NOAA_OI")  #start year 198901 end year 201809 -- Here monthly data is averaged from weekly values and thus the time between 1989 to 2018 is 346 steps and not 357
plt.plot(dates[0:468], ext_Hadisst1_SH[1308:1776], label="Hadisst1")  #start year 187001 end year 201805 but we plot from 1979 onwards till 201712
plt.plot(dates[0:468], ext_Hadisst2_SH[1548:2016], label="Hadisst2")  #start year 185001 end year 201803 but we plot from 1979 onwards till 201712

# format the ticks
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Southern Hemisphere")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh.pdf', format='pdf')
plt.close(fig)

# create month data to plot - September
fig = plt.figure()
ax = plt.axes()

ext_NSIDC0051_SH_SEP = ext_NSIDC0051_SH [8::12]
ext_OSISAF409a_SH_SEP = ext_OSISAF409a_SH [8::12]
ext_OSISAF450_SH_SEP = ext_OSISAF450_SH [8::12]
ext_NSIDC_G02202_v3_SH_SEP = ext_NSIDC_G02202_v3_SH [8::12]
ext_ESACCI_AMSR_SH_SEP = ext_ESACCI_AMSR_SH [8::12]
ext_ESACCI_SSMI_SH_SEP = ext_ESACCI_SSMI_SH [8::12]
ext_ICDC_ASI_AMSRE_SH_SEP = ext_ICDC_ASI_AMSRE_SH [8::12]
ext_ICDC_ASI_SSMI_SH_SEP = ext_ICDC_ASI_SSMI_SH [8::12]
ext_NOAA_OI_SH_SEP = ext_NOAA_OI_SH [8::12]
ext_Hadisst1_SH_SEP = ext_Hadisst1_SH [8::12]
ext_Hadisst2_SH_SEP = ext_Hadisst2_SH [8::12]
dates_SEP = dates [8::12]

plt.plot(dates_SEP[0:39], ext_NSIDC0051_SH_SEP, label="NSIDC0051")   
plt.plot(dates_SEP[0:37], ext_OSISAF409a_SH_SEP, label="OSISAF409a") 
plt.plot(dates_SEP[11:37], ext_OSISAF450_SH_SEP, label="OSISAF450") 
plt.plot(dates_SEP[10:39], ext_NSIDC_G02202_v3_SH_SEP, label="NSIDC_G02202_v3")  
plt.plot(dates_SEP[24:32], ext_ESACCI_AMSR_SH_SEP, label="ESACCI_AMSR") 
plt.plot(dates_SEP[13:30], ext_ESACCI_SSMI_SH_SEP, label="ESACCI_SSMI")  
plt.plot(dates_SEP[24:33], ext_ICDC_ASI_AMSRE_SH_SEP, label="ICDC_ASI_AMSRE")  
plt.plot(dates_SEP[13:40], ext_ICDC_ASI_SSMI_SH_SEP, label="ICDC_ASI_SSMI")  
#plt.plot(dates_SEP[10:39], ext_NOAA_OI_SH_SEP, label="NOAA_OI")  
plt.plot(dates_SEP[0:39], ext_Hadisst1_SH_SEP[109:148], label="Hadisst1") 
plt.plot(dates_SEP[0:39], ext_Hadisst2_SH_SEP[129:168], label="Hadisst2") 

# format the ticks
years = mdates.YearLocator()   # every year
#months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Southern Hemisphere - September")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh_SEP.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh_SEP.pdf', format='pdf')
plt.close(fig)

# create month data to plot - February
fig = plt.figure()
ax = plt.axes()

ext_NSIDC0051_SH_FEB = ext_NSIDC0051_SH [1::12]
ext_OSISAF409a_SH_FEB = ext_OSISAF409a_SH [1::12]
ext_OSISAF450_SH_FEB = ext_OSISAF450_SH [1::12]
ext_NSIDC_G02202_v3_SH_FEB = ext_NSIDC_G02202_v3_SH [1::12]
ext_ESACCI_AMSR_SH_FEB = ext_ESACCI_AMSR_SH [1::12]
ext_ESACCI_SSMI_SH_FEB = ext_ESACCI_SSMI_SH [1::12]
ext_ICDC_ASI_AMSRE_SH_FEB = ext_ICDC_ASI_AMSRE_SH [1::12]
ext_ICDC_ASI_SSMI_SH_FEB = ext_ICDC_ASI_SSMI_SH [1::12]
ext_NOAA_OI_SH_FEB = ext_NOAA_OI_SH [1::12]
ext_Hadisst1_SH_FEB = ext_Hadisst1_SH [1::12]
ext_Hadisst2_SH_FEB = ext_Hadisst2_SH [1::12]
dates_FEB = dates [1::12]

plt.plot(dates_FEB[0:39], ext_NSIDC0051_SH_FEB, label="NSIDC0051")   
plt.plot(dates_FEB[0:37], ext_OSISAF409a_SH_FEB, label="OSISAF409a") 
plt.plot(dates_FEB[11:37], ext_OSISAF450_SH_FEB, label="OSISAF450") 
plt.plot(dates_FEB[10:39], ext_NSIDC_G02202_v3_SH_FEB, label="NSIDC_G02202_v3")  
plt.plot(dates_FEB[24:32], ext_ESACCI_AMSR_SH_FEB, label="ESACCI_AMSR") 
plt.plot(dates_FEB[13:30], ext_ESACCI_SSMI_SH_FEB, label="ESACCI_SSMI")  
plt.plot(dates_FEB[24:33], ext_ICDC_ASI_AMSRE_SH_FEB, label="ICDC_ASI_AMSRE")  
plt.plot(dates_FEB[13:40], ext_ICDC_ASI_SSMI_SH_FEB, label="ICDC_ASI_SSMI")  
#plt.plot(dates_FEB[10:39], ext_NOAA_OI_SH_FEB, label="NOAA_OI")  
plt.plot(dates_FEB[0:39], ext_Hadisst1_SH_FEB[109:148], label="Hadisst1") 
plt.plot(dates_FEB[0:39], ext_Hadisst2_SH_FEB[129:168], label="Hadisst2") 

# format the ticks
years = mdates.YearLocator()   # every year
#months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
#ax.xaxis.set_minor_locator(months)

plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')

# Labels
plt.title("Sea Ice Extent in Southern Hemisphere - February")
plt.ylabel("Sea Ice Extent (million sq. kms)")
plt.xlabel("Time")
plt.legend()
#plt.grid(True)

plt.show()
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh_FEB.eps', format='eps')
fig.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Sea-Ice_Conc_Extent_sh_FEB.pdf', format='pdf')
plt.close(fig)

