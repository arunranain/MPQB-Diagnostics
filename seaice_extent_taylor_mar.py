#!/usr/bin/python
# Common sea ice diags
# Arun Rana
# Taylor Diagram - https://pypi.org/project/SkillMetrics/
# https://github.com/PeterRochford/SkillMetrics/wiki

import netCDF4
import numpy as np
import sys
import scipy.stats
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
from pandas import Series, DataFrame, Panel
from seaice_commondiags import compute_extent
from seaice_commondiags import compute_area
from netCDF4 import Dataset
import pickle
import skill_metrics as sm
from sys import version_info

#NSIDC0051 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc') #lat 448 lon 304 (25*25 resolution)
conc_NSIDC0051_NH = data['siconc'][158:360:12]
area_NSIDC0051_NH = data['areacello'][:]/100
mask_NSIDC0051_NH = data['sftof'][:]
ext_NSIDC0051_NH = compute_extent(conc_NSIDC0051_NH, area_NSIDC0051_NH, threshold = 15.0, mask = mask_NSIDC0051_NH)

#OSISAF409a - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF409a_NH = data['siconc'][158:360:12]
area_OSISAF409a_NH = data['areacello'][0,:,:]
mask_OSISAF409a_NH = data['sftof'][0,:,:]
ext_OSISAF409a_NH = compute_extent(conc_OSISAF409a_NH, area_OSISAF409a_NH, threshold = 15.0, mask = mask_OSISAF409a_NH) 

#OSISAF450 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF450_NH = data['siconc'][26:228:12]
area_OSISAF450_NH = data['areacello'][0,:,:]
mask_OSISAF450_NH = data['sftof'][0,:,:]
ext_OSISAF450_NH = compute_extent(conc_OSISAF450_NH, area_OSISAF450_NH, threshold = 15.0, mask = mask_OSISAF450_NH) 

#NSIDC_G02202_v3 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc')
conc_NSIDC_G02202_v3_NH = data['siconc'][38:240:12]
area_NSIDC_G02202_v3_NH = data['areacello'][0,:,:]
mask_NSIDC_G02202_v3_NH = data['sftof'][0,:,:]
ext_NSIDC_G02202_v3_NH = compute_extent(conc_NSIDC_G02202_v3_NH, area_NSIDC_G02202_v3_NH, threshold = 15.0, mask = mask_NSIDC_G02202_v3_NH) 

#ICDC_ASI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc')
conc_ICDC_ASI_SSMI_NH = data['siconc'][2:204:12]
area_ICDC_ASI_SSMI_NH = data['areacello'][0,:,:]
mask_ICDC_ASI_SSMI_NH = data['sftof'][0,:,:]
ext_ICDC_ASI_SSMI_NH = compute_extent(conc_ICDC_ASI_SSMI_NH, area_ICDC_ASI_SSMI_NH, threshold = 15.0, mask = mask_ICDC_ASI_SSMI_NH) 

#Hadisst1 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc_Hadisst1_NH = data['siconc'][1466:1668:12]
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

#Hadisst2 - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/native/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc_Hadisst2_NH = data['siconc'][1706:1908:12]
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

#ESACCI_SSMI - NH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc') 
conc_ESACCI_SSMI_NH = data['siconc'][2::12]
area_ESACCI_SSMI_NH = data['areacello'][0,:,:]
mask_ESACCI_SSMI_NH = data['sftof'][0,:,:]
ext_ESACCI_SSMI_NH = compute_extent(conc_ESACCI_SSMI_NH, area_ESACCI_SSMI_NH, threshold = 15.0, mask = mask_ESACCI_SSMI_NH) 

# create time data and append each ext time series for PCA
ext_all = (ext_ESACCI_SSMI_NH, ext_NSIDC0051_NH, ext_OSISAF409a_NH, ext_OSISAF450_NH, ext_NSIDC_G02202_v3_NH, ext_ICDC_ASI_SSMI_NH, ext_Hadisst1_NH, ext_Hadisst2_NH)
ext_all_obs = DataFrame({'ref' : ext_ESACCI_SSMI_NH, 'pred1' : ext_NSIDC0051_NH, 'pred2' : ext_OSISAF409a_NH, 'pred3' : ext_OSISAF450_NH, 'pred4' : ext_NSIDC_G02202_v3_NH, 'pred5' : ext_ICDC_ASI_SSMI_NH, 'pred6' : ext_Hadisst1_NH, 'pred7' : ext_Hadisst2_NH})
print ext_all_obs.shape
ext_all_obs.to_pickle('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Data/ext_all_obs_NH_MAR.pkl')

#Taylor Diagram - ESACCI
def load_obj(name):
    # Load object from file in pickle format
    if version_info[0] == 2:
        suffix = 'pkl'
    else:
        suffix = 'pkl3'

    with open(name + '.' + suffix, 'rb') as f:
        return pickle.load(f) # Python2 succeeds

class Container(object): 
    
    def __init__(self, pred1, pred2, pred3, pred4, pred5, pred6, ref):
        self.pred1 = pred1
        self.pred2 = pred2
        self.pred3 = pred3
        self.pred4 = pred4
        self.pred5 = pred5
        self.pred6 = pred6
        self.pred7 = pred7
        self.ref = ref

if __name__ == '__main__':
    
    # Close any previously open graphics windows
    # ToDo: fails to work within Eclipse
    plt.close('all')
        
    # Read data from pickle file
    data = load_obj('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Data/ext_all_obs_NH_MAR')

    # Calculate statistics for Taylor diagram
    # The first array element (e.g. taylor_stats1[0]) corresponds to the 
    # reference series while the second and subsequent elements
    # (e.g. taylor_stats1[1:]) are those for the predicted series.
    taylor_stats1 = sm.taylor_statistics(data.pred1,data.ref,'data')
    taylor_stats2 = sm.taylor_statistics(data.pred2,data.ref,'data')
    taylor_stats3 = sm.taylor_statistics(data.pred3,data.ref,'data')
    taylor_stats4 = sm.taylor_statistics(data.pred4,data.ref,'data')
    taylor_stats5 = sm.taylor_statistics(data.pred5,data.ref,'data')
    taylor_stats6 = sm.taylor_statistics(data.pred6,data.ref,'data')
    taylor_stats7 = sm.taylor_statistics(data.pred7,data.ref,'data')
    
    # Store statistics in arrays
    sdev = np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1], 
                     taylor_stats2['sdev'][1], taylor_stats3['sdev'][1],
		     taylor_stats4['sdev'][1], taylor_stats5['sdev'][1],
		     taylor_stats6['sdev'][1], taylor_stats7['sdev'][1]])
    crmsd = np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1], 
                      taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1],
		      taylor_stats4['crmsd'][1], taylor_stats5['crmsd'][1],
		      taylor_stats6['crmsd'][1], taylor_stats6['crmsd'][1]])
    ccoef = np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1], 
                      taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1],
		      taylor_stats4['ccoef'][1], taylor_stats5['ccoef'][1],
		      taylor_stats6['ccoef'][1], taylor_stats6['ccoef'][1]])

    # Specify labels for points in a cell array (M1 for model prediction 1,
    # etc.). Note that a label needs to be specified for the reference even
    # though it is not used.
    label = ['Non-Dimensional Observation', 'NSIDC0051', 'OSISAF409a', 'OSISAF450', 'NSIDC_G02202_v3', 'ICDC_ASI_SSMI', 'Hadisst1', 'Hadisst2']
    sm.taylor_diagram(sdev,crmsd,ccoef, markerLabel = label,
                      markerLabelColor = 'r', markerSize = 5,
                      markerColor = 'r', markerLegend = 'on',
                      colOBS = 'y', markerobs = '*',
                      tickRMS = (0,0.2,0.4,0.6,0.8,1), tickRMSangle = 110.0,
                      colRMS = 'm', styleRMS = ':', widthRMS = 0.3, 
                      titleRMS = 'off', tickSTD = (0,0.2,0.4,0.6,0.8,1), 
                      axismax = 1.0, colSTD = 'b', styleSTD = '-.', 
                      widthSTD = 0.3, titleSTD = 'off', tickCOR = (-1,-0.9,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,0.9,1),
                      colCOR = 'k', styleCOR = '--', widthCOR = 0.3, 
                      titleCOR = 'off')

    #Show Plot
    plt.show()
    plt.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Obs_Taylor_ESACCI_Extent_nh_MAR.eps', format='eps')
    plt.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Obs_Taylor_ESACCI_Extent_nh_MAR.pdf', format='pdf')
    plt.close()
del data

######################################################################################################################################################################################################################################################################
#NSIDC0051 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc') #lat 448 lon 304 (25*25 resolution)
conc_NSIDC0051_SH = data['siconc'][157:360:12]
area_NSIDC0051_SH = data['areacello'][:]/100
mask_NSIDC0051_SH = data['sftof'][:]
ext_NSIDC0051_SH = compute_extent(conc_NSIDC0051_SH, area_NSIDC0051_SH, threshold = 15.0, mask = mask_NSIDC0051_SH)

#OSISAF409a - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF409a_SH = data['siconc'][157:360:12]
area_OSISAF409a_SH = data['areacello'][0,:,:]
mask_OSISAF409a_SH = data['sftof'][0,:,:]
ext_OSISAF409a_SH = compute_extent(conc_OSISAF409a_SH, area_OSISAF409a_SH, threshold = 15.0, mask = mask_OSISAF409a_SH) 

#OSISAF450 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc') #lat 1120 lon 760 (10*10 resolution)
conc_OSISAF450_SH = data['siconc'][25:228:12]
area_OSISAF450_SH = data['areacello'][0,:,:]
mask_OSISAF450_SH = data['sftof'][0,:,:]
ext_OSISAF450_SH = compute_extent(conc_OSISAF450_SH, area_OSISAF450_SH, threshold = 15.0, mask = mask_OSISAF450_SH) 

#NSIDC_G02202_v3 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc')
conc_NSIDC_G02202_v3_SH = data['siconc'][37:240:12]
area_NSIDC_G02202_v3_SH = data['areacello'][0,:,:]
mask_NSIDC_G02202_v3_SH = data['sftof'][0,:,:]
ext_NSIDC_G02202_v3_SH = compute_extent(conc_NSIDC_G02202_v3_SH, area_NSIDC_G02202_v3_SH, threshold = 15.0, mask = mask_NSIDC_G02202_v3_SH) 

#ICDC_ASI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc')
conc_ICDC_ASI_SSMI_SH = data['siconc'][1:204:12]
area_ICDC_ASI_SSMI_SH = data['areacello'][0,:,:]
mask_ICDC_ASI_SSMI_SH = data['sftof'][0,:,:]
ext_ICDC_ASI_SSMI_SH = compute_extent(conc_ICDC_ASI_SSMI_SH, area_ICDC_ASI_SSMI_SH, threshold = 15.0, mask = mask_ICDC_ASI_SSMI_SH) 

#Hadisst1 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc')
conc_Hadisst1_SH = data['siconc'][1465:1668:12]
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

#Hadisst2 - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ukmo/hadisst2/processed/native/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc')
conc_Hadisst2_SH = data['siconc'][1705:1908:12]
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

#ESACCI_SSMI - SH
data = Dataset('/storepelican/CLIMDATA/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc') 
conc_ESACCI_SSMI_SH = data['siconc'][1::12]
area_ESACCI_SSMI_SH = data['areacello'][0,:,:]
mask_ESACCI_SSMI_SH = data['sftof'][0,:,:]
ext_ESACCI_SSMI_SH = compute_extent(conc_ESACCI_SSMI_SH, area_ESACCI_SSMI_SH, threshold = 15.0, mask = mask_ESACCI_SSMI_SH) 

# create time data and append each ext time series for PCA
ext_all = (ext_ESACCI_SSMI_SH, ext_NSIDC0051_SH, ext_OSISAF409a_SH, ext_OSISAF450_SH, ext_NSIDC_G02202_v3_SH, ext_ICDC_ASI_SSMI_SH, ext_Hadisst1_SH, ext_Hadisst2_SH)
ext_all_obs = DataFrame({'ref' : ext_ESACCI_SSMI_SH, 'pred1' : ext_NSIDC0051_SH, 'pred2' : ext_OSISAF409a_SH, 'pred3' : ext_OSISAF450_SH, 'pred4' : ext_NSIDC_G02202_v3_SH, 'pred5' : ext_ICDC_ASI_SSMI_SH, 'pred6' : ext_Hadisst1_SH, 'pred7' : ext_Hadisst2_SH})
print ext_all_obs.shape
ext_all_obs.to_pickle('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Data/ext_all_obs_SH_FEB.pkl')

#Taylor Diagram - ESACCI
def load_obj(name):
    # Load object from file in pickle format
    if version_info[0] == 2:
        suffix = 'pkl'
    else:
        suffix = 'pkl3'

    with open(name + '.' + suffix, 'rb') as f:
        return pickle.load(f) # Python2 succeeds

class Container(object): 
    
    def __init__(self, pred1, pred2, pred3, pred4, pred5, pred6, ref):
        self.pred1 = pred1
        self.pred2 = pred2
        self.pred3 = pred3
        self.pred4 = pred4
        self.pred5 = pred5
        self.pred6 = pred6
        self.pred7 = pred7
        self.ref = ref

if __name__ == '__main__':
    
    # Close any previously open graphics windows
    # ToDo: fails to work within Eclipse
    plt.close('all')
        
    # Read data from pickle file
    data = load_obj('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Data/ext_all_obs_SH_FEB')

    # Calculate statistics for Taylor diagram
    # The first array element (e.g. taylor_stats1[0]) corresponds to the 
    # reference series while the second and subsequent elements
    # (e.g. taylor_stats1[1:]) are those for the predicted series.
    taylor_stats1 = sm.taylor_statistics(data.pred1,data.ref,'data')
    taylor_stats2 = sm.taylor_statistics(data.pred2,data.ref,'data')
    taylor_stats3 = sm.taylor_statistics(data.pred3,data.ref,'data')
    taylor_stats4 = sm.taylor_statistics(data.pred4,data.ref,'data')
    taylor_stats5 = sm.taylor_statistics(data.pred5,data.ref,'data')
    taylor_stats6 = sm.taylor_statistics(data.pred6,data.ref,'data')
    taylor_stats7 = sm.taylor_statistics(data.pred7,data.ref,'data')
    
    # Store statistics in arrays
    sdev = np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1], 
                     taylor_stats2['sdev'][1], taylor_stats3['sdev'][1],
		     taylor_stats4['sdev'][1], taylor_stats5['sdev'][1],
		     taylor_stats6['sdev'][1], taylor_stats7['sdev'][1]])
    crmsd = np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1], 
                      taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1],
		      taylor_stats4['crmsd'][1], taylor_stats5['crmsd'][1],
		      taylor_stats6['crmsd'][1], taylor_stats6['crmsd'][1]])
    ccoef = np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1], 
                      taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1],
		      taylor_stats4['ccoef'][1], taylor_stats5['ccoef'][1],
		      taylor_stats6['ccoef'][1], taylor_stats6['ccoef'][1]])

    # Specify labels for points in a cell array (M1 for model prediction 1,
    # etc.). Note that a label needs to be specified for the reference even
    # though it is not used.
    label = ['Non-Dimensional Observation', 'NSIDC0051', 'OSISAF409a', 'OSISAF450', 'NSIDC_G02202_v3', 'ICDC_ASI_SSMI', 'Hadisst1', 'Hadisst2']
    sm.taylor_diagram(sdev,crmsd,ccoef, markerLabel = label,
                      markerLabelColor = 'r', markerSize = 5,
                      markerColor = 'r', markerLegend = 'on',
                      colOBS = 'y', markerobs = '*',
                      tickRMS = (0,0.2,0.4,0.6,0.8,1), tickRMSangle = 110.0,
                      colRMS = 'm', styleRMS = ':', widthRMS = 0.3, 
                      titleRMS = 'off', tickSTD = (0,0.2,0.4,0.6,0.8,1), 
                      axismax = 1.0, colSTD = 'b', styleSTD = '-.', 
                      widthSTD = 0.3, titleSTD = 'off', tickCOR = (-1,-0.9,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,0.9,1),
                      colCOR = 'k', styleCOR = '--', widthCOR = 0.3, 
                      titleCOR = 'off')

    #Show Plot
    plt.show()
    plt.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Obs_Taylor_ESACCI_Extent_sh_FEB.eps', format='eps')
    plt.savefig('/home/elic/arunr/git/UCL_BE_AR/Diagnostics/Results/Obs_Taylor_ESACCI_Extent_sh_FEB.pdf', format='pdf')
    plt.close()
del data

