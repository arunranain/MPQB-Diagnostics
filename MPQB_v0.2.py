#!/usr/bin/python
#
# Arun Rana and Francois Massonnet
# 01/07/2019

# Edit History
# 13/09/2019 - Plynomial fit changed to Linregress (instead of polyfit)
# 11/10/2019 - Projection changed to Robinson based on comment from Peter
# 25/11/2019 - In trend section changed conc to Var1, Var2, Var3, and VarRef after comments from Yassmina
# 25/11/2019 - in readme section below point 4, added changes that should be made to m.fillcontinents(color='gray') if it is land dataset after comments from Yassmina
"""
README:
The following script is meant MPQB Diagnostics and plotting. It has been divided into 2 parts i.e spatial comparison diagnostics (averaged over time length) and temporal comparison diagnostics (averagered over space).
The below space is used to define what/where you as ECV evaluator has to change in order to run it for your ECV. The following should be modified according to needs of your ECV:

1. Install skill_metrics module - $ pip install SkillMetrics
2. Input Data Information - provide links to your datasets that are interpolated with "Remap_v0.bash" or similar. You need to define Reference, Dataset1, Dataset2 and so on.
P.S. - The script input data part assumes your variable is 3 dimensional (time, lat, lon) and not otherewise. If it is more than 3 dimension please collapse/squeeze against the dimension that leaves you with time/lat/lon.
3. In all spatial plots please use your ECV specific/desired projection for mapping. In the following script we have used robinator. Other Projection are listed and linked in documentation.
4. I will leave other lat/lon information same unless compelling to keep the plots homogenised and same goes for filling continents and coastlines. Where you have land dataset you need to change m.fillcontinents(color='gray') to not fill land as grey but with your dataset.
5. In all spatial plots you need to adjust colormap according to ECV in cmap='yourcolor' in e.g. m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
6. Similarly you need to define the range for the displayed variable accoding to diagnostic presented by chnaging plt.clim range values in e.g. plt.clim(0, 100)
7. Please defibne the subplot titles and plot tiles as per need and convinience in all the plots.
8. Means and variability - Annual Mean - Position of the colorbar (cax1 = fig.add_axes([0.9, 0.2, 0.02, 0.6])) in the plot needs to be adjusted to linking and based on number of datasets in the plot.
9. Means and variability - grid-point correlation (and in all futher spatial plots) - Both position of the colorbars need to be fixed based on the plot/datasets in cax1 and cax2.
10. Time series of global/hemisphere annual/monthly mean (line plot) - Start and end year needs to be defined and changed. If daily dataset then we need to account for days as well. All the labels/titles/linecolor(ranges and position) needs to be defined.

"""
## Import all required modules*******************************************************************************************************************************************************************
import numpy as np
from scipy.stats import pearsonr
from scipy.stats.mstats import linregress
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import datetime as dt
from pandas import DataFrame
from netCDF4 import Dataset
import os
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
import pickle
import skill_metrics as sm
from sys import version_info
#from seaice_commondiags import compute_extent - This is specific function for calculation of SeaIce extent from concentration and is not needed for everyone 

## Input dataset information - All datasets are interpolated to 1*1 grid with remap script for comparison****************************************************************************************
Num_Obs = 4
Reference = Dataset('/PATH_TO_YOUR_DIR/interp/Reference.nc')
lat = Reference.variables['latitude'][:]
lon = Reference.variables['longitude'][:]
lon, lat = np.meshgrid(lon,lat)
VarRef = Reference.variables['varname']
del Reference
Dataset1 = Dataset('/PATH_TO_YOUR_DIR/interp/Dataset1.nc')
Var1 = Dataset1.variables['varname']
del Dataset1
Dataset2 = Dataset('/PATH_TO_YOUR_DIR/interp/Dataset2.nc')
Var2 = Dataset2.variables['varname']
del Dataset2
Dataset3 = Dataset('/PATH_TO_YOUR_DIR/interp/Dataset3.nc')
Var3 = Dataset3.variables['varname']
del Dataset3


## Spatial comparison MPQB diagnostics***********************************************************************************************************************************************************

# Means and variability - Annual Mean -----------------------------------------------------------------------------------------------------------------------------------------------------------
PlotVariable = np.nanmean(Var1, axis=0)
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c') #cyl - Equidistant cylindrical projection
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(Var2, axis=0)
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(Var3, axis=0)
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(VarRef, axis=0)
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image1, cax1, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureMeans.jpg', format='jpg')
plt.close(fig)

# Means and variability - grid-point correlation --------------------------------------------------------------------------------------------------------------------------------------------------
PlotVariable = np.zeros((Var1.shape[1], Var1.shape[2]))
for n in range(Var1.shape[1]):
   for m in range(Var1.shape[2]):
      PlotVariable[n, m] = pearsonr(Var1[:,n,m], VarRef[:,n,m])[0]
PlotVariable[PlotVariable==0.] = np.nan
PlotVariable = np.ma.array(PlotVariable,mask=np.isnan(PlotVariable))
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image2=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var2.shape[1], Var2.shape[2]))
for n in range(Var2.shape[1]):
   for m in range(Var2.shape[2]):
      PlotVariable[n, m] = pearsonr(Var2[:,n,m], VarRef[:,n,m])[0]
PlotVariable[PlotVariable==0.] = np.nan
PlotVariable = np.ma.array(PlotVariable,mask=np.isnan(PlotVariable))
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var3.shape[1], Var3.shape[2]))
for n in range(Var3.shape[1]):
   for m in range(Var3.shape[2]):
      PlotVariable[n, m] = pearsonr(Var3[:,n,m], VarRef[:,n,m])[0]
PlotVariable[PlotVariable==0.] = np.nan
PlotVariable = np.ma.array(PlotVariable,mask=np.isnan(PlotVariable))
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(VarRef, axis=0)
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image1, cax1, orientation='horizontal', label='% (Fraction)')
cax2 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image2, cax2, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureCorrelations.jpg', format='jpg')
plt.close(fig)

# Means and variability - grid-point RMSD -------------------------------------------------------------------------------------------------------------------------------------------------------
PlotVariable = np.zeros((Var1.shape[1], Var1.shape[2]))
for n in range(Var1.shape[1]):
   for m in range(Var1.shape[2]):
      PlotVariable[n, m] = np.sqrt(((Var1[:,n,m] - VarRef[:,n,m]) ** 2).mean())
PlotVariable[PlotVariable==0.] = np.nan
PlotVariable = np.ma.array(PlotVariable,mask=np.isnan(PlotVariable))
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image2=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Reds')
plt.clim(0, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var2.shape[1], Var2.shape[2]))
for n in range(Var2.shape[1]):
   for m in range(Var2.shape[2]):
      PlotVariable[n, m] = np.sqrt(((Var2[:,n,m] - VarRef[:,n,m]) ** 2).mean())
PlotVariable[PlotVariable==0.] = np.nan
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Reds')
plt.clim(0, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var3.shape[1], Var3.shape[2]))
for n in range(Var3.shape[1]):
   for m in range(Var3.shape[2]):
      PlotVariable[n, m] = np.sqrt(((Var3[:,n,m] - VarRef[:,n,m]) ** 2).mean())
PlotVariable[PlotVariable==0.] = np.nan
PlotVariable = np.ma.array(PlotVariable,mask=np.isnan(PlotVariable))
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Reds')
plt.clim(0, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(VarRef, axis=0)
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image1, cax1, orientation='horizontal', label='% (Fraction)')
cax2 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image2, cax2, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureRMSD.jpg', format='jpg')
plt.close(fig)

# Trends - grid-point linear trends - Linregress --------------------------------------------------------------------------------------------------------------------------------
PlotVariable = np.zeros((Var1.shape[1], Var1.shape[2]))
for n in range(Var1.shape[1]):
   for m in range(Var1.shape[2]):
      if np.sum(np.logical_not(np.isnan(Var1[:,n,m])))>=2:
         PlotVariable[n, m] = linregress(range(len(Var1[:,n,m])),Var1[:,n,m]).slope
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image2=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var2.shape[1], Var2.shape[2]))
for n in range(Var2.shape[1]):
   for m in range(Var2.shape[2]):
      if np.sum(np.logical_not(np.isnan(Var2[:,n,m])))>=2:
         PlotVariable[n, m] = linregress(range(len(Var2[:,n,m])),Var2[:,n,m]).slope
#fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((Var3.shape[1], Var3.shape[2]))
for n in range(Var3.shape[1]):
   for m in range(Var3.shape[2]):
      if np.sum(np.logical_not(np.isnan(Var3[:,n,m])))>=2:
         PlotVariable[n, m] = linregress(range(len(Var3[:,n,m])),Var3[:,n,m]).slope
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.zeros((VarRef.shape[1], VarRef.shape[2]))
for n in range(VarRef.shape[1]):
   for m in range(VarRef.shape[2]):
      if np.sum(np.logical_not(np.isnan(VarRef[:,n,m])))>=2:
         PlotVariable[n, m] = linregress(range(len(VarRef[:,n,m])),VarRef[:,n,m]).slope
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu')
plt.clim(-1, 1)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image1, cax1, orientation='horizontal', label='% (Fraction)')
cax2 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image2, cax2, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureTrend.jpg', format='jpg')
plt.close(fig)

# Differences - grid-point Absolute differences -----------------------------------------------------------------------------------------------------------------------------------------------
PlotVariable = np.nanmean(Var1, axis=0) - np.nanmean(VarRef, axis=0)
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image2=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(Var2, axis=0) - np.nanmean(VarRef, axis=0)
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(Var3, axis=0) - np.nanmean(VarRef, axis=0)
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(VarRef, axis=0)
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image1, cax1, orientation='horizontal', label='% (Fraction)')
cax2 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image2, cax2, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureRMSD.jpg', format='jpg')
plt.close(fig)

# Differences - grid-point Relative differences -----------------------------------------------------------------------------------------------------------------------------------------------
PlotVariable = (np.nanmean(Var1, axis=0) - np.nanmean(VarRef, axis=0))/(np.nanmean(VarRef, axis=0)) *100
fig = plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image2=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = (np.nanmean(Var2, axis=0) - np.nanmean(VarRef, axis=0))/(np.nanmean(VarRef, axis=0)) *100
plt.subplot(2,2,2)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = (np.nanmean(Var3, axis=0) - np.nanmean(VarRef, axis=0))/(np.nanmean(VarRef, axis=0)) *100
plt.subplot(2,2,3)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='RdBu_r')
plt.clim(-10, 10)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

PlotVariable = np.nanmean(VarRef, axis=0)
plt.subplot(2,2,4)
m = Basemap(projection='robin', lon_0=0,resolution='c')
m.drawcoastlines(color='lightgray', linewidth=0.2)
m.fillcontinents(color='gray')
m.drawparallels(np.arange(-90.,120.,30.), labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.), labels=[0,0,0,1] )
image1=m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
plt.clim(0, 100)
plt.title("YOUR_SUBPLOT_TITLE")
del PlotVariable

cax1 = fig.add_axes([0.7, 0.075, 0.175, 0.02])
plt.colorbar(image1, cax1, orientation='horizontal', label='% (Fraction)')
cax2 = fig.add_axes([0.9, 0.2, 0.02, 0.6])
plt.colorbar(image2, cax2, orientation='vertical', label='% (Fraction)')
plt.suptitle("YOUR_PLOT_TITLE")
fig.savefig('PATH_TO_YOUR_DIR/FigureRMSD.jpg', format='jpg')
plt.close(fig)

## Temporal comparison MPQB diagnostics********************************************************************************************************************************************************
# Average time series over lat/lon (space)
Avg_VarRef = np.nanmean(VarRef, axis=1)
Avg_VarRef = np.nanmean(Avg_VarRef, axis=1)
Avg_Var1 = np.nanmean(Var1, axis=1)
Avg_Var1 = np.nanmean(Avg_Var1, axis=1)
Avg_Var2 = np.nanmean(Var2, axis=1)
Avg_Var2 = np.nanmean(Avg_Var2, axis=1)
Avg_Var3 = np.nanmean(Var3, axis=1)
Avg_Var3 = np.nanmean(Avg_Var3, axis=1)

avg_all = (Avg_VarRef, Avg_Var1, Avg_Var2, Avg_Var3)
avg_all_obs = DataFrame({'ref' : Avg_VarRef, 'pred1' : Avg_Var1, 'pred2' : Avg_Var2, 'pred3' : Avg_Var3})
print (avg_all_obs.shape)
avg_all_obs.to_pickle('PATH_TO_YOUR_DIR/avg_all_obs.pkl')

# Time series of global/hemisphere annual/monthly mean (line plot) ----------------------------------------------------------------------------------------------------------------------------
dates = []
for years in range(1977, 2017):
    for months in range(1, 13):
        dates.append(dt.datetime(year=years, month=months, day=1))
fig = plt.figure()
ax = plt.axes()
plt.plot(dates[0:468], Avg_VarRef, label="REFERENCE")
plt.plot(dates[0:468], Avg_Var1, label="OBS1")
plt.plot(dates[0:468], Avg_Var2, label="OBS2")
plt.plot(dates[0:468], Avg_Var3, label="OBS3")

years = mdates.YearLocator()
yearsFmt = mdates.DateFormatter('%Y')
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
plt.setp(ax.get_xticklabels(), rotation=70, horizontalalignment='right')
plt.title("YOUR_PLOT_TITLE")
plt.ylabel("YOUR_PLOT_YLABEL")
plt.xlabel("Years")
plt.tight_layout(rect=[0, 0, 0.75, 1])
plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", borderaxespad=0, mode='expand', ncol=1)
plt.show()
fig.savefig('PATH_TO_YOUR_DIR/FigureTimeseries.jpg', format='jpg')
plt.close(fig)

# Taylor Diagram ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def load_obj(name):
    # Load object from file in pickle format
    if version_info[0] == 3:
        suffix = 'pkl'
    else:
        suffix = 'pkl2'
    with open(name + '.' + suffix, 'rb') as f:
        return pickle.load(f) # Python3 succeeds

class Container(object): 
    
    def __init__(self, pred1, pred2, pred3, ref):
        self.pred1 = pred1
        self.pred2 = pred2
        self.pred3 = pred3
        self.ref = ref

if __name__ == '__main__':
    
    # Close any previously open graphics windows
    plt.close('all')
        
    # Read data from pickle file
    data = load_obj('PATH_TO_YOUR_DIR/avg_all_obs.pkl')

    # Calculate statistics for Taylor diagram (The first array element (e.g. taylor_stats1[0]) corresponds to the reference series while the second and subsequent elements (e.g. taylor_stats1[1:]) are those for the predicted series.
    taylor_stats1 = sm.taylor_statistics(data.pred1,data.ref,'data')
    taylor_stats2 = sm.taylor_statistics(data.pred2,data.ref,'data')
    taylor_stats3 = sm.taylor_statistics(data.pred3,data.ref,'data')
    
    # Store statistics in arrays
    sdev = np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1], 
                     taylor_stats2['sdev'][1], taylor_stats3['sdev'][1]])
    crmsd = np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1], 
                      taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1]])
    ccoef = np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1], 
                      taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1]])

    # Specify labels for points in a cell array (M1 for model prediction 1, etc.). Note that a label needs to be specified for the reference even though it is not used.
    label = ['Non-Dimensional Observation', 'OBS1', 'OBS2', 'OBS3']
    sm.taylor_diagram(sdev,crmsd,ccoef, markerLabel = label,
                      markerLabelColor = 'r', markerSize = 5,
                      markerColor = 'r', markerLegend = 'on',
                      colOBS = 'y', markerobs = '*',
                      tickRMS = (0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), tickRMSangle = 110.0,
                      colRMS = 'm', styleRMS = ':', widthRMS = 0.3, 
                      titleRMS = 'on', tickSTD = (0,0.2,0.4,0.6,0.8,1,1.2,1.4), 
                      axismax = 1.4, colSTD = 'b', styleSTD = '-.', 
                      widthSTD = 0.3, titleSTD = 'on', tickCOR = (-1,-0.9,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,0.9,1),
                      colCOR = 'k', styleCOR = '--', widthCOR = 0.3, 
                      titleCOR = 'on')
    plt.show()
    plt.savefig('PATH_TO_YOUR_DIR/FigureTaylor.jpg', format='jpg')
    plt.close()
del data
