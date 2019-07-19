import itertools as it
import os
import sys

import cartopy.crs as ccrs
import dask
import matplotlib
matplotlib.use('agg') # related to https://github.com/ipython/ipython/issues/10627
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/home/crezees/projects/MPQB/c3s_511_mpqb/')
sys.path.append('/home/crezees/projects/trend_lim_framework/')
from diag1d import *
from intake import open_catalog
from mpqb_sharedutils import *
from sharedutils import *

# MPQB conventions, dimension coordinate names:
# - lat , lon, time
# - lon [-180,180]

# Some parameters
figtype = 'png' # png/pdf/... , used as filename extension in plt.savefig()
plotdir = '/home/crezees/projects/MPQB/c3s_511_mpqb/plotdir/'
refdataset = 'satellite_soil_moisture'
varname = 'sm'
plotkwargs = {
    # Specify colormaps
    # (Note that these are colorblind friendly: https://github.com/matplotlib/matplotlib/issues/7081/)
    'cmapname_diverging' : 'RdYlBu_r',
    'cmapname_values' : 'YlOrBr'
}

# When using Python Intake, specify catalogfile
catalogfile = '/home/crezees/projects/earthdata_reader/catalogs/mpqb_sm_processed.yml'

# Specify both the 'long names' (corresponding to the catalog file) 
# and 'short names' of the datasets.
dataset_shortnames = {
    'era_interim_land' : 'EI Land',
    'reanalysis_era5_single_levels' : 'ERA5',
    'satellite_soil_moisture' : 'ESA CCI',
    'reanalysis_uerra_europe_soil_levels_uerra_harmonie' : 'UERRA'
}

# Open the catalog
mpqb_cat = open_catalog(catalogfile)

# Read the data
mpqb_keys = list(mpqb_cat)
print("Loading :",'; '.join(mpqb_keys))

# For evaluators: this is the entry point into the framework. Here we read in the data
# using Python intake, but another approach is possible. The result needs to be 
# a dictionary with as keys the (long) names of the datasets and as 
# values the data as Xarray DataArrays.
mpqb_datasets = {key : mpqb_cat[key].to_dask()[varname] for key in mpqb_keys}

###########################################################################
###########################################################################
###########################################################################

# Assert that shapes of MPQB datasets are identical
try:
    assert(len(set([dset.shape for _,dset in mpqb_datasets.items()]))==1)
    print("Input shapes look fine")
except AssertionError:
    print("Warning: input shapes are not identical")
    for key,val in mpqb_datasets.items():
        print(key," has shape: ",val.shape)

# This creates a list of length 2 tuples specifying all combinations 
# for comparison (including ref to ref; which is handled as a seperate case)
#TODO make mpqb_datasets an ordered dictionary
comparisons = list(it.product([key for key in mpqb_datasets],[refdataset]))

all_datasets = [key for key in mpqb_datasets]

# Set the short names
for key,val in mpqb_datasets.items():
    mpqb_datasets[key].attrs['datasetname'] = dataset_shortnames[key]

print("Starting MPQB processing")
print("Starting pearsonr")
#TODO add for each of below the case that dsa=dsb 
# (this is the 'special' plot out of the multi panel showing the REF dataset)
for dsa,dsb in comparisons:
    if dsa!=dsb:
        corr,_ = mpqb_pearsonr(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(corr,cmaptype='diverging',**plotkwargs)
        savename = os.path.join(plotdir,'pearsonr_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting rmsd")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        rmsd = mpqb_rmsd(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(rmsd)
        savename = os.path.join(plotdir,'rmsd_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting absdiff")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        absdiff = mpqb_absdiff(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(absdiff,cmaptype='diverging')
        savename = os.path.join(plotdir,'absdiff_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting reldiff")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        reldiff = mpqb_reldiff(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(reldiff,cmaptype='diverging')
        savename = os.path.join(plotdir,'reldiff_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass