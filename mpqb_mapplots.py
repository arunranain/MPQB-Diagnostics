import itertools as it
import os
import sys

import cartopy.crs as ccrs
import dask
import matplotlib
matplotlib.use('agg') # related to https://github.com/ipython/ipython/issues/10627
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.append('/home/crezees/projects/MPQB/c3s_511_mpqb/')
sys.path.append('/home/crezees/projects/trend_lim_framework/')
from diag1d import *
from intake import open_catalog
from mpqb_sharedutils import *
from sharedutils import *
import argparse
import logging
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("cfgfile")
parser.add_argument("filestream")
args = parser.parse_args()

logging.basicConfig()
logging.root.setLevel(logging.INFO)
logger = logging.getLogger('mpqb-preprocessing')
logger.setLevel(logging.DEBUG)
logger.info("Using following cfg-file: {0}".format(args.cfgfile))

with open(args.cfgfile,'r') as handle:
    cfg = yaml.safe_load(handle)

# MPQB conventions, dimension coordinate names:
# - lat , lon, time
# - lon [-180,180]
# - only one variable on the file

# Some parameters
figtype = 'png' # png/pdf/... , used as filename extension in plt.savefig()

sessiondir = os.path.join(cfg['basedir'],cfg['sessionname'])
plotdir = os.path.join(sessiondir,'plots')
logger.info('Creating directory: {0}'.format(plotdir))
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
else:
    logger.warning("Dir exists already, be aware of overwriting: {0}".format(plotdir))
plotdir_stream = os.path.join(plotdir,args.filestream)
if not os.path.exists(plotdir_stream):
    os.mkdir(plotdir_stream)
else:
    logger.warning("Dir exists already, be aware of overwriting: {0}".format(plotdir_stream))

    
plotkwargs = {
    # Specify colormaps
    # (Note that these are colorblind friendly: https://github.com/matplotlib/matplotlib/issues/7081/)
    'cmapname_diverging' : 'RdYlBu_r',
    'cmapname_values' : 'YlOrBr'
}

dataset_shortnames = cfg['dataset_shortnames']
refdataset = cfg['reference']

# For evaluators: this is the entry point into the framework. Here we read in the data
# using Python intake, but another approach is possible. The result needs to be 
# a dictionary with as keys the (long) names of the datasets and as 
# values the data as Xarray DataArrays.
# Read in the data 
mpqb_datasets = {}
for dataset in cfg['datasets']:
    inputdir = os.path.join(sessiondir,args.filestream)
    inputdir_dataset = os.path.join(inputdir,dataset)
    logger.info("Reading in data from: {0}".format(inputdir_dataset))

    infilepattern = os.path.join(inputdir_dataset,'*.nc')
    mpqb_datasets[dataset] = xr.open_mfdataset(infilepattern)[cfg['varname']]
    
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
        savename = os.path.join(plotdir_stream,'pearsonr_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting rmsd")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        rmsd = mpqb_rmsd(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(rmsd)
        savename = os.path.join(plotdir_stream,'rmsd_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting absdiff")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        absdiff = mpqb_absdiff(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(absdiff,cmaptype='diverging')
        savename = os.path.join(plotdir_stream,'absdiff_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting reldiff")
for dsa,dsb in comparisons:
    if dsa!=dsb:
        reldiff = mpqb_reldiff(mpqb_datasets[dsa],mpqb_datasets[dsb])
        mpqb_plot(reldiff,cmaptype='diverging')
        savename = os.path.join(plotdir_stream,'reldiff_'+mpqb_datasets[dsa].attrs['datasetname'].replace(" ","")+'-'+mpqb_datasets[dsb].attrs['datasetname'].replace(" ","")+'.'+figtype)
        print(savename)
        plt.savefig(savename)
    else:
        pass

print("Starting mktrends")
for datasetname,dataset in mpqb_datasets.items():
    dataset_yearly = dataset.resample({'time' : 'Y'}).mean(dim='time').compute()
    # Put back the attributes.
    dataset_yearly.attrs = dataset.attrs
    mkresult = mpqb_mankendall(dataset_yearly)
    mpqb_plot(mkresult,cmaptype='diverging')
    savename = os.path.join(plotdir_stream,'mktrend_'+datasetname+'.'+figtype)
    print(savename)
    plt.savefig(savename)

    
#Copyright (C) 2019, Bas Crezee, ETH Zürich, Institut for Atmospheric and Climate Science
#Distributed under GNU General Public License GPL-3.0
#
#This program is free software: you can redistribute it and/or modify it under the terms 
#of the GNU General Public License as published by the Free Software Foundation, either 
#version 3 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#See the GNU General Public License for more details.
# 
#Please acknowledge ETH Zurich, Institute for Atmospheric and Climate Science, Bas Crezee in any further use of the software.#
#
#See https://www.gnu.org/licenses/
