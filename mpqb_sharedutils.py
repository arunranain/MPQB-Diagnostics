from sharedutils import *
from diag1d import *
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def mpqb_pearsonr(da1,da2):
    rval,pval = parallel_apply_along_axis(pearsonr1d,0,(da1.values,da2.values))
    template = da1[0:1].mean('time')
    rval_da = template.copy()
    pval_da = template.copy()
    rval_da.values = rval
    pval_da.values = pval
    rval_da.name = 'corr'
    pval_da.name = 'p-value'
    rval_da.attrs['datasetname'] = da1.attrs['datasetname']+'-'+da2.attrs['datasetname']
    pval_da.attrs['datasetname'] = da1.attrs['datasetname']+'-'+da2.attrs['datasetname']
    return rval_da,pval_da

def mpqb_rmsd(da1,da2):
    rmsdval = parallel_apply_along_axis(rmsd1d,0,(da1.values,da2.values))
    template = da1[0:1].mean('time')
    rmsd_da = template.copy()
    rmsd_da.values = rmsdval
    rmsd_da.name = 'RMSD of '+rmsd_da.name
    rmsd_da.attrs['datasetname'] = da1.attrs['datasetname']+'-'+da2.attrs['datasetname']
    return rmsd_da

def mpqb_absdiff(da1,da2):
    absdiffval = parallel_apply_along_axis(absdiffaxismean1d,0,(da1.values,da2.values))
    template = da1[0:1].mean('time')
    absdiff_da = template.copy()
    absdiff_da.values = absdiffval
    absdiff_da.name = 'Abs difference in '+absdiff_da.name
    absdiff_da.attrs['datasetname'] = da1.attrs['datasetname']+'-'+da2.attrs['datasetname']
    return absdiff_da

def mpqb_reldiff(da1,da2):
    reldiffval = parallel_apply_along_axis(reldiffaxismean1d,0,(da1.values,da2.values))
    template = da1[0:1].mean('time')
    reldiff_da = template.copy()
    reldiff_da.values = reldiffval
    reldiff_da.name = 'Rel. difference in '+reldiff_da.name
    reldiff_da.attrs['datasetname'] = da1.attrs['datasetname']+'-'+da2.attrs['datasetname']
    return reldiff_da

def mpqb_mankendall(da):
    mk = parallel_apply_along_axis(mannkendall1d,0,da.values)
    template = da[0:1].mean('time')
    mk_da = template.copy()
    mk_da.values = mk
    mk_da.name = 'Mann-Kendall test on '+mk_da.name
    mk_da.attrs['datasetname'] = da.attrs['datasetname']
    return mk_da

def mpqb_plot(dataset,cmaptype='values',robust=True,**kwargs):
    # Get right cmap
    # This is quite an ugly implementation, but it works.
    if 'cmapname_diverging' in kwargs and cmaptype=='diverging':
        mpqbcmap = matplotlib.cm.get_cmap(kwargs.pop('cmapname_diverging'))
        _ = kwargs.pop('cmapname_values')
    elif 'cmapname_values' in kwargs and cmaptype=='values':
        mpqbcmap = matplotlib.cm.get_cmap(kwargs.pop('cmapname_values'))
        # Also pop the other
        _ = kwargs.pop('cmapname_diverging')
    else:
        print("Warning, taking default cmap, please specify")
        mpqbcmap = matplotlib.cm.get_cmap()
    mpqbcmap.set_bad(color='grey')
    
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_subplot(111, projection=ccrs.Robinson())
    dataset.plot(ax=ax, transform=ccrs.PlateCarree(),robust=robust,cmap=mpqbcmap,**kwargs)
    ax.coastlines()
    if 'datasetname' in dataset.attrs:
        plt.title(dataset.attrs['datasetname'])
    
    

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
