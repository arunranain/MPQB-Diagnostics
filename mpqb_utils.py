import xarray as xr
import matplotlib.pyplot as plt
import xesmf as xe
import intake
import matplotlib
import cartopy.crs as ccrs
import os
import glob

def homogenize_mpqb_dataset(indat,name=None):
    if name=='era_interim_land':
        indat = indat.rename({'longitude' : 'lon', 'latitude': 'lat', 'swvl1' : 'sm'})
        indat = indat.assign_coords(lon=(((indat.lon + 180) % 360) - 180))
    elif name=='reanalysis_era5_single_levels':
        indat = indat.rename({'longitude' : 'lon', 'latitude': 'lat', 'swvl1' : 'sm'})
        indat = indat.assign_coords(lon=(((indat.lon + 180) % 360) - 180))
    elif name=='reanalysis_uerra_europe_soil_levels_uerra_harmonie':
        indat = indat.rename({'latitude' : 'lat', 'longitude': 'lon', 'vsw' : 'sm'})
        indat = indat.assign_coords(lon=(((indat.lon + 180) % 360) - 180))
    elif name=='satellite_soil_moisture':
        pass
    else:
        raise ValueError
    return indat

def regrid_mpqb(input_ds,varname,dest_grid,regionalgrid=None,regridmethod='bilinear',**kwargs):
    '''
    Function for regridding datasets for the MPQB. If regionalgrid
    is True, a small offset is added, such that it is possible to mask
    any datapoints outside the regional grid. See github issue # 
    
    Parameters
    ----------
      input_ds : xarray.Dataset
         the dataset to be regridded
      varname : str
         the name of the variable to be selected for regridding
      destgrid : xarray.Dataset
         the destination grid 
      regionalgrid : bool
         boolean to indicate if this is a regional grid or not
      regridmethod : passed on to xe.Regridder
    
    '''
    input_da = input_ds[varname]
    if regionalgrid:
        print("Start regridding for regional grid")
        constantval = 10
        input_da = input_da + constantval
        assert(int((input_da==0.).sum())==0) # Make sure that there are no zero's in the data, since they
        # will be masked out...
        regridder = xe.Regridder(input_ds,dest_grid,regridmethod,**kwargs)
        da_out = regridder(input_da)
        da_out = da_out.where(da_out!=0.)
        da_out = da_out - constantval
    else:
        print("Start regridding for global grid")
        regridder = xe.Regridder(input_ds,dest_grid,regridmethod,**kwargs)
        da_out = regridder(input_da)
    return da_out
        