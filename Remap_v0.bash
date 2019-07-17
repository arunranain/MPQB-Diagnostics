#!/bin/bash
# 
# Arun Rana and Francois Massonnet
# 01/07/2019
# Remapping the nc files to 1*1 degree grid for intercomparison in MPQB

set -o nounset
set -o errexit
set -x

# Hadisst1 -  We remap all the variables that are needed for the comparison (In the following example only 1 variables to be remapped)
filein=${YOUR_DIR}/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc
rootout=${YOUR_DIR}/obs/ice/siconc/ukmo/hadisst1/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc

# Hadisst1 -  We remap all the variables that are needed for the comparison (In the following example we have 3 different variables to be remapped)
filein=${PATH_TO_YOUR_DIR}/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc
rootout=${PATH_TO_YOUR_DIR}/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncatted -a coordinates,areacello,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,areacello $filein regridded2.nc
ncatted -a coordinates,sftof,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,sftof $filein regridded3.nc
ncks -F -A -v siconc regridded.nc regridded4.nc
ncks -F -A -v areacello regridded2.nc regridded4.nc
ncks -F -A -v sftof     regridded3.nc regridded4.nc
ncrcat -F regridded4.nc ${rootout}/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc

