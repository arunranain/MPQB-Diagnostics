#!/bin/bash
# 
# Arun Rana
# Remapping the nc files to 1*1 degree grid same as Hadisst1

set -o nounset
set -o errexit
set -x

# NSIDC0051 - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/
mkdir -p $rootout
cdo sinfov $filein                         #getting information about the variables and checking if the generic co-ordinates are defined on desired variable
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein  #CDO can't find the generic cordinates because these variables are not assigned (checking above) and thus use the NCO command to add the missing attribute.  
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_r1i1p1_mon_197901-201712_nh-psn25.nc

# OSI-409-a - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_OSI-409a_r1i1p1_197901-201512_nh.nc

# OSI-450 - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_OSI-450_r1i1p1_199001-201512_nh.nc

# NSIDC_G02202_v3 - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_nh.nc

# ESACCI_AMSR - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/AMSR/processed/native/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_NH.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/AMSR/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_NH.nc

# ESACCI_SSMI - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/SSMI/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_NH.nc

# ICDC_ASI_AMSRE - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-AMSRE/processed/native/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_n.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-AMSRE/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_n.nc

# ICDC_ASI_SSMI - NH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_arc.nc

# NOAA_OI
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NOAA/NOAA-OI/processed/native/siconc_SImon_NOAA_OI_r1i1p1_198912-201809_reg1.0.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NOAA/NOAA-OI/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_NOAA_OI_r1i1p1_198912-201809_reg1.0.nc

# Hadisst1 -  We are doing mask and area as well so we can use to calculate extent
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ukmo/hadisst1/processed/native/siconc_SImon_ukmo_HadISST1_r1i1p1_187001-201805_reg1.0.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ukmo/hadisst1/processed/interp/
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

# Hadisst2
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ukmo/hadisst2/processed/native/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ukmo/hadisst2/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ukmo_HadISST2_r1i1p1_185001-201803_reg1.0.nc

# NSIDC0051 - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/NSIDC-0051/processed/native/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/NSIDC-0051/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_r1i1p1_mon_197901-201712_sh-pss25.nc

# OSI-409-a - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/native/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-409-a/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_OSI-409a_r1i1p1_197901-201512_sh.nc

# OSI-450 - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/OSI-SAF/OSI-450/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_OSI-450_r1i1p1_199001-201512_sh.nc

# NSIDC_G02202_v3 - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/G02202_v3/processed/native/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/NSIDC/G02202_v3/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_NSIDC_G02202_v3_r1i1p1_198901-201712_sh.nc

# ESACCI_AMSR - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/AMSR/processed/native/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_SH.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/AMSR/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ESACCI_AMSR_r1i1p1_200301-201012_SH.nc

# ESACCI_SSMI - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/SSMI/processed/native/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ESACCI/SSMI/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ESACCI_SSMI_r1i1p1_199201-200812_SH.nc

# ICDC_ASI_AMSRE - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-AMSRE/processed/native/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_s.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-AMSRE/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"lon lat" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ICDC_ASI-AMSRE_r1i1p1_200301-201112_s.nc

# ICDC_ASI_SSMI - SH
filein=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-SSMI/processed/native/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc
rootout=${TECLIM_CLIMATE_DATA}/obs/ice/siconc/ICDC/ASI-SSMI/processed/interp/
mkdir -p $rootout
cdo sinfov $filein
ncatted -a coordinates,siconc,c,c,"longitude latitude" $filein
cdo -L remapbil,r360x180 -selname,siconc $filein regridded.nc
ncrcat -F regridded.nc ${rootout}/siconc_SImon_ICDC_SSMI_r1i1p1_199201-201812_ant.nc

