#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

# You need a copy of NSIDC files interp on ERAinterim grid

#nsidc_dir = '/fsnet/data/meom/DATA_SET/ICE_NSIDC/interp_ERAint/'
nsidc_dir = '/Volumes/P1/Data/NSIDC/interp_ERAint/'

# Build the list of files we need, from 1979 to 1998
filelist = []
for year in np.arange(1979,1998+1):
	filelist.append(nsidc_dir + 'NSIDC_IceConcentration-ERAinterim_y' + str(year) + '.nc')

# Read the data for Ice Concentration
mfid_in = nc.MFDataset(filelist,'r')
icefrac_in = mfid_in.variables['ice_conc'][:]
lon        = mfid_in.variables['lon'][0,:]
lat        = mfid_in.variables['lat'][:,0]
mfid_in.close()

# Read lan sea mask
file_lsm = '../../masks/lsm_ERAinterim_roms.nc'
fid_lsm = nc.Dataset(file_lsm,'r')
lsm = fid_lsm.variables['lsm'][:]
fid_lsm.close()

# Compute monthly mean
nt,ny,nx = icefrac_in.shape
spval = -9999.

icefrac_tmp = np.empty((ny,nx))
icefrac_out = np.empty((12,ny,nx))
for month in np.arange(12):
	icefrac_tmp[:,:] = icefrac_in[month::12,:,:].mean(axis=0)
	icefrac_tmp = icefrac_tmp[::-1,:]       # flip upside down
	icefrac_tmp[:(ny/2),:] = 0.             # remove southern hemisphere
	icefrac_tmp[np.where(lsm == 0)] = spval # mask with lsm
	icefrac_out[month,:,:] = icefrac_tmp[:,:]

lat_out = lat[::-1]
lon_out = lon

days_in_month = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
cumul_days = days_in_month.cumsum() - (days_in_month / 2.)

# write monthly climato for Ice concentration
fid = nc.Dataset('./NSIDC_ice_fraction_monthly_clim_y1979-1998.nc', 'w', format='NETCDF3_CLASSIC')
fid.description = 'Ice Fraction monthly climatology from NSIDC'
# dimensions
fid.createDimension('lat', ny)
fid.createDimension('lon', nx)
fid.createDimension('time', None)
# variables
latitudes  = fid.createVariable('lat',  'f8', ('lat',))
longitudes = fid.createVariable('lon',  'f8', ('lon',))
times      = fid.createVariable('time', 'f8', ('time',))
variable   = fid.createVariable('ifrac','f4', ('time','lat','lon',),fill_value=spval)

# attributes
longitudes.units = "degrees_east"
longitudes.valid_min = lon.min()
longitudes.valid_max = lon.max()
longitudes.long_name = "longitude"

latitudes.units = "degrees_north"
latitudes.valid_min = lat.min()
latitudes.valid_max = lat.max()
latitudes.long_name = "latitude"

times[:] = cumul_days
times.units ="days since 1900-01-01"

variable.long_name = 'Ice Fraction'
variable.coordinates = "lon lat"
variable.missing_value = spval

# data
latitudes[:]    = lat_out
longitudes[:]   = lon_out
variable[:,:,:] = icefrac_out

# close
fid.close()

