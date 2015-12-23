#!/usr/bin/env python

import netCDF4 as nc 
import numpy as np

# read input data
fid_in = nc.Dataset('./satiabppolesclim.nc','r')
lon = fid_in.variables['lon'][:]
lat = fid_in.variables['lat'][:]
data = fid_in.variables['data'][0,:,:]
fid_in.close()

ny, nx = data.shape
mask = np.ones((ny,nx))
mask[np.where(data.mask == True)] = 0

# write mask
fid = nc.Dataset('../../data/masks/lsm_POLES.nc', 'w', format='NETCDF3_CLASSIC')
fid.description = 'Land Sea Mask for POLES'
# dimensions
fid.createDimension('lat', ny)
fid.createDimension('lon', nx)
# variables
latitudes  = fid.createVariable('lat', 'f8', ('lat',))
longitudes = fid.createVariable('lon', 'f8', ('lon',))
variable   = fid.createVariable('lsm', 'i4', ('lat','lon',))

# attributes
longitudes.units = "degrees_east"
longitudes.valid_min = lon.min()
longitudes.valid_max = lon.max()
longitudes.long_name = "longitude"

latitudes.units = "degrees_north"
latitudes.valid_min = lat.min()
latitudes.valid_max = lat.max()
latitudes.long_name = "latitude"

variable.long_name = 'Land Sea Mask'
variable.coordinates = "lon lat"

# data
latitudes[:]    = lat
longitudes[:]   = lon
variable[:,:] = mask

# close
fid.close()
