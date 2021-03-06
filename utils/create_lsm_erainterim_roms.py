#!/usr/bin/env python

import netCDF4 as nc

def readnc(filein,varin):
	fid = nc.Dataset(filein,'r')
	out = fid.variables[varin][:]
	fid.close()
	return out

def write_lsm(lon_array,lat_array,var):
                fid = nc.Dataset('../data/masks/lsm_erainterim_roms.nc','w', format='NETCDF3_CLASSIC')
                fid.description = 'ERAinterim post-processing (raphael.dussin@gmail.com)'
                # dimensions
                fid.createDimension('lat', lat_array.shape[0])
                fid.createDimension('lon', lon_array.shape[0])
                # variables
                latitudes  = fid.createVariable('lat', 'f8', ('lat',))
                longitudes = fid.createVariable('lon', 'f8', ('lon',))
                variable   = fid.createVariable('lsm', 'i4', ('lat','lon',))
                # data
                latitudes[:]    = lat_array
                longitudes[:]   = lon_array
                variable[:,:] = var

                # attributes
                longitudes.units = "degrees_east"
                longitudes.valid_min = lon_array.min()
                longitudes.valid_max = lon_array.max()
                longitudes.long_name = "longitude"

                latitudes.units = "degrees_north"
                latitudes.valid_min = lat_array.min()
                latitudes.valid_max = lat_array.max()
                latitudes.long_name = "latitude"

                variable.long_name = 'land sea mask'
                variable.coordinates = "lon lat"
                variable.valid_range = var.min() , var.max()

                # close
                fid.close()
                return None



filelsm_nemo = '../data/masks/lsm_erainterim.nc'

lon = readnc(filelsm_nemo,'lon0')
lat = readnc(filelsm_nemo,'lat0')
lsm = readnc(filelsm_nemo,'lsm')

write_lsm(lon,lat[::-1],lsm[::-1])
