#!/usr/bin/env python

import netCDF4 as nc 
import numpy as np
import lib_brodeau as lolo
import drown_kara as drown

model = 'ROMS'

if model == 'ROMS':
	file_lsm = '../../masks/lsm_ERAinterim_roms.nc'
	var_name_t2 = 'Tair'
else:
	file_lsm = '../../masks/lsm_ERAinterim.nc'
	var_name_t2 = 't2'

file_t2_ERAinterim = '../../data/DFS5.2/temperature/t2_ERAinterim_monthly_clim_y1979-1998.nc'
file_t2_POLES = '../../data/DFS5.2/temperature/t2_POLES-ERAinterim_monthly_clim_1979-1998.nc'

# calendar stuff
days_in_month          = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
days_in_month_leapyear = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
cumul_days = days_in_month.cumsum() - (days_in_month / 2.)

spval = -9999.

# Read data from ERAinterim climato
fid_ERAint = nc.Dataset(file_t2_ERAinterim,'r')
lon = fid_ERAint.variables['lon'][:]
lat = fid_ERAint.variables['lat'][:]
t2_ERAint = fid_ERAint.variables[var_name_t2][:]
fid_ERAint.close()

# Read data from POLES climato
fid_POLES = nc.Dataset(file_t2_POLES,'r')
t2_POLES = fid_POLES.variables['t2'][:]
fid_POLES.close()

# Read data for lsm
fid_lsm = nc.Dataset(file_lsm,'r')
lsm = fid_lsm.variables['lsm'][:]
fid_lsm.close()

nt, ny, nx = t2_POLES.shape

idx_70N = (np.abs(lat - 70)).argmin()
idx_65N = (np.abs(lat - 65)).argmin()

print idx_70N, lat[idx_70N]
print idx_65N, lat[idx_65N]
#
offset_orig = np.empty((ny,nx))
offset_tmp  = np.empty((ny,nx))

offset_out = np.empty((nt,ny,nx))

for month in np.arange(nt):
	# compute original offset
	offset_orig[:,:] = t2_POLES[month,:,:] - t2_ERAint[month,:,:] 
	# make a working copy
	offset_tmp[:,:] = offset_orig[:,:].copy()
	# drown
	tmp = drown.cslf(offset_tmp.transpose(),lsm.transpose())
	offset_tmp = tmp.transpose()
	# smooth the resulting field
	for nsmooth in np.arange(8):
		tmp2 = lolo.smoother(lsm.transpose(),offset_tmp.transpose())
		offset_tmp = tmp2.transpose()

	# remove everything south of 70N
	offset_tmp[:idx_70N,:] = 0.

	#linear blend
	for jj in np.arange(idx_65N,idx_70N):
		wgt_era = float(idx_70N -jj) / float(idx_70N - idx_65N)
		wgt_poles = 1. - wgt_era
		offset_tmp[jj,:] = offset_orig[jj,:] * wgt_poles
		#offset_tmp[month,jj,:] = offset[month,jj,:] / 1.6

	# mask land
	offset_tmp[np.where(lsm == 0)] = spval

	offset_out[month,:,:] = offset_tmp[:,:] 


# remove extreme positive values
offset_out[np.where(offset_out > 2.)] = 2.0

# write monthly climato for t2 from ERAinterim
fid = nc.Dataset('../../data/DFS5.2/temperature/monthly_offset.nc', 'w', format='NETCDF3_CLASSIC')
fid.description = 'offset POLES - ERAinterim'
# dimensions
fid.createDimension('lat', ny)
fid.createDimension('lon', nx)
fid.createDimension('time', None)
# variables
latitudes  = fid.createVariable('lat',  'f8', ('lat',))
longitudes = fid.createVariable('lon',  'f8', ('lon',))
times      = fid.createVariable('time', 'f8', ('time',))
variable   = fid.createVariable(var_name_t2,   'f4', ('time','lat','lon',),fill_value=spval)

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

variable.long_name = 'offset POLES - ERAinterim'
variable.coordinates = "lon lat"
variable.missing_value = spval

# data
latitudes[:]    = lat
longitudes[:]   = lon
variable[:,:,:] = offset_out

# close
fid.close()

