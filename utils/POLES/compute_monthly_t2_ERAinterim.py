#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import calendar

# Data from ERAinterim
model = 'ROMS'

if model == 'ROMS':
	dir_ERAinterim = '/Volumes/P1/Data/ERAinterim/'
	file_lsm = '../../data/masks/lsm_ERAinterim_roms.nc'
	fileroot_ERAinterim = dir_ERAinterim + 'drowned_t2_ERAinterim_<YYYY>_ROMS.nc'
	var_name_t2 = 'Tair'
else:
	dir_ERAinterim = ''
	file_lsm = '../../data/masks/lsm_ERAinterim.nc'
	fileroot_ERAinterim = dir_ERAinterim + 'drowned_t2_ERAinterim_<YYYY>.nc'
	var_name_t2 = 't2'

# years 1979 to 1998
first_year = 1979
last_year  = 1998

# calendar stuff
days_in_month          = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
days_in_month_leapyear = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
cumul_days = days_in_month.cumsum() - (days_in_month / 2.)

# Read coordinates
fid_lsm = nc.Dataset(file_lsm,'r')
lsm = fid_lsm.variables['lsm'][:]
lon = fid_lsm.variables['lon'][:]
lat = fid_lsm.variables['lat'][:]
fid_lsm.close()

ny, nx = lsm.shape
nt = 12 ; nfpd = 8
spval= -9999.

# Alloc arrays
t2_monthly_sum  = np.zeros((nt,ny,nx))
t2_monthly_mean = np.zeros((nt,ny,nx))
nframes_sum = np.zeros((nt))

# Add all values for each month
for year in np.arange(first_year,last_year+1):
	thisyear_file = fileroot_ERAinterim.replace('<YYYY>',str(year))
	fid_in = nc.Dataset(thisyear_file,'r')
	# deal with leap years
	if calendar.isleap(year):
		ndays = days_in_month_leapyear
	else:
		ndays = days_in_month
	# loop on month
	for month in np.arange(12):
		ffm = ndays[:month].sum() * nfpd      # first frame of month
		lfm = ndays[:month+1].sum() * nfpd    # last frame of month
		nframes_sum[month] = nframes_sum[month] + (lfm - ffm)
		t2_monthly_sum[month,:,:] = t2_monthly_sum[month,:,:] + fid_in.variables[var_name_t2][ffm:lfm,:,:].sum(axis=0)
	fid_in.close()

# Compute mean from sum
for month in np.arange(12):
	t2_monthly_mean[month,:,:] = t2_monthly_sum[month,:,:] / nframes_sum[month]

# write monthly climato for t2 from ERAinterim
fid = nc.Dataset('../../data/DFS5.2/temperature/t2_ERAinterim_monthly_clim_y1979-1998.nc', 'w', format='NETCDF3_CLASSIC')
fid.description = 'Air temperature monthly climatology from ERAinterim'
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

variable.long_name = 'Air Temperature at 2m'
variable.coordinates = "lon lat"
variable.missing_value = spval

# data
latitudes[:]    = lat
longitudes[:]   = lon
variable[:,:,:] = t2_monthly_mean

# close
fid.close()
