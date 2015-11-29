import netCDF4 as nc
import numpy as np
import calendar
import datetime as dt
import humidity_toolbox
import os 

class DFS52_processing():
	''' A Class to create the DFS5.2 dataset from ERAinterim '''

	def __init__(self,dict_input,dict_datafiles):
		# constants
		self.rho_w = 1000.
		self.nsec_per_day = 86400.
		self.reftime = dt.datetime(1900,1,1,0,0)
		self.spval = 1.0e+15
		self.dataset = 'DFS5.2'
		# read inputs
		self.dict_input = dict_input
		for key in dict_input:
                	exec('self.' + key + '=dict_input[key]')
		for key in dict_datafiles:
                	exec('self.' + key + '=dict_datafiles[key]')
		# set time related stuff
		if calendar.isleap(self.year):
			self.ndays = 366
		else:
			self.ndays = 365
		# forcing set freq
		if self.freq == '3h':
			self.nframes_per_day = 8
		elif self.freq == '6h':
			self.nframes_per_day = 4

		self.nframes = self.ndays * self.nframes_per_day
		
		print self.year, 'has', self.nframes, 'frames'

		# variable names
		if self.target_model == 'ROMS':
			self.name_t2     = 'Tair'        ; self.name_time_t2     = 'tair_time'
			self.name_q2     = 'Qair'        ; self.name_time_q2     = 'qair_time'
			self.name_u10    = 'Uwind'       ; self.name_time_u10    = 'wind_time'
			self.name_v10    = 'Vwind'       ; self.name_time_v10    = 'wind_time'
			self.name_radsw  = 'swrad'       ; self.name_time_radsw  = 'srf_time'
			self.name_radlw  = 'lwrad_down'  ; self.name_time_radlw  = 'lrf_time'
			self.name_precip = 'rain'        ; self.name_time_precip = 'rain_time'
			self.name_snow   = 'rain'        ; self.name_time_snow   = 'rain_time'
			self.name_msl    = 'Pair'        ; self.name_time_msl    = 'pair_time'
		else:
			self.name_t2     = 't2'          ; self.name_time_t2     = 'time'
			self.name_q2     = 'q2'          ; self.name_time_q2     = 'time'
			self.name_u10    = 'u10'         ; self.name_time_u10    = 'time'
			self.name_v10    = 'v10'         ; self.name_time_v10    = 'time'
			self.name_radsw  = 'radsw'       ; self.name_time_radsw  = 'time'
			self.name_radlw  = 'radlw'       ; self.name_time_radlw  = 'time'
			self.name_precip = 'precip'      ; self.name_time_precip = 'time'
			self.name_snow   = 'snow'        ; self.name_time_snow   = 'time'
			self.name_msl    = 'msl'         ; self.name_time_msl    = 'time'
			self.name_tcc    = 'tcc'         ; self.name_time_tcc    = 'time'
		return None

	def __call__(self):
		#if self.dict_input.has_key('file_precip'):
		#	print 'Processing precip file...'
		#	self.process_precip_to_daily()
		#if self.dict_input.has_key('file_snow'):
		#	print 'Processing snow file...'
		#	self.process_snow_to_daily()
		if self.dict_input.has_key('file_radlw'):
			print 'Create DFS5.2 longwave file...'
			self.process_radlw_file()
		if self.dict_input.has_key('file_radsw'):
			print 'Processing shortwave file...'
			self.process_radsw_file()
		#if self.dict_input.has_key('file_d2'):
		#	print 'Create specific humidity file...'
		#	self.create_q2_file()
		#if self.dict_input.has_key('file_t2'):
		#	print 'Rewrite t2 file...'
		#	self.process_t2_file()
		#if self.dict_input.has_key('file_msl'):
		#	print 'Rewrite msl file...'
		#	self.process_msl_file()
		#if self.dict_input.has_key('file_tcc'):
		#	print 'Rewrite tcc file...'
		#	self.process_tcc_file()
		if self.dict_input.has_key('file_u10'):
			print 'Create DFS5.2 u10 file...'
			self.process_u10_file()
		if self.dict_input.has_key('file_v10'):
			print 'Create DFS5.2 v10 file...'
			self.process_v10_file()
		return None

	#------------------ Meta functions ------------------------------------------

	def process_radlw_file(self):
		''' Create longwave radiation file '''
		radlw_tmp = np.empty((self.ny,self.nx))
		radlw_min = +1.0e+36
		radlw_max = -1.0e+36
		nframes = self.nframes / self.nframes_per_day # daily file
                # open file
                fid_radlw = self._opennc(self.file_radlw)
                # read coordinates and time
                lon = self._readnc(fid_radlw,'lon')
                lat = self._readnc(fid_radlw,'lat')
                time = self._readnc(fid_radlw,self.name_time_radlw)
		time_min = time.min()
		time_max = time.max()
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'lwrad_down','time_dim':'lrf_time','time_var':'lrf_time','long name':'Downwelling longwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'radlw','time_dim':'time','time_var':'time','long name':'Downwelling longwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '.nc'}
		# faster to loop twice than write attributes after data
		for kt in np.arange(0,nframes):
			radlw_tmp[:,:]   = self._readnc_oneframe(fid_radlw,self.name_radlw,kt)
			# attributes
			radlw_min_tmp = radlw_tmp.min()
			radlw_max_tmp = radlw_tmp.max()
			radlw_min = np.minimum(radlw_min,radlw_min_tmp)
			radlw_max = np.maximum(radlw_max,radlw_max_tmp)

		my_dict['var_valid_min'] = radlw_min
		my_dict['var_valid_max'] = radlw_max
		my_dict['time_valid_min'] = time_min
		my_dict['time_valid_max'] = time_max
                # open output file to write
	        fidout,timeout,varout = self._create_ncfile(lon,lat,time,my_dict)
		# open radlw factor
		fid_factor = self._opennc(self.longwave_factor)
		radlw_factor = self._readnc(fid_factor,'radlw')
		self._closenc(fid_factor)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			radlw_factor = radlw_factor[::-1,:]

		# run the computation
		correct_heat_budget = 1.0088 # attempt to close heat budget 
		for kt in np.arange(0,nframes):
			radlw_tmp[:,:]   = self._readnc_oneframe(fid_radlw,self.name_radlw,kt)
			timeout[kt]    = time[kt]
			varout[kt,:,:] = radlw_tmp[:,:] * radlw_factor[:,:] * correct_heat_budget
                # close file
                self._closenc(fid_radlw)
		self._finalize_ncfile(fidout)
		radlw_tmp = None
		return None

	def process_radsw_file(self):
		''' Create shortwave radiation file '''
		radsw_tmp = np.empty((self.ny,self.nx))
		radsw_min = +1.0e+36
		radsw_max = -1.0e+36
		nframes = self.nframes / self.nframes_per_day # daily file
                # open file
                fid_radsw = self._opennc(self.file_radsw)
                # read coordinates and time
                lon = self._readnc(fid_radsw,'lon')
                lat = self._readnc(fid_radsw,'lat')
                time = self._readnc(fid_radsw,self.name_time_radsw)
		time_min = time.min()
		time_max = time.max()
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'swrad','time_dim':'srf_time','time_var':'srf_time','long name':'Shortwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'radsw','time_dim':'time','time_var':'time','long name':'Shortwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '.nc'}
		# faster to loop twice than write attributes after data
		for kt in np.arange(0,nframes):
			radsw_tmp[:,:]   = self._readnc_oneframe(fid_radsw,self.name_radsw,kt)
			# attributes
			radsw_min_tmp = radsw_tmp.min()
			radsw_max_tmp = radsw_tmp.max()
			radsw_min = np.minimum(radsw_min,radsw_min_tmp)
			radsw_max = np.maximum(radsw_max,radsw_max_tmp)

		my_dict['var_valid_min'] = radsw_min
		my_dict['var_valid_max'] = radsw_max
		my_dict['time_valid_min'] = time_min
		my_dict['time_valid_max'] = time_max
                # open output file to write
	        fidout,timeout,varout = self._create_ncfile(lon,lat,time,my_dict)
		# open radsw factor
		fid_factor = self._opennc(self.shortwave_factor)
		radsw_factor = self._readnc(fid_factor,'radsw')
		self._closenc(fid_factor)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			radsw_factor = radsw_factor[::-1,:]

		# run the computation
		for kt in np.arange(0,nframes):
			radsw_tmp[:,:]   = self._readnc_oneframe(fid_radsw,self.name_radsw,kt)
			timeout[kt]    = time[kt]
			varout[kt,:,:] = radsw_tmp[:,:] * radsw_factor[:,:]
                # close file
                self._closenc(fid_radsw)
		self._finalize_ncfile(fidout)
		radsw_tmp = None
		return None

	def process_u10_file(self):
		''' Create zonal wind file '''
		u10_tmp = np.empty((self.ny,self.nx))
		u10_min = +1.0e+36
		u10_max = -1.0e+36
                # open file
                fid_u10 = self._opennc(self.file_u10)
                # read coordinates and time
                lon = self._readnc(fid_u10,'lon')
                lat = self._readnc(fid_u10,'lat')
                time = self._readnc(fid_u10,self.name_time_u10)
		time_min = time.min()
		time_max = time.max()
                # output file informations
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Uwind','time_dim':'wind_time','time_var':'wind_time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'u10','time_dim':'time','time_var':'time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '.nc'}
		# faster to loop twice than write attributes after data
		for kt in np.arange(0,self.nframes):
			u10_tmp[:,:]   = self._readnc_oneframe(fid_u10,self.name_u10,kt)
			# attributes
			u10_min_tmp = u10_tmp.min()
			u10_max_tmp = u10_tmp.max()
			u10_min = np.minimum(u10_min,u10_min_tmp)
			u10_max = np.maximum(u10_max,u10_max_tmp)

		my_dict['var_valid_min'] = u10_min
		my_dict['var_valid_max'] = u10_max
		my_dict['time_valid_min'] = time_min
		my_dict['time_valid_max'] = time_max
                # open output file to write
	        fidout,timeout,varout = self._create_ncfile(lon,lat,time,my_dict)
		# open background velocity
		fid_bgd = self._opennc(self.u10_background)
		u10_background = self._readnc(fid_bgd,'u10')
		self._closenc(fid_bgd)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			u10_background = u10_background[::-1,:]

		# run the computation
		for kt in np.arange(0,self.nframes):
			u10_tmp[:,:]   = self._readnc_oneframe(fid_u10,self.name_u10,kt)
			timeout[kt]    = time[kt]
			varout[kt,:,:] = u10_tmp[:,:] + u10_background[:,:]
                # close file
                self._closenc(fid_u10)
		self._finalize_ncfile(fidout)
		u10_tmp = None
		return None

	def process_v10_file(self):
		''' Create meridional wind file '''
		v10_tmp = np.empty((self.ny,self.nx))
		v10_min = +1.0e+36
		v10_max = -1.0e+36
                # open file
                fid_v10 = self._opennc(self.file_v10)
                # read coordinates and time
                lon = self._readnc(fid_v10,'lon')
                lat = self._readnc(fid_v10,'lat')
                time = self._readnc(fid_v10,self.name_time_v10)
		time_min = time.min()
		time_max = time.max()
                # output file informations
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'vwind','time_dim':'wind_time','time_var':'wind_time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'v10','time_dim':'time','time_var':'time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '.nc'}
		# faster to loop twice than write attributes after data
		for kt in np.arange(0,self.nframes):
			v10_tmp[:,:]   = self._readnc_oneframe(fid_v10,self.name_v10,kt)
			# attributes
			v10_min_tmp = v10_tmp.min()
			v10_max_tmp = v10_tmp.max()
			v10_min = np.minimum(v10_min,v10_min_tmp)
			v10_max = np.maximum(v10_max,v10_max_tmp)

		my_dict['var_valid_min'] = v10_min
		my_dict['var_valid_max'] = v10_max
		my_dict['time_valid_min'] = time_min
		my_dict['time_valid_max'] = time_max
                # open output file to write
	        fidout,timeout,varout = self._create_ncfile(lon,lat,time,my_dict)
		# open background velocity
		fid_bgd = self._opennc(self.v10_background)
		v10_background = self._readnc(fid_bgd,'v10')
		self._closenc(fid_bgd)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			v10_background = v10_background[::-1,:]

		# run the computation
		for kt in np.arange(0,self.nframes):
			v10_tmp[:,:]   = self._readnc_oneframe(fid_v10,self.name_v10,kt)
			timeout[kt]    = time[kt]
			varout[kt,:,:] = v10_tmp[:,:] + v10_background[:,:]
                # close file
                self._closenc(fid_v10)
		self._finalize_ncfile(fidout)
		v10_tmp = None
		return None









	#------------------ NetCDF functions ------------------------------------------
	def _opennc(self,myfile):
		fid = nc.Dataset(myfile,'r')
		return fid

	def _closenc(self,fid):
		fid.close()
		return None

        def _readnc_frames(self,fid,myvar,myframe_start,myframe_end):
                ''' read data from netcdf '''
                out = fid.variables[myvar][myframe_start:myframe_end,:,:].squeeze()
                return out

        def _readnc_oneframe(self,fid,myvar,myframe):
                ''' read data from netcdf '''
                out = fid.variables[myvar][myframe,:,:].squeeze()
                return out

        def _readnc(self,fid,myvar):
                ''' read data from netcdf '''
                out = fid.variables[myvar][:].squeeze()
                return out

	def _create_ncfile(self,lon_array,lat_array,time,dict_wrt):
        	fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
	        fid.description = 'DFS5.2 (raphael.dussin@gmail.com)'
	        # dimensions
	        fid.createDimension('lat', lat_array.shape[0])
	        fid.createDimension('lon', lon_array.shape[0])
	        fid.createDimension(dict_wrt['time_dim'], None)
	        # variables
	        latitudes  = fid.createVariable('lat', 'f8', ('lat',))
	        longitudes = fid.createVariable('lon', 'f8', ('lon',))
	        times      = fid.createVariable(dict_wrt['time_var'], 'f8', (dict_wrt['time_dim'],))
	        variable   = fid.createVariable(dict_wrt['varname'], 'f4', (dict_wrt['time_dim'],'lat','lon',),fill_value=self.spval)
	        # data
	        latitudes[:]    = lat_array
	        longitudes[:]   = lon_array
	
		# attributes
		longitudes.units = "degrees_east" 
		longitudes.valid_min = lon_array.min()
		longitudes.valid_max = lon_array.max()
		longitudes.long_name = "longitude" 

		latitudes.units = "degrees_north" 
		latitudes.valid_min = lat_array.min()
		latitudes.valid_max = lat_array.max()
		latitudes.long_name = "latitude" 

		times.units = "days since " + self.reftime.isoformat()
                times.calendar = "LEAP"
		times.valid_min = dict_wrt['time_valid_min']
                times.valid_max = dict_wrt['time_valid_max']

		variable.long_name = dict_wrt['long name']
		variable.units = dict_wrt['units']
		variable.coordinates = "lon lat" 
		variable.time = dict_wrt['time_var']
		variable.missing_value = self.spval
		variable.valid_range = dict_wrt['var_valid_min'], dict_wrt['var_valid_max']
		return fid, times, variable

	def _finalize_ncfile(self,fid):
	        # close
	        fid.close()
	        return None
