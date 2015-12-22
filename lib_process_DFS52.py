import netCDF4 as nc
import numpy as np
import calendar
import datetime as dt
import humidity_toolbox
import os 
import lib_ioncdf as ioncdf
import matplotlib.pylab as plt

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
			#self.process_radlw_file()
		if self.dict_input.has_key('file_radsw'):
			print 'Create DFS5.2 shortwave file...'
			#self.process_radsw_file()
		#if self.dict_input.has_key('file_d2'):
		#	print 'Create specific humidity file...'
		#	self.create_q2_file()
		if self.dict_input.has_key('file_t2'):
			print 'Create DFS5.2 t2 file...'
			self.process_t2_file()
		#if self.dict_input.has_key('file_msl'):
		#	print 'Rewrite msl file...'
		#	self.process_msl_file()
		#if self.dict_input.has_key('file_tcc'):
		#	print 'Rewrite tcc file...'
		#	self.process_tcc_file()
		if self.dict_input.has_key('file_u10'):
			print 'Create DFS5.2 u10 file...'
			#self.process_u10_file()
		if self.dict_input.has_key('file_v10'):
			print 'Create DFS5.2 v10 file...'
			#self.process_v10_file()
		return None

	#------------------ Meta functions ------------------------------------------

	def process_radlw_file(self):
		''' Create longwave radiation file '''
		nframes   = self.nframes / self.nframes_per_day # daily file
		radlw_tmp = np.empty((self.ny,self.nx))
		radlw_out = np.empty((nframes,self.ny,self.nx))
                # open input file
                fid_radlw = ioncdf.opennc(self.file_radlw)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_radlw,'lon')
                lat  = ioncdf.readnc(fid_radlw,'lat')
                time = ioncdf.readnc(fid_radlw,self.name_time_radlw)
		# open radlw factor
		fid_factor = ioncdf.opennc(self.longwave_factor)
		radlw_factor = ioncdf.readnc(fid_factor,'radlw')
		ioncdf.closenc(fid_factor)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			radlw_factor = radlw_factor[::-1,:]
		# run the computation
		correct_heat_budget = 1.0088 # attempt to close heat budget 
		for kt in np.arange(0,nframes):
			radlw_tmp[:,:]    = ioncdf.readnc_oneframe(fid_radlw,self.name_radlw,kt)
			radlw_out[kt,:,:] = radlw_tmp[:,:] * radlw_factor[:,:] * correct_heat_budget
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'lwrad_down','time_dim':'lrf_time','time_var':'lrf_time','long name':'Downwelling longwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'radlw','time_dim':'time','time_var':'time','long name':'Downwelling longwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = radlw_out.min()
		my_dict['var_valid_max']  = radlw_out.max()
                # close input file and write output
                ioncdf.closenc(fid_radlw)
	        ioncdf.write_ncfile(lon,lat,time,radlw_out,my_dict)
		# clear arrays
		radlw_tmp = None ; radlw_out = None
		return None

	def process_radsw_file(self):
		''' Create shortwave radiation file '''
		nframes   = self.nframes / self.nframes_per_day # daily file
		radsw_tmp = np.empty((self.ny,self.nx))
		radsw_out = np.empty((nframes,self.ny,self.nx))
                # open file
                fid_radsw = ioncdf.opennc(self.file_radsw)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_radsw,'lon')
                lat  = ioncdf.readnc(fid_radsw,'lat')
                time = ioncdf.readnc(fid_radsw,self.name_time_radsw)
		# open radsw factor
		fid_factor = ioncdf.opennc(self.shortwave_factor)
		radsw_factor = ioncdf.readnc(fid_factor,'radsw')
		ioncdf.closenc(fid_factor)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			radsw_factor = radsw_factor[::-1,:]
		# run the computation
		for kt in np.arange(0,nframes):
			radsw_tmp[:,:]    = ioncdf.readnc_oneframe(fid_radsw,self.name_radsw,kt)
			radsw_out[kt,:,:] = radsw_tmp[:,:] * radsw_factor[:,:]
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'swrad','time_dim':'srf_time','time_var':'srf_time','long name':'Shortwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'radsw','time_dim':'time','time_var':'time','long name':'Shortwave radiation',\
                        'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['var_valid_min']  = radsw_out.min()
		my_dict['var_valid_max']  = radsw_out.max()
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
                # close input file and write output
                ioncdf.closenc(fid_radsw)
	        ioncdf.write_ncfile(lon,lat,time,radsw_out,my_dict)
		radsw_tmp = None ; radsw_out = None
		return None

	def process_u10_file(self):
		''' Create zonal wind file '''
		u10_tmp = np.empty((self.ny,self.nx))
		u10_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_u10 = ioncdf.opennc(self.file_u10)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_u10,'lon')
                lat  = ioncdf.readnc(fid_u10,'lat')
                time = ioncdf.readnc(fid_u10,self.name_time_u10)
		# open background velocity
		fid_bgd = ioncdf.opennc(self.u10_background)
		u10_background = ioncdf.readnc(fid_bgd,'u10')
		ioncdf.closenc(fid_bgd)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			u10_background = u10_background[::-1,:]
		# run the computation
		for kt in np.arange(0,self.nframes):
			u10_tmp[:,:]    = ioncdf.readnc_oneframe(fid_u10,self.name_u10,kt)
			u10_out[kt,:,:] = u10_tmp[:,:] + u10_background[:,:]
                # output file informations
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Uwind','time_dim':'wind_time','time_var':'wind_time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'u10','time_dim':'time','time_var':'time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['var_valid_min']  = u10_out.min()
		my_dict['var_valid_max']  = u10_out.max()
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
                # close input file and write output
                ioncdf.closenc(fid_u10)
	        ioncdf.write_ncfile(lon,lat,time,u10_out,my_dict)
		# clear arrays
		u10_tmp = None ; u10_out = None
		return None

	def process_v10_file(self):
		''' Create meridional wind file '''
		v10_tmp = np.empty((self.ny,self.nx))
		v10_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_v10 = ioncdf.opennc(self.file_v10)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_v10,'lon')
                lat  = ioncdf.readnc(fid_v10,'lat')
                time = ioncdf.readnc(fid_v10,self.name_time_v10)
		# open background velocity
		fid_bgd = ioncdf.opennc(self.v10_background)
		v10_background = ioncdf.readnc(fid_bgd,'v10')
		ioncdf.closenc(fid_bgd)
		# flip upside down for ROMS
		if self.target_model == 'ROMS':
			v10_background = v10_background[::-1,:]
		# run the computation
		for kt in np.arange(0,self.nframes):
			v10_tmp[:,:]    = ioncdf.readnc_oneframe(fid_v10,self.name_v10,kt)
			v10_out[kt,:,:] = v10_tmp[:,:] + v10_background[:,:]
                # output file informations
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Vwind','time_dim':'wind_time','time_var':'wind_time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'v10','time_dim':'time','time_var':'time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['var_valid_min']  = v10_out.min()
		my_dict['var_valid_max']  = v10_out.max()
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
                # close input file and write output
                ioncdf.closenc(fid_v10)
	        ioncdf.write_ncfile(lon,lat,time,v10_out,my_dict)
		# clear arrays
		v10_tmp = None ; v10_out = None
		return None

	def process_t2_file(self):
		''' Create Air temperature file '''
		t2_tmp = np.empty((self.ny,self.nx))
		t2_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_t2 = ioncdf.opennc(self.file_t2)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_t2,'lon')
                lat  = ioncdf.readnc(fid_t2,'lat')
                time = ioncdf.readnc(fid_t2,self.name_time_t2)
	
		#------- Antarctic correction --------
		# make sure that we use latitude from south to north
		if self.target_model == 'ROMS':
			lat_s2n = lat
		else:
			lat_s2n = lat[::-1]

		value_south = -2.0 # we remove 2 degrees C
		lattrans1   = -60 # start linear transition to correction at 60S
		lattrans2   = -75 # correction is at full value at 75S
		jtrans1 = (np.abs(lat_s2n-lattrans1)).argmin()
		jtrans2 = (np.abs(lat_s2n-lattrans2)).argmin()

		correction_south = np.zeros((self.ny,self.nx))
		correction_south[0:jtrans2,:] = value_south
		for jj in np.arange(jtrans2,jtrans1):
			correction_south[jj,:] = value_south * float(jj - jtrans1) / float(jtrans2 - jtrans1)

		if not self.target_model == 'ROMS':
			correction_south = correction_south[::-1,:]

		# run the computation
		for kt in np.arange(0,self.nframes):
			t2_tmp[:,:]    = ioncdf.readnc_oneframe(fid_t2,self.name_t2,kt)
			t2_out[kt,:,:] = t2_tmp[:,:] + correction_south[:,:]










                # output file informations
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Tair','time_dim':'tair_time','time_var':'tair_time','long name':'Air Temperature at 2m',\
	                'units':'degC','fileout':self.output_dir + 't2_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'t2','time_dim':'time','time_var':'time','long name':'Air Temperature at 2m',\
	                'units':'degC','fileout':self.output_dir + 't2_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['var_valid_min']  = t2_out.min()
		my_dict['var_valid_max']  = t2_out.max()
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
                # close input file and write output
                ioncdf.closenc(fid_t2)
	        ioncdf.write_ncfile(lon,lat,time,t2_out,my_dict)
		# clear arrays
		t2_tmp = None ; t2_out = None
		return None
