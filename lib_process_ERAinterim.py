import netCDF4 as nc
import numpy as np
import calendar
import datetime as dt
import humidity_toolbox
import os 

class ERAinterim_processing():
	''' A Class to perform operation on ERAinterim files 
	such as decumulation, means, unit conversions,... '''

	def __init__(self,dict_input):
		# constants
		self.rho_w = 1000.
		self.nsec_per_day = 86400.
		self.reftime = dt.datetime(1900,1,1,0,0)
		self.spval = 1.0e+15
		self.dataset = 'ERAinterim'
		# read inputs
		self.dict_input = dict_input
		for key in dict_input:
                	exec('self.' + key + '=dict_input[key]')
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

		return None

	def __call__(self):
		if self.dict_input.has_key('file_precip'):
			print 'Processing precip file...'
			self.process_precip_to_daily()
		if self.dict_input.has_key('file_snow'):
			print 'Processing snow file...'
			self.process_snow_to_daily()
		if self.dict_input.has_key('file_radlw'):
			print 'Processing longwave file...'
			self.process_radlw_to_daily()
		if self.dict_input.has_key('file_radsw'):
			print 'Processing shortwave file...'
			self.process_radsw_to_daily()
		if self.dict_input.has_key('file_d2'):
			print 'Create specific humidity file...'
			self.create_q2_file()
		if self.dict_input.has_key('file_t2'):
			print 'Rewrite t2 file...'
			self.process_t2_file()
		if self.dict_input.has_key('file_msl'):
			print 'Rewrite msl file...'
			self.process_msl_file()
		if self.dict_input.has_key('file_tcc'):
			print 'Rewrite tcc file...'
			self.process_tcc_file()
		if self.dict_input.has_key('file_u10'):
			print 'Rewrite u10 file...'
			self.process_u10_file()
		if self.dict_input.has_key('file_v10'):
			print 'Rewrite v10 file...'
			self.process_v10_file()
		return None

	#------------------ Meta functions ------------------------------------------

	def process_precip(self):
		''' processing the precipitation file :
		decumulating fields and output on native frequency '''
		precip_out = np.empty((self.nframes,self.ny,self.nx))
		# open file
		fid_precip = self._opennc(self.file_precip)
		# read coordinates and time
		lon = self._readnc(fid_precip,'lon')
		lat = self._readnc(fid_precip,'lat')
		time = self._readnc(fid_precip,'time')
		# run the decumulation
		for kt in np.arange(0,self.nframes,self.ncumul):
			tmp = self._readnc_frames(fid_precip,'TP',kt,kt+self.ncumul)
			tmp_decumul = self.decumul(tmp,self.ncumul)
			precip_out[kt:kt+self.ncumul,:,:] = tmp_decumul[:,:,:].copy()
		# close file
		self._closenc(fid_precip)
		# units conversion (from m to m/s then to kg/m2/s)
		cumul_time = self.nsec_per_day * self.ncumul / self.nframes_per_day 
		precip_out = precip_out * self.rho_w / cumul_time
		# remove negative values
		precip_out[np.where(precip_out < 0.)] = 0.
		# write file
		if self.target_model == 'ROMS':
			my_dict = {'varname':'rain','time_dim':'rain_time','time_var':'rain_time','long name':'Total Precipitation',\
			'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
			self._write_ncfile(lon,lat,time,precip_out,my_dict)
		elif self.target_model == 'NEMO':
			my_dict = {'varname':'precip','time_dim':'time','time_var':'time','long name':'Total Precipitation',\
			'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '.nc'}
			self._write_ncfile(lon,lat,time,precip_out,my_dict)
		precip_out = None
		return None
		
	def process_precip_to_daily(self):
		''' processing the precipitation file :
		sum all cumulated fields and divide by a day'''
		precip_out = np.empty((self.ndays,self.ny,self.nx))
		time = np.empty((self.ndays))
		# open file
		fid_precip = self._opennc(self.file_precip)
		# read coordinates and time
		lon = self._readnc(fid_precip,'lon')
		lat = self._readnc(fid_precip,'lat')
		# run the decumulation
		for kt in np.arange(0,self.ndays):
			tmp = np.zeros((self.ny,self.nx))
			nvalues = self.nframes_per_day / self.ncumul # number of values to read
			for kc in np.arange(nvalues):
				# frames to read
				# ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
				kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
				tmp = tmp + self._readnc_oneframe(fid_precip,'TP',kframe-1) # C indexing hence -1

			precip_out[kt,:,:] = tmp.copy()
			this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
			time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
		# close file
		self._closenc(fid_precip)
		# units conversion (from m to m/s then to kg/m2/s)
		cumul_time = self.nsec_per_day 
		precip_out = precip_out * self.rho_w / cumul_time
		# remove negative values
		precip_out[np.where(precip_out < 0.)] = 0.
		# write file
		if self.target_model == 'ROMS':
			my_dict = {'varname':'rain','time_dim':'rain_time','time_var':'rain_time','long name':'Total Precipitation',\
			'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '_daily_ROMS.nc'}
			self._write_ncfile(lon,lat[::-1],time,precip_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
			my_dict = {'varname':'precip','time_dim':'time','time_var':'time','long name':'Total Precipitation',\
			'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '_daily.nc'}
			self._write_ncfile(lon,lat,time,precip_out,my_dict)
		precip_out = None
		return None

        def process_snow_to_daily(self):
                ''' processing the snowfall file :
                sum all cumulated fields and divide by a day'''
                snow_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_snow = self._opennc(self.file_snow)
                # read coordinates and time
                lon = self._readnc(fid_snow,'lon')
                lat = self._readnc(fid_snow,'lat')
                # run the decumulation
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + self._readnc_oneframe(fid_snow,'SF',kframe-1) # C indexing hence -1
                        snow_out[kt,:,:] = tmp.copy()
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_snow)
                # units conversion (from m to m/s then to kg/m2/s)
                cumul_time = self.nsec_per_day
                snow_out = snow_out * self.rho_w / cumul_time
                # remove negative values
                snow_out[np.where(snow_out < 0.)] = 0.
                # write file
		if self.target_model == 'ROMS':
                	my_dict = {'varname':'rain','time_dim':'rain_time','time_var':'rain_time','long name':'Snow Fall',\
 	               'units':'kg.m-2.s-1','fileout':self.output_dir + 'snow_' + self.dataset + '_' + str(self.year) + '_daily_ROMS.nc'}
			self._write_ncfile(lon,lat[::-1],time,snow_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
                	my_dict = {'varname':'snow','time_dim':'time','time_var':'time','long name':'Snow Fall',\
 	               'units':'kg.m-2.s-1','fileout':self.output_dir + 'snow_' + self.dataset + '_' + str(self.year) + '_daily.nc'}
			self._write_ncfile(lon,lat,time,snow_out,my_dict)
		snow_out = None
                return None

        def process_radlw_to_daily(self):
                ''' processing the longwave radiation file :
                sum all cumulated fields and divide by a day'''
                radlw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radlw = self._opennc(self.file_radlw)
                # read coordinates and time
                lon = self._readnc(fid_radlw,'lon')
                lat = self._readnc(fid_radlw,'lat')
                # run the decumulation
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + self._readnc_oneframe(fid_radlw,'STRD',kframe-1) # C indexing hence -1
                        radlw_out[kt,:,:] = tmp.copy()
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_radlw)
                # units conversion (from W.m-2.s to W.m-2)
                cumul_time = self.nsec_per_day
                radlw_out = radlw_out / cumul_time
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'lwrad_down','time_dim':'lrf_time','time_var':'lrf_time','long name':'Downwelling longwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '_daily_ROMS.nc'}
       		        self._write_ncfile(lon,lat[::-1],time,radlw_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'radlw','time_dim':'time','time_var':'time','long name':'Downwelling longwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radlw_' + self.dataset + '_' + str(self.year) + '_daily.nc'}
       		        self._write_ncfile(lon,lat,time,radlw_out,my_dict)
		radlw_out = None
                return None

        def process_radsw_to_daily(self):
                ''' processing the shortwave radiation file :
                sum all 3h fields and divide by a day'''
                radsw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radsw = self._opennc(self.file_radsw)
                # read coordinates and time
                lon = self._readnc(fid_radsw,'lon')
                lat = self._readnc(fid_radsw,'lat')
                # run the summation
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + self._readnc_oneframe(fid_radsw,'SSRD',kframe-1) # C indexing hence -1
                        radsw_out[kt,:,:] = tmp.copy()
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_radsw)
                # units conversion (from W.m-2.s to W.m-2)
                cumul_time = self.nsec_per_day
                radsw_out = radsw_out / cumul_time
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'swrad','time_dim':'srf_time','time_var':'srf_time','long name':'Shortwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '_daily_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,radsw_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'radsw','time_dim':'time','time_var':'time','long name':'Shortwave radiation',\
	                'units':'W.m-2','fileout':self.output_dir + 'radsw_' + self.dataset + '_' + str(self.year) + '_daily.nc'}
	                self._write_ncfile(lon,lat,time,radsw_out,my_dict)
		radsw_out = None
                return None

	def create_q2_file(self):
		''' Create specific humidity file '''
		q2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_d2 = self._opennc(self.file_d2)
                fid_msl = self._opennc(self.file_msl)
                # read coordinates and time
                lon = self._readnc(fid_d2,'lon')
                lat = self._readnc(fid_d2,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			d2 = self._readnc_oneframe(fid_d2,'D2M',kt)
			msl = self._readnc_oneframe(fid_msl,'MSL',kt)
			q2_tmp = humidity_toolbox.q2_from_d2_and_msl(d2,msl)
			q2_out[kt,:,:] = q2_tmp.copy()
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_d2)
                self._closenc(fid_msl)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Qair','time_dim':'qair_time','time_var':'qair_time','long name':'Specific humidity at 2m',\
	                'units':'kg/kg','fileout':self.output_dir + 'q2_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,q2_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'q2','time_dim':'time','time_var':'time','long name':'Specific humidity at 2m',\
	                'units':'kg/kg','fileout':self.output_dir + 'q2_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,q2_out,my_dict)
		q2_out = None
                return None

	def process_t2_file(self):
		''' Rewrite temperature file according to model's needs '''
		t2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_t2 = self._opennc(self.file_t2)
                # read coordinates and time
                lon = self._readnc(fid_t2,'lon')
                lat = self._readnc(fid_t2,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			t2_out[kt,:,:] = self._readnc_oneframe(fid_t2,'T2M',kt) - 273.15
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_t2)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Tair','time_dim':'tair_time','time_var':'tair_time','long name':'Air Temperature at 2m',\
	                'units':'degC','fileout':self.output_dir + 't2_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,t2_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'t2','time_dim':'time','time_var':'time','long name':'Air Temperature at 2m',\
	                'units':'degC','fileout':self.output_dir + 't2_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,t2_out,my_dict)
		t2_out = None
                return None

	def process_msl_file(self):
		''' Rewrite pressure file according to model's needs '''
		msl_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_msl = self._opennc(self.file_msl)
                # read coordinates and time
                lon = self._readnc(fid_msl,'lon')
                lat = self._readnc(fid_msl,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			msl_out[kt,:,:] = self._readnc_oneframe(fid_msl,'MSL',kt)
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_msl)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Pair','time_dim':'pair_time','time_var':'pair_time','long name':'Mean sea-level pressure',\
	                'units':'Pa','fileout':self.output_dir + 'msl_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,msl_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'msl','time_dim':'time','time_var':'time','long name':'Mean sea-level pressure',\
	                'units':'Pa','fileout':self.output_dir + 'msl_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,msl_out,my_dict)
		msl_out = None
		return None

	def process_tcc_file(self):
		''' Rewrite pressure file according to model's needs '''
		tcc_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_tcc = self._opennc(self.file_tcc)
                # read coordinates and time
                lon = self._readnc(fid_tcc,'lon')
                lat = self._readnc(fid_tcc,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			tcc_out[kt,:,:] = self._readnc_oneframe(fid_tcc,'MSL',kt)
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_tcc)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'cloud','time_dim':'cloud_time','time_var':'cloud_time','long name':'Total cloud cover',\
	                'units':'N/A','fileout':self.output_dir + 'tcc_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,tcc_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'tcc','time_dim':'time','time_var':'time','long name':'Total cloud cover',\
	                'units':'N/A','fileout':self.output_dir + 'tcc_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,tcc_out,my_dict)
		tcc_out = None
		return None

	def process_u10_file(self):
		''' Rewrite zonal wind file according to model's needs '''
		u10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_u10 = self._opennc(self.file_u10)
                # read coordinates and time
                lon = self._readnc(fid_u10,'lon')
                lat = self._readnc(fid_u10,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			u10_out[kt,:,:] = self._readnc_oneframe(fid_u10,'U10M',kt)
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_u10)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Uwind','time_dim':'wind_time','time_var':'wind_time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,u10_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'u10','time_dim':'time','time_var':'time','long name':'Zonal wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'u10_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,u10_out,my_dict)
		u10_out = None
		return None

	def process_v10_file(self):
		''' Rewrite meridional wind file according to model's needs '''
		v10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_v10 = self._opennc(self.file_v10)
                # read coordinates and time
                lon = self._readnc(fid_v10,'lon')
                lat = self._readnc(fid_v10,'lat')
		# run the computation
		for kt in np.arange(0,self.nframes):
			v10_out[kt,:,:] = self._readnc_oneframe(fid_v10,'V10M',kt)
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                self._closenc(fid_v10)
                # write file
		if self.target_model == 'ROMS':
	                my_dict = {'varname':'Vwind','time_dim':'wind_time','time_var':'wind_time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
	                self._write_ncfile(lon,lat[::-1],time,v10_out[:,::-1,:],my_dict)
		elif self.target_model == 'NEMO':
	                my_dict = {'varname':'v10','time_dim':'time','time_var':'time','long name':'Meridional wind speed at 10m',\
	                'units':'m/s','fileout':self.output_dir + 'v10_' + self.dataset + '_' + str(self.year) + '.nc'}
	                self._write_ncfile(lon,lat,time,v10_out,my_dict)
		v10_out = None
		return None

	#------------------ compute functions ------------------------------------------

	def decumul(self,data,nchunk):
		''' decumul a chunk of N cumulated values '''
		data_out = data.copy()
		for kt in np.arange(1,nchunk):
			data_out[kt,:,:] = data[kt,:,:] - data[kt-1,:,:]
		return data_out

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

	def _write_ncfile(self,lon_array,lat_array,time,var,dict_wrt):
        	fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
	        fid.description = 'ERAinterim post-processing (raphael.dussin@gmail.com)'
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
	        times[:]        = time
	        variable[:,:,:] = var
	
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
		times.valid_min = time.min()
                times.valid_max = time.max()
                times.calendar = "LEAP"

		variable.long_name = dict_wrt['long name']
		variable.units = dict_wrt['units']
		variable.coordinates = "lon lat" 
		variable.time = dict_wrt['time_var']
		variable.missing_value = self.spval
		variable.valid_range = var.min() , var.max()

	        # close
	        fid.close()
	        return None

class ERAinterim_drown():

        def __init__(self,dict_input):
                # read inputs
                for key in dict_input:
                        exec('self.' + key + '=dict_input[key]')

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
			self.drownexe    = 'mask_drown_field_roms.x'
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
			self.drownexe    = 'mask_drown_field.x'
                return None

	def __call__(self):

		for var in self.listvar:
			exec('filein = self.file_' + var)
			exec('varin = self.name_' + var)
			exec('timevar = self.name_time_' + var)
			command = self.create_drown_command(filein,varin,timevar)
			os.system(command)
		return None

	def _fileout_name(self,filein):
		filetmp = filein.replace('/',' ').split()[-1]
		fileout = 'drowned_' + filetmp
		return fileout

	def create_drown_command(self,filein,varin,timevar):
		fileout = self._fileout_name(filein)
		if self.target_model == 'ROMS':
			command = self.sosie_dir + self.drownexe + ' -D -i ' + filein + ' -v ' + varin + ' -t ' + timevar + \
			' -w ' + timevar + ' -m ' + self.lsm_file + ' -o ' + self.output_dir + fileout
		else:
			command = self.sosie_dir + self.drownexe + ' -D -i ' + filein + ' -v ' + varin + ' -t ' + timevar + \
			' -m ' + self.lsm_file + ' -o ' + self.output_dir + fileout
		return command

