import numpy as np
import calendar
import datetime as dt
import os 
from interimR import humidity_toolbox
import interimR.lib_ioncdf as ioncdf
from interimR.mod_drown import mod_drown as drwn

class ERAinterim_processing():
	''' A Class to perform operation on ERAinterim files 
	such as decumulation, means, unit conversions,... '''

	def __init__(self,dict_input,drown=True):
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

		self.drown = drown
		if self.drown:
			fid_lsm = ioncdf.opennc(self.lsm_file)
			self.lsm = ioncdf.readnc(fid_lsm,'lsm')
			ioncdf.closenc(fid_lsm)
			self.drownstring = 'drowned_'
		else:
			self.drownstring = ''
		
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

	def process_precip_to_daily(self):
		''' processing the precipitation file :
		sum all cumulated fields and divide by a day'''
		precip_out = np.empty((self.ndays,self.ny,self.nx))
		time = np.empty((self.ndays))
		# open file
		fid_precip = ioncdf.opennc(self.file_precip)
		# read coordinates and time
		lon = ioncdf.readnc(fid_precip,'lon')
		lat = ioncdf.readnc(fid_precip,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the decumulation
		cumul_time = self.nsec_per_day 
		for kt in np.arange(0,self.ndays):
			tmp = np.zeros((self.ny,self.nx))
			nvalues = self.nframes_per_day / self.ncumul # number of values to read
			for kc in np.arange(nvalues):
				# frames to read
				# ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
				kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
				tmp = tmp + ioncdf.readnc_oneframe(fid_precip,'TP',kframe-1) # C indexing hence -1
			precip_out[kt,:,:] = tmp.copy() * self.rho_w / cumul_time # units conversion (from m to m/s then to kg/m2/s)
			if self.target_model == 'ROMS':
				precip_out[kt,:,:] = precip_out[kt,::-1,:]
			if self.drown:
				precip_out[kt,:,:] = self.drown_wrapper(precip_out[kt,:,:])
			this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
			time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
		# close file
		ioncdf.closenc(fid_precip)
		# remove negative values
		precip_out[np.where(precip_out < 0.)] = 0.
		# write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_precip
		my_dict['time_dim']       = self.name_time_precip
		my_dict['time_var']       = self.name_time_precip
		my_dict['units']          = 'kg.m-2.s-1'
                my_dict['long name']      = 'Total Precipitation'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = precip_out.min()
                my_dict['var_valid_max']  = precip_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_precip + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,precip_out,my_dict)
		precip_out = None
		return None

        def process_snow_to_daily(self):
                ''' processing the snowfall file :
                sum all cumulated fields and divide by a day'''
                snow_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_snow = ioncdf.opennc(self.file_snow)
                # read coordinates and time
                lon = ioncdf.readnc(fid_snow,'lon')
                lat = ioncdf.readnc(fid_snow,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
                # run the decumulation
                cumul_time = self.nsec_per_day
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + ioncdf.readnc_oneframe(fid_snow,'SF',kframe-1) # C indexing hence -1
                        snow_out[kt,:,:] = tmp.copy() * self.rho_w / cumul_time # units conversion (from m to m/s then to kg/m2/s)
			if self.target_model == 'ROMS':
				snow_out[kt,:,:] = snow_out[kt,::-1,:]
			if self.drown:
				snow_out[kt,:,:] = self.drown_wrapper(snow_out[kt,:,:])
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_snow)
                # remove negative values
                snow_out[np.where(snow_out < 0.)] = 0.
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_snow
		my_dict['time_dim']       = self.name_time_snow
		my_dict['time_var']       = self.name_time_snow
		my_dict['units']          = 'kg.m-2.s-1'
                my_dict['long name']      = 'Snow Fall'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = snow_out.min()
                my_dict['var_valid_max']  = snow_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_snow + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,snow_out,my_dict)
		snow_out = None
                return None

        def process_radlw_to_daily(self):
                ''' processing the longwave radiation file :
                sum all cumulated fields and divide by a day'''
                radlw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radlw = ioncdf.opennc(self.file_radlw)
                # read coordinates and time
                lon = ioncdf.readnc(fid_radlw,'lon')
                lat = ioncdf.readnc(fid_radlw,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
                # run the decumulation
                cumul_time = self.nsec_per_day
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + ioncdf.readnc_oneframe(fid_radlw,'STRD',kframe-1) # C indexing hence -1
                        radlw_out[kt,:,:] = tmp.copy() / cumul_time # units conversion (from W.m-2.s to W.m-2)
			if self.target_model == 'ROMS':
				radlw_out[kt,:,:] = radlw_out[kt,::-1,:]
			if self.drown:
				radlw_out[kt,:,:] = self.drown_wrapper(radlw_out[kt,:,:])
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_radlw)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_radlw
		my_dict['time_dim']       = self.name_time_radlw
		my_dict['time_var']       = self.name_time_radlw
		my_dict['units']          = 'W/m2'
                my_dict['long name']      = 'Downwelling longwave radiation'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = radlw_out.min()
                my_dict['var_valid_max']  = radlw_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_radlw + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,radlw_out,my_dict)
		radlw_out = None
                return None

        def process_radsw_to_daily(self):
                ''' processing the shortwave radiation file :
                sum all 3h fields and divide by a day'''
                radsw_out = np.empty((self.ndays,self.ny,self.nx))
                time = np.empty((self.ndays))
                # open file
                fid_radsw = ioncdf.opennc(self.file_radsw)
                # read coordinates and time
                lon = ioncdf.readnc(fid_radsw,'lon')
                lat = ioncdf.readnc(fid_radsw,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
                # run the summation
                cumul_time = self.nsec_per_day
                for kt in np.arange(0,self.ndays):
                        tmp = np.zeros((self.ny,self.nx))
                        nvalues = self.nframes_per_day / self.ncumul # number of values to read
                        for kc in np.arange(nvalues):
                                # frames to read
                                # ERAinterim (nframes_per_day = 8) : for kt = 0 (day 1) we read frames 4 and 8
                                kframe = (kt * self.nframes_per_day) + (kc+1) * self.ncumul
                                tmp = tmp + ioncdf.readnc_oneframe(fid_radsw,'SSRD',kframe-1) # C indexing hence -1
                        radsw_out[kt,:,:] = tmp.copy() / cumul_time # units conversion (from W.m-2.s to W.m-2)
			if self.target_model == 'ROMS':
				radsw_out[kt,:,:] = radsw_out[kt,::-1,:]
			if self.drown:
				radsw_out[kt,:,:] = self.drown_wrapper(radsw_out[kt,:,:])
                        this_day = dt.datetime(self.year,1,1,12,0) + dt.timedelta(days=int(kt))
                        time[kt] = (this_day - self.reftime).days + (this_day - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_radsw)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_radsw
		my_dict['time_dim']       = self.name_time_radsw
		my_dict['time_var']       = self.name_time_radsw
		my_dict['units']          = 'W/m2'
                my_dict['long name']      = 'Shortwave radiation'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = radsw_out.min()
                my_dict['var_valid_max']  = radsw_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_radsw + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,radsw_out,my_dict)
		radsw_out = None
                return None

	def create_q2_file(self):
		''' Create specific humidity file '''
		q2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_d2 = ioncdf.opennc(self.file_d2)
                fid_msl = ioncdf.opennc(self.file_msl)
                # read coordinates and time
                lon = ioncdf.readnc(fid_d2,'lon')
                lat = ioncdf.readnc(fid_d2,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			d2 = ioncdf.readnc_oneframe(fid_d2,'D2M',kt)
			msl = ioncdf.readnc_oneframe(fid_msl,'MSL',kt)
			q2_tmp = humidity_toolbox.q2_from_d2_and_msl(d2,msl)
			q2_out[kt,:,:] = q2_tmp.copy()
			if self.target_model == 'ROMS':
				q2_out[kt,:,:] = q2_out[kt,::-1,:]
			if self.drown:
				q2_out[kt,:,:] = self.drown_wrapper(q2_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_d2)
                ioncdf.closenc(fid_msl)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_q2
		my_dict['time_dim']       = self.name_time_q2
		my_dict['time_var']       = self.name_time_q2
		my_dict['units']          = 'kg/kg'
                my_dict['long name']      = 'Specific humidity at 2m'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = q2_out.min()
                my_dict['var_valid_max']  = q2_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_q2 + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,q2_out,my_dict)
		q2_out = None
                return None

	def process_t2_file(self):
		''' Rewrite temperature file according to model's needs '''
		t2_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_t2 = ioncdf.opennc(self.file_t2)
                # read coordinates and time
                lon = ioncdf.readnc(fid_t2,'lon')
                lat = ioncdf.readnc(fid_t2,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			t2_out[kt,:,:] = ioncdf.readnc_oneframe(fid_t2,'T2M',kt)
			if self.target_model == 'ROMS':
				t2_out[kt,:,:] = t2_out[kt,::-1,:] - 273.15
			if self.drown:
				t2_out[kt,:,:] = self.drown_wrapper(t2_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_t2)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_t2
		my_dict['time_dim']       = self.name_time_t2
		my_dict['time_var']       = self.name_time_t2
                my_dict['long name']      = 'Air Temperature at 2m'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = t2_out.min()
                my_dict['var_valid_max']  = t2_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
                                            self.drownstring + self.name_t2 + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
	                model_dependent = {'units':'degC','description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
	                model_dependent = {'units':'K','description':my_dict['description'] + '\nNEMO-ready ERAinterim forcing'}
		my_dict.update(model_dependent)
                ioncdf.write_ncfile(lon,lat,time,t2_out,my_dict)
		t2_out = None
                return None

	def process_msl_file(self):
		''' Rewrite pressure file according to model's needs '''
		msl_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_msl = ioncdf.opennc(self.file_msl)
                # read coordinates and time
                lon = ioncdf.readnc(fid_msl,'lon')
                lat = ioncdf.readnc(fid_msl,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			msl_out[kt,:,:] = ioncdf.readnc_oneframe(fid_msl,'MSL',kt)
			if self.target_model == 'ROMS':
				msl_out[kt,:,:] = msl_out[kt,::-1,:]
			if self.drown:
				msl_out[kt,:,:] = self.drown_wrapper(msl_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_msl)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_msl
		my_dict['time_dim']       = self.name_time_msl
		my_dict['time_var']       = self.name_time_msl
		my_dict['units']          = 'Pa'
                my_dict['long name']      = 'Mean Sea level Pressure'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = msl_out.min()
                my_dict['var_valid_max']  = msl_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_msl + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,msl_out,my_dict)
		msl_out = None
		return None

	def process_tcc_file(self):
		''' Rewrite pressure file according to model's needs '''
		tcc_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_tcc = ioncdf.opennc(self.file_tcc)
                # read coordinates and time
                lon = ioncdf.readnc(fid_tcc,'lon')
                lat = ioncdf.readnc(fid_tcc,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			tcc_out[kt,:,:] = ioncdf.readnc_oneframe(fid_tcc,'TCC',kt)
			if self.target_model == 'ROMS':
				tcc_out[kt,:,:] = tcc_out[kt,::-1,:]
			if self.drown:
				tcc_out[kt,:,:] = self.drown_wrapper(tcc_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_tcc)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_tcc
		my_dict['time_dim']       = self.name_time_tcc
		my_dict['time_var']       = self.name_time_tcc
		my_dict['units']          = ''
                my_dict['long name']      = 'Total cloud cover'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = tcc_out.min()
                my_dict['var_valid_max']  = tcc_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_tcc + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,tcc_out,my_dict)
		tcc_out = None
		return None

	def process_u10_file(self):
		''' Rewrite zonal wind file according to model's needs '''
		u10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_u10 = ioncdf.opennc(self.file_u10)
                # read coordinates and time
                lon = ioncdf.readnc(fid_u10,'lon')
                lat = ioncdf.readnc(fid_u10,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			u10_out[kt,:,:] = ioncdf.readnc_oneframe(fid_u10,'U10M',kt)
			if self.target_model == 'ROMS':
				u10_out[kt,:,:] = u10_out[kt,::-1,:]
			if self.drown:
				u10_out[kt,:,:] = self.drown_wrapper(u10_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_u10)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_u10
		my_dict['time_dim']       = self.name_time_u10
		my_dict['time_var']       = self.name_time_u10
		my_dict['units']          = 'm/s'
                my_dict['long name']      = 'Zonal wind speed at 10m'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = u10_out.min()
                my_dict['var_valid_max']  = u10_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_u10 + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,u10_out,my_dict)
		u10_out = None
		return None

	def process_v10_file(self):
		''' Rewrite meridional wind file according to model's needs '''
		v10_out = np.empty((self.nframes,self.ny,self.nx))
                time = np.empty((self.nframes))
                # open file
                fid_v10 = ioncdf.opennc(self.file_v10)
                # read coordinates and time
                lon = ioncdf.readnc(fid_v10,'lon')
                lat = ioncdf.readnc(fid_v10,'lat')
		if self.target_model == 'ROMS':
			lat = lat [::-1]
		# run the computation
		for kt in np.arange(0,self.nframes):
			v10_out[kt,:,:] = ioncdf.readnc_oneframe(fid_v10,'V10M',kt)
			if self.target_model == 'ROMS':
				v10_out[kt,:,:] = v10_out[kt,::-1,:]
			if self.drown:
				v10_out[kt,:,:] = self.drown_wrapper(v10_out[kt,:,:])
			this_time = dt.datetime(self.year,1,1,0,0) + dt.timedelta(seconds=int(kt)*86400/self.nframes_per_day)
                        time[kt] = (this_time - self.reftime).days + (this_time - self.reftime).seconds / 86400.
                # close file
                ioncdf.closenc(fid_v10)
                # write file
		my_dict = {}
		my_dict['description']    = 'file processed by interimR. contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_v10
		my_dict['time_dim']       = self.name_time_v10
		my_dict['time_var']       = self.name_time_v10
		my_dict['units']          = 'm/s'
                my_dict['long name']      = 'Meridional wind speed at 10m'
                my_dict['spval']          = self.spval
                my_dict['reftime']        = self.reftime
                my_dict['time_valid_min'] = time.min()
                my_dict['time_valid_max'] = time.max()
                my_dict['var_valid_min']  = v10_out.min()
                my_dict['var_valid_max']  = v10_out.max()
		my_dict['fileout']        = self.processed_nc_dir + \
		                            self.drownstring + self.name_v10 + '_' + self.dataset + '_' + str(self.year) + '.nc'
		if self.target_model == 'ROMS':
			model_dependent = {'description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
		elif self.target_model == 'NEMO':
			model_dependent = {'description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
		my_dict.update(model_dependent)
		ioncdf.write_ncfile(lon,lat,time,v10_out,my_dict)
		v10_out = None
		return None

	#------------------ compute functions ------------------------------------------

	def decumul(self,data,nchunk):
		''' decumul a chunk of N cumulated values '''
		data_out = data.copy()
		for kt in np.arange(1,nchunk):
			data_out[kt,:,:] = data[kt,:,:] - data[kt-1,:,:]
		return data_out

	#------------------ wrapper drown -----------------------------------------------

	def drown_wrapper(self,data_in):
		''' A wrapper for drown to make it cleaner in the main functions '''
		tmp = drwn.drown(0,data_in[:,:].transpose(),self.lsm.transpose(),200,40)
		data_out = tmp.transpose()
		return data_out


