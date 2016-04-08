import numpy as np
import calendar
import datetime as dt
import os 
from interimR import humidity_toolbox
import interimR.lib_ioncdf as ioncdf
from interimR.mod_drown import mod_drown as drwn

class DFS52_processing():
	''' A Class to create the DFS5.2 dataset from ERAinterim '''

	def __init__(self,dict_input,dict_datafiles,drown=True):
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
			self.days_in_month = np.array([31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
		else:
			self.ndays = 365
			self.days_in_month = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
		# forcing set freq
		if self.freq == '3h':
			self.nframes_per_day = 8
		elif self.freq == '6h':
			self.nframes_per_day = 4

		self.nframes = self.ndays * self.nframes_per_day
		
		print self.year, 'has', self.nframes, 'frames'

		self.drown = drown
		if self.drown:
			self.drownstring = 'drowned_'
		else:
			self.drownstring = ''

		fid_lsm = ioncdf.opennc(self.mask)
		self.lsm = ioncdf.readnc(fid_lsm,'lsm')
		ioncdf.closenc(fid_lsm)

		print self.name_t2
		print self.name_time_t2
		return None

	def __call__(self):
		if self.dict_input.has_key('file_precip'):
			print 'Rewrite precip file...'
			self.process_precip_file()
		if self.dict_input.has_key('file_snow'):
			print 'Rewrite snow file...'
			self.process_snow_file()
		if self.dict_input.has_key('file_radlw'):
			print 'Create DFS5.2 longwave file...'
			self.process_radlw_file()
		if self.dict_input.has_key('file_radsw'):
			print 'Create DFS5.2 shortwave file...'
			self.process_radsw_file()
		if self.dict_input.has_key('file_t2'):
			print 'Create DFS5.2 t2 file...'
			self.process_t2_file()
		if self.dict_input.has_key('file_q2'):
			print 'Create DFS5.2 q2 file...'
			self.process_q2_file()
		if self.dict_input.has_key('file_msl'):
			print 'Rewrite msl file...'
			self.process_msl_file()
		if self.dict_input.has_key('file_tcc'):
			print 'Rewrite tcc file...'
			self.process_tcc_file()
		if self.dict_input.has_key('file_u10'):
			print 'Create DFS5.2 u10 file...'
			self.process_u10_file()
		if self.dict_input.has_key('file_v10'):
			print 'Create DFS5.2 v10 file...'
			self.process_v10_file()
		return None

	#------------------ Meta functions ------------------------------------------

	def process_tcc_file(self):
		''' Create tcc file '''
		tcc_tmp = np.empty((self.ny,self.nx))
		tcc_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_tcc = ioncdf.opennc(self.file_tcc)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_tcc,'lon')
                lat  = ioncdf.readnc(fid_tcc,'lat')
                time = ioncdf.readnc(fid_tcc,self.name_time_tcc)
		# copy
		for kt in np.arange(0,self.nframes):
			tcc_tmp[:,:]    = ioncdf.readnc_oneframe(fid_tcc,self.name_tcc,kt)
			tcc_out[kt,:,:] = (tcc_tmp[:,:]) * self.lsm[:,:]
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'Cloud','time_dim':'cloud_time','time_var':'cloud_time','long name':'Total Cloud Cover',\
                        'units':'','fileout':self.output_dir + 'tcc_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'tcc','time_dim':'time','time_var':'time','long name':'Total Cloud Cover',\
                        'units':'','fileout':self.output_dir + 'tcc_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = tcc_out.min()
		my_dict['var_valid_max']  = tcc_out.max()
                # close input file and write output
                ioncdf.closenc(fid_tcc)
	        ioncdf.write_ncfile(lon,lat,time,tcc_out,my_dict)
		# clear arrays
		tcc_tmp = None ; tcc_out = None
		return None

	def process_msl_file(self):
		''' Create msl file '''
		msl_tmp = np.empty((self.ny,self.nx))
		msl_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_msl = ioncdf.opennc(self.file_msl)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_msl,'lon')
                lat  = ioncdf.readnc(fid_msl,'lat')
                time = ioncdf.readnc(fid_msl,self.name_time_msl)
		# copy
		for kt in np.arange(0,self.nframes):
			msl_tmp[:,:]    = ioncdf.readnc_oneframe(fid_msl,self.name_msl,kt)
			msl_out[kt,:,:] = (msl_tmp[:,:]) * self.lsm[:,:]
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'Pair','time_dim':'pair_time','time_var':'pair_time','long name':'Mean Sea Level Pressure',\
                        'units':'Pa','fileout':self.output_dir + 'msl_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'msl','time_dim':'time','time_var':'time','long name':'Mean Sea Level Pressure',\
                        'units':'Pa','fileout':self.output_dir + 'msl_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = msl_out.min()
		my_dict['var_valid_max']  = msl_out.max()
                # close input file and write output
                ioncdf.closenc(fid_msl)
	        ioncdf.write_ncfile(lon,lat,time,msl_out,my_dict)
		# clear arrays
		msl_tmp = None ; msl_out = None
		return None

	def process_snow_file(self):
		''' Create snow fall file '''
		nframes   = self.nframes / self.nframes_per_day # daily file
		snow_tmp = np.empty((self.ny,self.nx))
		snow_out = np.empty((nframes,self.ny,self.nx))
                # open input file
                fid_snow = ioncdf.opennc(self.file_snow)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_snow,'lon')
                lat  = ioncdf.readnc(fid_snow,'lat')
                time = ioncdf.readnc(fid_snow,self.name_time_snow)
		# copy
		for kt in np.arange(0,nframes):
			snow_tmp[:,:]    = ioncdf.readnc_oneframe(fid_snow,self.name_snow,kt)
			snow_out[kt,:,:] = (snow_tmp[:,:]) * self.lsm[:,:]
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'rain','time_dim':'rain_time','time_var':'rain_time','long name':'Snow Fall',\
                        'units':'kg.m-2.s-1','fileout':self.output_dir + 'snow_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'snow','time_dim':'time','time_var':'time','long name':'Snow Fall',\
                        'units':'kg.m-2.s-1','fileout':self.output_dir + 'snow_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = snow_out.min()
		my_dict['var_valid_max']  = snow_out.max()
                # close input file and write output
                ioncdf.closenc(fid_snow)
	        ioncdf.write_ncfile(lon,lat,time,snow_out,my_dict)
		# clear arrays
		snow_tmp = None ; snow_out = None
		return None

	def process_precip_file(self):
		''' Create precip file '''
		nframes   = self.nframes / self.nframes_per_day # daily file
		precip_tmp = np.empty((self.ny,self.nx))
		precip_out = np.empty((nframes,self.ny,self.nx))
                # open input file
                fid_precip = ioncdf.opennc(self.file_precip)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_precip,'lon')
                lat  = ioncdf.readnc(fid_precip,'lat')
                time = ioncdf.readnc(fid_precip,self.name_time_precip)
		# copy
		for kt in np.arange(0,nframes):
			precip_tmp[:,:]    = ioncdf.readnc_oneframe(fid_precip,self.name_precip,kt)
			precip_out[kt,:,:] = (precip_tmp[:,:]) * self.lsm[:,:]
                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'rain','time_dim':'rain_time','time_var':'rain_time','long name':'Total Precipitation',\
                        'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'precip','time_dim':'time','time_var':'time','long name':'Total Precipitation',\
                        'units':'kg.m-2.s-1','fileout':self.output_dir + 'precip_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = precip_out.min()
		my_dict['var_valid_max']  = precip_out.max()
                # close input file and write output
                ioncdf.closenc(fid_precip)
	        ioncdf.write_ncfile(lon,lat,time,precip_out,my_dict)
		# clear arrays
		precip_tmp = None ; precip_out = None
		return None

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
			radlw_out[kt,:,:] = (radlw_tmp[:,:] * radlw_factor[:,:] * correct_heat_budget) * self.lsm[:,:]
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
			radsw_out[kt,:,:] = (radsw_tmp[:,:] * radsw_factor[:,:]) * self.lsm[:,:]
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
			u10_out[kt,:,:] = (u10_tmp[:,:] + u10_background[:,:]) * self.lsm[:,:]
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
			v10_out[kt,:,:] = (v10_tmp[:,:] + v10_background[:,:]) * self.lsm[:,:]
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
		elif self.target_model == 'NEMO':
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

		#------- Arctic correction --------
		correction_north = np.zeros((self.ny,self.nx))

               	fid_off = ioncdf.opennc(self.t2_poles_offset)
		monthly_offset_poles = ioncdf.readnc(fid_off,'Tair')
		ioncdf.closenc(fid_off)
		fid_ifra = ioncdf.opennc(self.ice_nsidc)
		monthly_ice_fraction = ioncdf.readnc(fid_ifra,'ifrac')
		ioncdf.closenc(fid_ifra)

		# poles correction north of 70 N
		lattrans3   = 65 # correction poles transition starts at 65, full at 70N
		jtrans3 = (np.abs(lat_s2n-lattrans3)).argmin()
		val_oce = -0.7

		# run the computation
		for kt in np.arange(0,self.nframes):
			t2_tmp[:,:]    = ioncdf.readnc_oneframe(fid_t2,self.name_t2,kt)
			#-----  arctic correction ------
			# interp monthly offset to daily value
			this_day = 1 + (kt / self.nframes_per_day)
			this_day_weights = self.time_interp_weights(this_day)
			this_day_offset  = self.interp_data_to_day(monthly_offset_poles,this_day_weights)
			this_day_ifrac   = self.interp_data_to_day(monthly_ice_fraction,this_day_weights)
			# compute artic correction
			correction_north[jtrans3:,:] = this_day_offset[jtrans3:,:]*this_day_ifrac[jtrans3:,:] + \
			                               val_oce*(1. -this_day_ifrac[jtrans3:,:])
			if not self.target_model == 'ROMS':
				correction_north[:,:] = correction_north[::-1,:]
			#------ add northern and southern hemisphere correction ------
			t2_out[kt,:,:] = (t2_tmp[:,:] + correction_north[:,:] + correction_south[:,:]) * self.lsm[:,:]
			if self.drown:
				t2_out[kt,:,:] = self.drown_wrapper(t2_out[kt,:,:])
			
                # output file informations
		my_dict = {}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['varname']        = self.name_t2
                my_dict['time_dim']       = self.name_time_t2
                my_dict['time_var']       = self.name_time_t2
                my_dict['long name']      = 'Air Temperature at 2m'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['var_valid_min']  = t2_out.min()
		my_dict['var_valid_max']  = t2_out.max()
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['fileout']        = self.processed_nc_dir + \
                                            self.drownstring + self.name_t2 + '_' + self.dataset + '_' + str(self.year) + '.nc'
                if self.target_model == 'ROMS':
                        model_dependent = {'units':'degC','description':my_dict['description'] + '\nROMS-ready ERAinterim forcing'}
                elif self.target_model == 'NEMO':
                        model_dependent = {'units':'K','description':my_dict['description'] + '\nNEMO-ready ERAinerim forcing'}
                my_dict.update(model_dependent)
                # close input file and write output
                ioncdf.closenc(fid_t2)
	        ioncdf.write_ncfile(lon,lat,time,t2_out,my_dict)
		# clear arrays
		t2_tmp = None ; t2_out = None
		# save this for q2
		self.file_t2_new = my_dict['fileout']
		return None

	def process_q2_file(self):
		''' Create q2 file '''
		t2_old = np.empty((self.ny,self.nx))
		t2_new = np.empty((self.ny,self.nx))
		msl    = np.empty((self.ny,self.nx))
		q2_tmp = np.empty((self.ny,self.nx))
		q2_out = np.empty((self.nframes,self.ny,self.nx))
                # open input file
                fid_q2 = ioncdf.opennc(self.file_q2)
                # read coordinates and time
                lon  = ioncdf.readnc(fid_q2,'lon')
                lat  = ioncdf.readnc(fid_q2,'lat')
                time = ioncdf.readnc(fid_q2,self.name_time_q2)

		# We need original t2 and msl from ERAinterim
		fid_t2old = ioncdf.opennc(self.file_t2)
		fid_msl   = ioncdf.opennc(self.file_msl)
		fid_t2new = ioncdf.opennc(self.file_t2_new)
		# run the computation
		for kt in np.arange(0,self.nframes):
			# read all the fields
			# ROMS has t2 in degC, others in Kelvin - qsat_from_t2_and_msl expects Celcius
			if self.target_model == 'ROMS':
				t2_old[:,:]    = ioncdf.readnc_oneframe(fid_t2old,self.name_t2,kt)
				t2_new[:,:]    = ioncdf.readnc_oneframe(fid_t2new,self.name_t2,kt)
			else:
				t2_old[:,:]    = ioncdf.readnc_oneframe(fid_t2old,self.name_t2,kt) - 273.15
				t2_new[:,:]    = ioncdf.readnc_oneframe(fid_t2new,self.name_t2,kt) - 273.15
			msl[:,:]       = ioncdf.readnc_oneframe(fid_msl,self.name_msl,kt)
			q2_tmp[:,:]    = ioncdf.readnc_oneframe(fid_q2,self.name_q2,kt)
			# compute humidity at saturation
			q_sat_new = humidity_toolbox.qsat_from_t2_and_msl(t2_new,msl)
			q_sat_old = humidity_toolbox.qsat_from_t2_and_msl(t2_old,msl)
			# compute new specific humidity
			q2_out[kt,:,:] = (q2_tmp[:,:] * q_sat_new / q_sat_old) * self.lsm[:,:]

                # output file informations
		if self.target_model == 'ROMS':
                        my_dict = {'varname':'Qair','time_dim':'qair_time','time_var':'qair_time','long name':'Specific Humidity',\
                        'units':'kg/kg','fileout':self.output_dir + 'q2_' + self.dataset + '_' + str(self.year) + '_ROMS.nc'}
                elif self.target_model == 'NEMO':
                        my_dict = {'varname':'q2','time_dim':'time','time_var':'time','long name':'Specific Humidity',\
                        'units':'kg/kg','fileout':self.output_dir + 'q2_' + self.dataset + '_' + str(self.year) + '.nc'}
		my_dict['description'] = 'DFS 5.2 (MEOM/LGGE) contact : raphael.dussin@gmail.com'
		my_dict['spval']          = self.spval
		my_dict['reftime']        = self.reftime
		my_dict['time_valid_min'] = time.min()
		my_dict['time_valid_max'] = time.max()
		my_dict['var_valid_min']  = q2_out.min()
		my_dict['var_valid_max']  = q2_out.max()
                # close input file and write output
                ioncdf.closenc(fid_q2)
	        ioncdf.write_ncfile(lon,lat,time,q2_out,my_dict)
		# clear arrays
		q2_tmp = None ; q2_out = None
		return None




	#------------------- Time interpolation ---------------------------------------------------------
        def interp_data_to_day(self,data_monthly,my_weights):
                # time interpolation of bias file
                nt,ny,nx = data_monthly.shape
                if ( nt != 12):
                        print 'Error : number of frames'
                data_thisday = np.zeros((ny,nx))
                for kt in np.arange(12):
                        data_thisday = data_thisday + data_monthly[kt,:,:] * my_weights[kt]
                return data_thisday

        def time_interp_weights(self,nday):
                ''' compute weights to apply to a monthly data file
                nday is numbers of days since begining of the year '''
                # init to zero
                weights = np.zeros((12))
                # cumul of days at the middle of the month
                cumul_days = self.days_in_month.cumsum() - (self.days_in_month / 2.)
                # find which month is the lowest bound of the interval
                # uses the property of python that -1 is last index
                nmonth = -1
                for km in np.arange(12):
                        if ( nday >= cumul_days[km] ):
                                nmonth = nmonth + 1
                # distance from current day to lower bound of interval
                dist_to_lower_bound = np.mod(nday - cumul_days[nmonth],self.ndays)
                # distance between lower and upper bound (uses periodicity properties for nmonths and ndays)
                interval = np.mod(cumul_days[np.mod(nmonth+1,12)] - cumul_days[np.mod(nmonth,12)],self.ndays)
                # weight for lower and upper bound
                wgt_1 = 1.0 - (float(dist_to_lower_bound) / interval )
                wgt_2 = 1.0 - wgt_1
                # change only values for lower and upper bound (all other remain zero)
                weights[np.mod(nmonth,12)] = wgt_1
                weights[np.mod(nmonth+1,12)] = wgt_2
                return weights

        #------------------ wrapper drown -----------------------------------------------

        def drown_wrapper(self,data_in):
                ''' A wrapper for drown to make it cleaner in the main functions '''
                tmp = drwn.drown(0,data_in[:,:].transpose(),self.lsm.transpose(),200,40)
                data_out = tmp.transpose()
                return data_out
