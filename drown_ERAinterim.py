#!/usr/bin/env python

from lib_process_ERAinterim import ERAinterim_drown
import numpy as np
import ConfigParser
import sys 

try:
	this_user = sys.argv[1]
except:
	print 'You must provide a user options entry for this script' ; exit()

# open the parser and read the info on the run
config = ConfigParser.ConfigParser()
config.read('./user.opts')

dict_ctl = {}

for item in config.options(this_user):
	dict_ctl[item] = config.get(this_user,item)

fyear = int(dict_ctl['first_year'])
lyear = int(dict_ctl['last_year'])

#--------------------------------------------------------------------------------------------------

for year in np.arange(fyear,lyear+1):

	if dict_ctl['target_model'] == 'ROMS':
		my_inputs = {'file_t2':     dict_ctl['processed_nc_dir'] + 't2_ERAinterim_'     + str(year) + '_ROMS.nc', \
		             'file_q2':     dict_ctl['processed_nc_dir'] + 'q2_ERAinterim_'     + str(year) + '_ROMS.nc', \
		             'file_u10':    dict_ctl['processed_nc_dir'] + 'u10_ERAinterim_'    + str(year) + '_ROMS.nc', \
		             'file_v10':    dict_ctl['processed_nc_dir'] + 'v10_ERAinterim_'    + str(year) + '_ROMS.nc', \
		             'file_radlw':  dict_ctl['processed_nc_dir'] + 'radlw_ERAinterim_'  + str(year) + '_daily_ROMS.nc', \
		             'file_radsw':  dict_ctl['processed_nc_dir'] + 'radsw_ERAinterim_'  + str(year) + '_daily_ROMS.nc', \
		             'file_precip': dict_ctl['processed_nc_dir'] + 'precip_ERAinterim_' + str(year) + '_daily_ROMS.nc', \
		             'file_snow':   dict_ctl['processed_nc_dir'] + 'snow_ERAinterim_'   + str(year) + '_daily_ROMS.nc', \
		             'file_msl':    dict_ctl['processed_nc_dir'] + 'msl_ERAinterim_'    + str(year) + '_ROMS.nc', \
		             'output_dir':  dict_ctl['processed_nc_dir']                                                , \
		             'sosie_dir':   dict_ctl['sosie_dir']                                                       , \
		             'lsm_file':    dict_ctl['lsm_file']                                                        , \
		             'target_model':dict_ctl['target_model']                                                    , \
		             'year':year, 'ncumul':4, 'nx':512, 'ny':256, 'freq':'3h'}
	else:
		my_inputs = {'file_t2':     dict_ctl['processed_nc_dir'] + 't2_ERAinterim_'     + str(year) + '.nc', \
		             'file_q2':     dict_ctl['processed_nc_dir'] + 'q2_ERAinterim_'     + str(year) + '.nc', \
		             'file_u10':    dict_ctl['processed_nc_dir'] + 'u10_ERAinterim_'    + str(year) + '.nc', \
		             'file_v10':    dict_ctl['processed_nc_dir'] + 'v10_ERAinterim_'    + str(year) + '.nc', \
		             'file_radlw':  dict_ctl['processed_nc_dir'] + 'radlw_ERAinterim_'  + str(year) + '_daily.nc', \
		             'file_radsw':  dict_ctl['processed_nc_dir'] + 'radsw_ERAinterim_'  + str(year) + '_daily.nc', \
		             'file_precip': dict_ctl['processed_nc_dir'] + 'precip_ERAinterim_' + str(year) + '_daily.nc', \
		             'file_snow':   dict_ctl['processed_nc_dir'] + 'snow_ERAinterim_'   + str(year) + '_daily.nc', \
		             'file_msl':    dict_ctl['processed_nc_dir'] + 'msl_ERAinterim_'    + str(year) + '.nc', \
		             'output_dir':  dict_ctl['processed_nc_dir']                                           , \
		             'sosie_dir':   dict_ctl['sosie_dir']                                                  , \
		             'lsm_file':    dict_ctl['lsm_file']                                                   , \
		             'target_model':dict_ctl['target_model']                                               , \
		             'year':year, 'ncumul':4, 'nx':512, 'ny':256, 'freq':'3h'}


	go = ERAinterim_drown(my_inputs)
	go()
