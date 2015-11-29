#!/usr/bin/env python

from lib_process_DFS52 import DFS52_processing
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

# read dfs.datafiles
datafiles = ConfigParser.ConfigParser()
datafiles.read('./dfs.datafiles')

dict_datafiles = {}

for item in datafiles.options(this_user):
        dict_datafiles[item] = datafiles.get(this_user,item)

#--------------------------------------------------------------------------------------------------

if dict_ctl['target_model'] == 'ROMS':
	suffix = '_ROMS.nc'
else:
	suffix = '.nc'

for year in np.arange(fyear,lyear+1):

	my_inputs = {'file_t2':     dict_ctl['processed_nc_dir'] + 'drowned_t2_ERAinterim_'     + str(year) + suffix, \
	             'file_q2':     dict_ctl['processed_nc_dir'] + 'drowned_q2_ERAinterim_'     + str(year) + suffix, \
	             'file_u10':    dict_ctl['processed_nc_dir'] + 'drowned_u10_ERAinterim_'    + str(year) + suffix, \
	             'file_v10':    dict_ctl['processed_nc_dir'] + 'drowned_v10_ERAinterim_'    + str(year) + suffix, \
	             'file_radlw':  dict_ctl['processed_nc_dir'] + 'drowned_radlw_ERAinterim_'  + str(year) + suffix, \
	             'file_radsw':  dict_ctl['processed_nc_dir'] + 'drowned_radsw_ERAinterim_'  + str(year) + suffix, \
	             'file_precip': dict_ctl['processed_nc_dir'] + 'drowned_precip_ERAinterim_' + str(year) + suffix, \
	             'file_snow':   dict_ctl['processed_nc_dir'] + 'drowned_snow_ERAinterim_'   + str(year) + suffix, \
	             'file_msl':    dict_ctl['processed_nc_dir'] + 'drowned_msl_ERAinterim_'    + str(year) + suffix, \
	             'file_tcc':    dict_ctl['processed_nc_dir'] + 'drowned_tcc_ERAinterim_'    + str(year) + suffix, \
	             'output_dir':  dict_ctl['dfs52_output_dir']                                          , \
	             'target_model':dict_ctl['target_model']                                              , \
	             'year':year, 'ncumul':4, 'nx':512, 'ny':256, 'freq':'3h'}

	go = DFS52_processing(my_inputs,dict_datafiles)
	go()
