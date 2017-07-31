#!/usr/bin/env python
#
# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.
#
# To run this example, you need an API key 
# available from https://api.ecmwf.int/v1/key/
#
# R.Dussin 2014 - adapt script to download ERAinterim

from ecmwfapi import ECMWFDataServer
import numpy as np
import ConfigParser
import os 
import sys

try:
        this_user = sys.argv[1]
except:
        print 'Please provide a user.opts entry to this script' ; exit()

# open the parser and read the info on the run
config = ConfigParser.ConfigParser()
config.read('./user.opts')

dict_ctl = {}

for item in config.options(this_user):
        dict_ctl[item] = config.get(this_user,item)

# Define the list of interesting variables with their associated parameter on MARS server
ecmwf_param = {'u10' : '165.128', 'v10' : '166.128', 'd2' : '168.128' , 't2' : '167.128' , \
                'msl' : '151.128', 'snow' : '144.128' , 'radsw' : '169.128' , 'radlw' : '175.128' , 'precip' : '228.128' , 'tcc' : '164.128'}

# Choose years to download
fyear = int(dict_ctl['first_year'])
lyear = int(dict_ctl['last_year'])

server = ECMWFDataServer()

for year in np.arange(fyear,lyear+1):

	for key in ecmwf_param.keys():

		print 'working on variable', key, ' for year ', str(year)
		filegrib = dict_ctl['original_grib_dir'] + key + '_ERAinterim_' + str(year) + '.grib'
		filenc   = dict_ctl['original_nc_dir']   + key + '_ERAinterim_' + str(year) + '.nc'

		server.retrieve({
		    'stream'    : "oper",
		    'levtype'   : "sfc",
	            'resol'      : "av",
		    'param'     : ecmwf_param[key],
		    'dataset'   : "interim",
		    'step'      : "3/6/9/12",
		    'time'      : "00/12",
		    'date'      : str(year) + "-01-01/to/" + str(year) + "-12-31",
		    'type'      : "fc",
		    'class'     : "ei",
		    'target'    : filegrib
		})

		# comand to run to convert to NC : 
		convert_to_nc='cdo -R -t ecmwf -f nc -r copy ' + filegrib + ' ' + filenc
		os.system(convert_to_nc)


