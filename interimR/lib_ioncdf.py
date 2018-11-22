import netCDF4 as nc


def opennc(myfile):
    fid = nc.Dataset(myfile,'r')
    return fid

def closenc(fid):
    fid.close()
    return None

def readnc_frames(fid,myvar,myframe_start,myframe_end):
    ''' read data from netcdf '''
    out = fid.variables[myvar][myframe_start:myframe_end,:,:].squeeze()
    return out

def readnc_oneframe(fid,myvar,myframe):
    ''' read data from netcdf '''
    out = fid.variables[myvar][myframe,:,:].squeeze()
    return out

def readnc(fid,myvar):
    ''' read data from netcdf '''
    out = fid.variables[myvar][:].squeeze()
    return out

def create_ncfile(lon_array,lat_array,time,dict_wrt):
    fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
    fid.description = dict_wrt['description']
    # dimensions
    fid.createDimension('lat', lat_array.shape[0])
    fid.createDimension('lon', lon_array.shape[0])
    fid.createDimension(dict_wrt['time_dim'], None)
    # variables
    latitudes  = fid.createVariable('lat', 'f8', ('lat',))
    longitudes = fid.createVariable('lon', 'f8', ('lon',))
    times      = fid.createVariable(dict_wrt['time_var'], 'f8', (dict_wrt['time_dim'],))
    variable   = fid.createVariable(dict_wrt['varname'], 'f4', (dict_wrt['time_dim'],'lat','lon',), \
                                fill_value=dict_wrt['spval'])
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

    times.units = "days since " + dict_wrt['reftime'].isoformat()
    times.calendar = "LEAP"
    times.valid_min = dict_wrt['time_valid_min']
    times.valid_max = dict_wrt['time_valid_max']

    variable.long_name = dict_wrt['long name']
    variable.units = dict_wrt['units']
    variable.coordinates = "lon lat"
    variable.time = dict_wrt['time_var']
    variable.missing_value = dict_wrt['spval']
    variable.valid_range = dict_wrt['var_valid_min'], dict_wrt['var_valid_max']
    return fid, times, variable

def finalize_ncfile(fid):
    fid.close()
    return None

def write_ncfile(lon_array,lat_array,time,var,dict_wrt):
    fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
    fid.description = dict_wrt['description']
    # dimensions
    fid.createDimension('lat', lat_array.shape[0])
    fid.createDimension('lon', lon_array.shape[0])
    fid.createDimension(dict_wrt['time_dim'], None)
    # variables
    latitudes  = fid.createVariable('lat', 'f8', ('lat',))
    longitudes = fid.createVariable('lon', 'f8', ('lon',))
    times      = fid.createVariable(dict_wrt['time_var'], 'f8', (dict_wrt['time_dim'],))
    variable   = fid.createVariable(dict_wrt['varname'], 'f4', (dict_wrt['time_dim'],'lat','lon',),\
                                    fill_value=dict_wrt['spval'])

    # attributes
    longitudes.units = "degrees_east"
    longitudes.valid_min = lon_array.min()
    longitudes.valid_max = lon_array.max()
    longitudes.long_name = "longitude"

    latitudes.units = "degrees_north"
    latitudes.valid_min = lat_array.min()
    latitudes.valid_max = lat_array.max()
    latitudes.long_name = "latitude"

    times.units = "days since " + dict_wrt['reftime'].isoformat()
    times.valid_min = time.min()
    times.valid_max = time.max()
    times.calendar = "LEAP"

    variable.long_name = dict_wrt['long name']
    variable.units = dict_wrt['units']
    variable.coordinates = "lon lat"
    variable.time = dict_wrt['time_var']
    variable.missing_value = dict_wrt['spval']
    variable.valid_range = var.min() , var.max()

    # data
    latitudes[:]    = lat_array
    longitudes[:]   = lon_array
    times[:]        = time
    variable[:,:,:] = var

    # close
    fid.close()
    return None

def write_lsm(lon_array,lat_array,var,dict_wrt):
    fid = nc.Dataset(dict_wrt['fileout'], 'w', format='NETCDF3_CLASSIC')
    fid.description = dict_wrt['description']
    # dimensions
    fid.createDimension(dict_wrt['dim_lat'], lat_array.shape[0])
    fid.createDimension(dict_wrt['dim_lon'], lon_array.shape[0])
    # variables
    latitudes  = fid.createVariable(dict_wrt['var_lat'], 'f8', (dict_wrt['dim_lat'],))
    longitudes = fid.createVariable(dict_wrt['var_lon'], 'f8', (dict_wrt['dim_lon'],))
    variable   = fid.createVariable(dict_wrt['var_lsm'], 'i4', (dict_wrt['dim_lat'],dict_wrt['dim_lon'],))

    # attributes
    longitudes.units = "degrees_east"
    longitudes.valid_min = lon_array.min()
    longitudes.valid_max = lon_array.max()
    longitudes.long_name = "longitude"

    latitudes.units = "degrees_north"
    latitudes.valid_min = lat_array.min()
    latitudes.valid_max = lat_array.max()
    latitudes.long_name = "latitude"

    variable.long_name = dict_wrt['long name']
    variable.units = dict_wrt['units']
    variable.coordinates = dict_wrt['dim_lat'] + ' ' + dict_wrt['dim_lon']
    variable.valid_range = var.min() , var.max()

    # data
    latitudes[:]    = lat_array
    longitudes[:]   = lon_array
    variable[:,:]   = var

    # close
    fid.close()
    return None
