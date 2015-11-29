import os 

class drown():

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
