#!/bin/bash

# Download from web
wget http://research.jisao.washington.edu/data_sets/iabppoles/satiabppolesclim.nc

# create land sea mask (needed by sosie)
python create_mask_POLES.py 

# interpolate to ERAinterim grid
sosie.x -f namelist.POLES
