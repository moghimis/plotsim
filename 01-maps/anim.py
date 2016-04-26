#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Plot Maps code

import os,sys
import netCDF4
import netcdftime

file_2_run='maps_assim_plot.py'

ncf='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/03_wave_limit_test/01_wave_hs0.5_t14/run_1000/08_forward/nri_his.nc'
nc=netCDF4.Dataset(ncf)
ncvar=nc.variables
timename_sim='ocean_time'
utime_sim=netcdftime.utime(ncvar[timename_sim].units)
times=ncvar[timename_sim][:]
sdates_sim=utime_sim.num2date(times) 

os.system('rm log.txt')
for i in range(len(sdates_sim[15:])):
    comm='python '+file_2_run+'  '+sdates_sim[i+15].isoformat()+' >> log.txt'
    print '>>>>>>>>>>>>> > > > > ',comm
    os.system(comm)
    