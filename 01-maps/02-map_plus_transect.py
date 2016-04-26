import netCDF4
import netcdftime
import numpy as np
from   datetime import datetime
import sys,os
from math import ceil
#import plotting stuff
from octant import plotting as ocpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import string
import okean


#time_now=datetime.now().isoformat()[5:16]
global rotate2nri
#import base_info
global base_dir
global inp_dir

rotate2nri = True
assim_ebb  = False
normalize  = False
vec        = True
plot_mask  = True
plot_contour      = True


#ftypes=['pdf','png']
ftypes=['png']

#For syn tests
if False:
    real_data = False
    #for syn assim paper
    prior1  = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5_10iteration_cross_sections_flood_wave/run_1000/08_forward/nri_his.nc'
    cw5flo = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5_10iteration_cross_sections_flood_wave/run_1004_flood_wave_low_err/08_forward/nri_his.nc'
    cw5ebb = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5iteration_cross_sections/run_1004_ebb_wave/08_forward/nri_his.nc'
    wav5   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5iteration_cross_sections/run_1004_wave/08_forward/nri_his.nc'
    flo5   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5iteration_cross_sections/run_1004_flood/08_forward/nri_his.nc'
    ebb5   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/10_only_ebb_iterate_alex_added_noise_err_reduction5/run_1005/08_forward/nri_his.nc'
    true   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/01_syn_4paper_5_10iteration_cross_sections_flood_wave/run_1020/08_forward/nri_his.nc'
    files=[prior1  ,   cw5flo,    cw5ebb,  wav5,    flo5, ebb5 , true]
    varn =['prior', 'CW5-flood', 'CW5-ebb','Wav5','Flood5','Ebb5','True']
    #itr >  0            1           2        3      4        5     6        
    itr=5

#For tru_tru tests
if False:
#     real_data = False
#     #for syn assim paper
#     prior1  = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1005/08_forward/nri_his.nc'
#     cw5flo = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1000_flood_wav/08_forward/nri_his.nc'
#     cw5ebb = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1000_ebb_wave/08_forward/nri_his.nc'
#     flo5   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1000_flood/08_forward/nri_his.nc'
#     ebb5   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1000_ebb/08_forward/nri_his.nc'
#     true   = '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/run_1005/08_forward/nri_his.nc'
#     files=[ prior1,      cw5flo,    cw5ebb,  flo5   , ebb5 , true]
#     varn =['prior', 'CW1-flood', 'CW1-ebb',  'Flood1','Ebb1','True']
#     #itr >  0            1           2        3      4        5     6        
#     itr=5
    real_data = False
    #for syn assim paper
    prior1  = '/home/server/pi/homes/moghimi/work/00-projs/01-muri/00-progs/cowast/project2/07-assim/scr_generations/set4_all/plotsim/cvh/assim_tru_tru_new_covh/run_1000/08_forward/nri_his.nc'
    cw5flo = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/99_assim_tru_tru_old/run_1000_flood/08_forward/nri_his.nc'
    cw5ebb = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/99_assim_tru_tru_old/run_1000_ebb/08_forward/nri_his.nc'
    flo5   = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/802_tru_tru_simple_alpha01_hcrit0.25/run_1002_flood/08_forward/nri_his.nc'
    ebb5   = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/802_tru_tru_simple_alpha01_hcrit0.25/run_1001_ebb/08_forward/nri_his.nc'
    flo6   = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/700_tru_tru_simple_alpha_2/run_1002/08_forward/nri_his.nc'
    ebb6   = '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/700_tru_tru_simple_alpha_2/run_1001/08_forward/nri_his.nc'
    ###################
    true   = 'none'
    files=[ prior1,   cw5flo,    cw5ebb,         flo5      , ebb5       ,  flo6          , ebb6      , prior1]
    varn =['True' , 'old-flood', 'old-ebb',  'cutoff-flood','cutoff-ebb',  'alpha2-flood','alpha2-ebb','True']
    ntr = [1,2,3,4,5,6]
    ntr = [5,6]
    #vars = ['h','uv']
    vars = ['h','uv','h_dif','uv_dif']
    #files=[ prior1,  flo5   , ebb5 , prior1]
    #varn =['True','new-flood','new-ebb','True']
    #itr >  0            1           2        3      4        5     6        
    itr=3

if False:
    #for DC meeting with real data
    real_data = True
    base_dir='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run7_real_paper/04_long_run_1itr_comparison_full_data/'
    #base_dir='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run7_real_paper/04_long_run_1itr_comparison_inside_cur_data_wave_all/'
    prior1       = base_dir + 'run_1000/08_forward/nri_his.nc'
    wave         = base_dir + 'run_1001/08_forward/nri_his.nc' 
    swift        = base_dir + 'run_1002/08_forward/nri_his.nc' 
    ebb          = base_dir + 'run_1003/08_forward/nri_his.nc' 
    swf_wav      = base_dir + 'run_1004/08_forward/nri_his.nc' 
    ebb_wav      = base_dir + 'run_1005/08_forward/nri_his.nc' 
    ebb_wav_swf  = base_dir + 'run_1006/08_forward/nri_his.nc' 
    true         = base_dir + 'run_1007/08_forward/nri_his.nc' 
    files=[prior1  ,  wave,  swift ,  ebb,    swf_wav, ebb_wav , ebb_wav_swf,  true]
    varn =['prior', 'Wav1', 'SWF1','Ebb1',  'SWF_Wav','Ebb_Wav','Ebb_Wav_SWF','True']
    #itr >  0          1      2       3        4        5          6            7     
    #Choose itr to be ploted
    ntr = [1,2,3,4,5,6]

if False:
    #for new cov compare ebb flood with hmin const = 1.5 m
    real_data = False
    base_dir='/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/690_compare_new_syn_ebb_flood/'
    prior1       = base_dir + 'run_1000/08_forward/nri_his.nc'
    wave         = base_dir + 'run_1001/08_forward/nri_his.nc' 
    swift        = base_dir + 'run_1002/08_forward/nri_his.nc' 
    ebb          = base_dir + 'none' 
    swf_wav      = base_dir + 'none' 
    ebb_wav      = base_dir + 'none' 
    ebb_wav_swf  = base_dir + 'none' 
    true         = base_dir + 'run_1020/08_forward/nri_his.nc' 
    files=[prior1  ,  wave,  swift ,  ebb,  true]
    varn =['prior', 'CW0-Ebb', 'CW0-Flood' ,'True']
    #itr >  0           1          2           3        4        5          6            7     
    #Choose itr to be ploted
    ntr = [1,2]
    vars = ['h','uv','h_dif','uv_dif']

if False:
    #for submitted paper
    #for new cov only ebb with hmin const = 1.5 m
    real_data = False
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases_hmin_create_bathy_were_zero/06_wave_current_1st_itr_final_paper/690_compare_new_covhh_syn_ebb_flood/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his.nc' 
    true    = base_dir + 'run_1020/08_forward/nri_his.nc' 
    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6    , ebb7      ,true]
    varn =['prior','Ebb','Flood','WAV1-Tm01=4s','WAV-Tm01=7s','WAV-Tm01=14s','CW-ebb','CW-flood', 'True']
    #itr >  0           1          2           3        4        5          6            7     
    #Choose itr to be ploted
    ntr = [1,2,3,4,5,6,7]
    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['h','uv']

if False:
    #for new cov only ebb with hmin const = 1.5 m
    real_data = False
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data/02_compare_1st_itr/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his.nc' 
    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6       , ebb7 ]
    varn =['prior','SAR', 'SWF', 'WAV',  'SAR-WAV','SWF-WAV','SAR-SWF-WAV', 'True']
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [1,2,3,4,5,6]
    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['uv']

if False:
    #for new cov only ebb with hmin const = 1.5 m
    real_data = True
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_10days_v2/02_compare_1st_itr/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his.nc' 
    ebb8    = base_dir + 'run_1008/08_forward/nri_his.nc' 
    ebb9    = base_dir + 'run_1009/08_forward/nri_his.nc' 
    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6         , ebb7   , ebb8   , ebb9  ]
    varn =['prior','SAR', 'SWF', 'WAV',  'SAR-WAV','SWF-WAV','SAR-SWF-WAV' ,'SAR-SWF-WAV_0.5Err' ,'SAR-SWF-WAV_0.25Err', 'True']
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [0,1,2,3,4,5,6,7,8]
    #ntr = [0,1]

    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['h','uv']


if False:
    #for new cov only ebb with hmin const = 1.5 m
    real_data = True
    #2D
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_2D_3days_01/01_main_ebb_2d_prior_true/'
    ntr = [0,1,2,3,4,5,6,7]

    #3D
    #base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_01/01_main_ebb_3d_prior_true/'
    #ntr = [0,1,2,3,4,7]


    prior1  = base_dir + 'run_1000/08_forward/nri_his_0001.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his_0001.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his_0001.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his_0001.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his_0001.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his_0001.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his_0001.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his_0001.nc' 

    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6         , ebb7   ]
    varn =[
            'Prior',
            'SAR-RAD1',
            'SAR-RAD2',
            'SAR-RAD3',
            'SAR-RAD4',
            'SAR-RAD5',
            'SAR-RAD6',
            'True',
            'SAR'
            ]
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [0,1,2,3,4,5,6,7]
    #ntr = [0,1]

    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['h','uv']
    #vars = ['h']


if False:
    real_data = True
    
    #base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_01/02_compare_1st_itr/'
    #base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_01/02_compare_2nd_itr/'
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_01/02_compare_3rd_itr/'
    #base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_01/02_compare_4th_itr/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his_0001.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his_0001.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his_0001.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his_0001.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his_0001.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his_0001.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his_0001.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his_0001.nc' 

    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6  , ebb7   ]
    varn =[
            'Prior',
            'SAR',
            'SWF',
            'RAD',
            'SAR-RAD',
            'SWF-RAD',
            'SAR-SWF-RAD',
            'True',
            ]
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [0,1,2,3,4,5,6]
    #ntr = [0,1]

    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['h','uv']
    #vars = ['h']

if False:
    real_data = False
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/wave_current/run6_syn_using_3D_runs/Compare_3d_4_paper/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his_0001.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his_0001.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his_0001.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his_0001.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his_0001.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his_0001.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his_0001.nc' 
    ebb7    = base_dir + 'run_1007/08_forward/nri_his_0001.nc' 
    ebb8    = base_dir + 'run_1008/08_forward/nri_his_0001.nc' 

    files=[prior1 ,ebb1   ,ebb2  ,ebb3   , ebb4    , ebb5   , ebb6 , ebb7  , ebb8]
    varn=['Prior',
          'CUR-ebb',
          'CUR-flood',
          'WAV-Tm01=4s',
          'WAV-Tm01=7s',
          'WAV-Tm01=14s',
          'CW-ebb',
          'CW-flood',
          'True']
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [0,1,2,3,4,5,6,7]
    #ntr = [0,1]

    #vars = ['h','uv','h_dif','uv_dif']
    #vars = ['h','uv']
    vars = ['h']


if True:
    real_data = True
    
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/real_data_3D_3days_02_new_prior/zz_comp_itrations/'
    prior1  = base_dir + 'run_1000/08_forward/nri_his_0001.nc'
    ebb1    = base_dir + 'run_1001/08_forward/nri_his_0001.nc' 
    ebb2    = base_dir + 'run_1002/08_forward/nri_his_0001.nc' 
    ebb3    = base_dir + 'run_1003/08_forward/nri_his_0001.nc' 
    ebb4    = base_dir + 'run_1004/08_forward/nri_his_0001.nc' 
    ebb5    = base_dir + 'run_1005/08_forward/nri_his_0001.nc' 
    ebb6    = base_dir + 'run_1006/08_forward/nri_his_0001.nc' 
    #ebb7    = base_dir + 'run_1007/08_forward/nri_his_0001.nc' 

    files=[prior1 ,ebb1   ,ebb2  ,ebb3 ]  #, ebb4]#    , ebb5   , ebb6  ]#, ebb7   ]
    varn =[
            'Prior',
            'Itr1',
            'Itr2',
            'Itr3',
            'Itr4',
            'Itr5',
            'True',
            ]
    #itr >  0           1          2           3        4        5          6     
    #Choose itr to be ploted
    ntr = [0,1,2,3] #,4]
    #ntr = [0,1]

    #vars = ['h','uv','h_dif','uv_dif']
    vars = ['h','uv']
    #vars = ['h']




surface_vel       = True

#Choose itr to be ploted

#### Functions ####
def read_nc_file(ncf,date_assim,surface_vel):
    #Read in his_file
    #base=base_dir+'/run_'+str(1000+itr)+'/08_forward/'
    #ncf=base+'nri_his.nc'
    #ncf=file
    #read data from nc files
    #if copy1: os.system('cp '+ncf+'   nri_his_'+str(1000+itr)+'.nc')
    nc=netCDF4.Dataset(ncf)
    ncvar=nc.variables
    timename_sim='ocean_time'
    utime_sim=netcdftime.utime(ncvar[timename_sim].units)
    times=ncvar[timename_sim][:]
    sdates_sim=utime_sim.num2date(times) 
    nt=len(sdates_sim)
    #print sdates_sim
    #To find closest time step to assimilation time
    sec_assim =utime_sim.date2num(date_assim)
    diff=np.abs(sec_assim-times)
    ind=np.where(diff==diff.min())
    n=np.array(ind).item()
    print '     >','chosen time step from netCDF file > ',n,sdates_sim[n]
    
    maskr=ncvar['mask_rho'][:]
    xr=ncvar['x_rho'][:]
    yr=ncvar['y_rho'][:]
    
    
    if surface_vel:
        ubar=ncvar['u']  [n,-1,:]
        vbar=ncvar['v']  [n,-1,:]      
    else:
        ubar=ncvar['ubar'][n,:]
        vbar=ncvar['vbar'][n,:]
    
    
    z=ncvar['zeta'][n,:]
    h=ncvar['h'][:]
    u=np.zeros_like(z)
    v=np.zeros_like(z)
    
    #to avoid mask areas come to the averaging later on.
    ubar[ubar>5]=0
    vbar[vbar>5]=0
    #interpolate to rho points (aprox.)
    u[:,1:-1]=(ubar[:,:-1]+ubar[:,1:])/2.
    u[:,0]   = ubar[:,0]
    u[:,-1]  = ubar[:,-1]
    
    v[1:-1,:]=(vbar[:-1,:]+vbar[1:,:])/2.
    v[0   ,:]=vbar [0  ,:]
    v[-1  ,:]=vbar [-1 ,:]   
    
    m_cut=maskr [j1:j2:k,i1:i2:k]
    x_cut=xr    [j1:j2:k,i1:i2:k]
    y_cut=yr    [j1:j2:k,i1:i2:k]
    h_cut=h     [j1:j2:k,i1:i2:k]
    u_cut=u     [j1:j2:k,i1:i2:k]
    v_cut=v     [j1:j2:k,i1:i2:k]
    z_cut=z     [j1:j2:k,i1:i2:k]
    
    #x_cut=np.ma.masked_where(m_cut==0,x_cut)
    #y_cut=np.ma.masked_where(m_cut==0,y_cut)
    h_cut=np.ma.masked_where(m_cut==0,h_cut)
    u_cut=np.ma.masked_where(m_cut==0,u_cut)
    v_cut=np.ma.masked_where(m_cut==0,v_cut)
    z_cut=np.ma.masked_where(m_cut==0,z_cut)
    
    
    u_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)) ,u_cut )
    v_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)) ,v_cut )
    
    
    nc.close()
    
    
    if rotate2nri:
        x_cut,y_cut=okean.calc.rot2d(x_cut,y_cut,ang=-np.pi/2.0)
        u_cut,v_cut=okean.calc.rot2d(u_cut,v_cut,ang=-np.pi/2.0)

    return x_cut,y_cut,h_cut,u_cut,v_cut,z_cut

def statatistics(val_da,val_mo):
    data_me=val_da
    data_mo=val_mo
    
    # Calculate statistics
    mean_me = data_me.mean()
    mean_mo = data_mo.mean()
    sd1 = data_me.std()
    sd2 = data_mo.std()
    delta = data_mo-data_me
    R2 = 1.-((data_me-data_mo)**2).sum()/((data_me-mean_me)**2).sum()
    
    # Print statistics
    unit='m'
    bias=mean_mo-mean_me
    #print 'Bias = %.3f %s' % (bias,unit)
    #bias_txt='Bias = %.3f %s' % (bias,unit)
    rmse=np.sqrt((delta**2).mean())
    #print 'RMSE = %.3f %s' % (rmse,unit)
    #rmse_txt='RMSE = %.3f %s' % (rmse,unit)
    nrmse1 = rmse / (data_me.max() - data_me.min())
    #print 'NRMSE1 =%s' % (nrmse1   )
    nrmse2 = rmse/ abs(data_me.mean() )
    #print 'NRMSE2 =%s' % (nrmse2 )
    nrmse3 = np.sqrt( (delta**2).sum() / (data_me**2).sum() )
    #print 'NRMSE3 =%s' % (nrmse3 )
    #index of agreement
    #nrmse3 = 1- (abs(delta) **2).sum() /((abs(data_mo-mean_me)+abs(data_me-mean_me) )**2).sum()
    big_error=(np.abs(delta)).max()
    #print 'Biggest error = ', ( big_error )
    mae = np.abs(delta).mean()
    #print 'MAE = %s %s' % (mae,unit)
    cor1 = (((data_me-mean_me)*(data_mo-mean_mo)).mean()/sd1/sd2)
    #print 'Correlation = %s' % (cor1)
    #print 'Coefficient of determination (R2) = %s' % R2
    return bias,rmse

   
def interp3(x,y,b,xnew,ynew,method):
    xf=x.flatten()
    yf=y.flatten()
    bf=b.flatten()
    interpm=method
    if interpm=='tri':
        xnewf=xnew.flatten()
        ynewf=ynew.flatten()
        #print 'tri interp ...'
        from delaunay import  triangulate
        tri=triangulate.Triangulation(xf, yf)
        interp_b=tri.nn_extrapolator(bf)
        bnewf = interp_b(xnewf,ynewf)
        bnew=bnewf.reshape(xnew.shape)
    elif interpm=='csa' :
        #print 'csa interp ...'
        import octant.csa as csa
        csa_interp = csa.CSA(xf, yf,bf)
        bnew = csa_interp(xnew,ynew)
    elif interpm=='grd' :
        #print 'griddata interp ...'
        bnew=pl.griddata(xf,yf,bf,xnew,ynew)
    return bnew

######################################################
################   FUNCS     #########################
######################################################

#sys.path.append('/home/server/pi/homes/moghimi/Desktop/cowast/projects/my_py_libs/plot')
#import plot_settings as ps
from pynmd.plotting import plot_settings as ps 


plt.close('all')
########  read directories  ##################
sys.path.append('../../pysim')

# base directory for this set or iterations
#base_dir=base_info.base_dir
#inp_dir=base_info.inp_dir
#prior=base_info.prior
   
global i1,i2,j1,j2,k
#slice information        
i1 = 35
i2 = 165
j1 = 40
j2 = 175

k  = 1

#cmap info
cmapj=plt.cm.jet_r
cmapo=ocpl.jetWoGn(reverse=True)

for itr in ntr:
    if real_data:
        sar_nc  = netCDF4.Dataset(base_dir+'/inp/obs/sar/uASAR.nc')
        sar_ncv = sar_nc.variables
        x_sar   = sar_ncv['x'][:].squeeze()
        y_sar   = sar_ncv['y'][:].squeeze()
        u_sar   = sar_ncv['u'][:].squeeze()
        v_sar   = sar_ncv['v'][:].squeeze()
        u_sar   = np.ma.masked_where(u_sar > 4 , u_sar)
        v_sar   = np.ma.masked_where(u_sar > 4 , v_sar)
        m_sar   = u_sar.mask
        ut = netcdftime.utime(sar_ncv['time'].units)
        date_sar = ut.num2date(sar_ncv['time'][0])
        print 'Date_sar > ', date_sar.isoformat()
        if rotate2nri:
            x_sar,y_sar = okean.calc.rot2d(x_sar,y_sar,ang=-np.pi/2.0)
            u_sar,v_sar = okean.calc.rot2d(u_sar,v_sar,ang=-np.pi/2.0)
    else:
        print 'assim > '     , varn [itr]
        print 'assim file > ', files[itr]
        if real_data:
            date_sar=datetime(2012,05,10,22,22,31)
        else:
            if 'ebb' in varn[itr]:
                date_sar=datetime(2012,05,10,20)   #ebb
            else:
                date_sar=datetime(2012,05,11,1,30)   #flood        
    
    ###################################################################################
    x_cut,y_cut,h_cut,u_cut,v_cut,z_cut=read_nc_file(files[itr] ,date_sar,surface_vel)
    x_tru,y_tru,h_tru,u_tru,v_tru,z_tru=read_nc_file(files[-1]  ,date_sar,surface_vel)
    
    from_ensemble = False
    if from_ensemble:
        mem_adj = base_dir+'/run_'+str(1000+itr)+'/04_mem_adj/'
        ncf = mem_adj + 'cur_members_prior.nc'
        nc=netCDF4.Dataset(ncf)
        ncvar=nc.variables
        
        x2  = ncvar['x_rho'][:]
        y2  = ncvar['y_rho'][:]
        u_mems_mean  = ncvar['u_mems'][:].mean(axis=2)  
        v_mems_mean  = ncvar['v_mems'][:].mean(axis=2)  
        h_mems_mean  = ncvar['h_mems'][:].mean(axis=2)  
    
        if rotate2nri:
            x2r ,y2r  = okean.calc.rot2d(x2 ,y2 ,ang=-np.pi/2.0)   
            u_mems_meanr, v_mems_meanr =   okean.calc.rot2d(u_mems_mean, v_mems_mean,ang=-np.pi/2.0)
    
        x_pri = x2r              [j1:j2:k,i1:i2:k]
        y_pri = y2r              [j1:j2:k,i1:i2:k]
        h_pri = h_mems_mean      [j1:j2:k,i1:i2:k]
        u_pri = u_mems_meanr     [j1:j2:k,i1:i2:k]
        v_pri = v_mems_meanr     [j1:j2:k,i1:i2:k]
        z_pri = np.ones_like(h_pri) * -100.0
    else:
        x_pri,y_pri,h_pri,u_pri,v_pri,z_pri=read_nc_file(files[0] ,date_sar,surface_vel)
    
    
    
    #Part two
    all_trans = {}
    #sys.path.append('../transect')
    import data1
    for case in ['pri','asi','tru']:
        if   case=='pri':  xr,yr,h,u,v,z=x_pri,y_pri,h_pri,u_pri,v_pri,z_pri
        elif case=='asi':  xr,yr,h,u,v,z=x_cut,y_cut,h_cut,u_cut,v_cut,z_cut
        elif case=='tru':  xr,yr,h,u,v,z=x_tru,y_tru,h_tru,u_tru,v_tru,z_tru
        
        for trans1 in ['main','secn','crBB','crA1','crA2','crA3']:    
            print '>>>> cut for > ',trans1
            if trans1=='main':
                #main_ch_xf
                xtrans=data1.main_ch_x
                ytrans=data1.main_ch_y
                cdfname='trans_main.nc'
            elif trans1=='secn':    
                #second_ch_xf
                xtrans=data1.Second_ch_x
                ytrans=data1.Second_ch_y
                cdfname='trans_second.nc'
            elif trans1=='crBB':    
                #second_ch_xf
                xtrans=data1.crossBB_x
                ytrans=data1.crossBB_y
                cdfname='trans_crossBB.nc'
            elif trans1=='crA1':    
                #second_ch_xf
                xtrans=data1.crossA1_x
                ytrans=data1.crossA1_y
                cdfname='trans_crossA1.nc'
            elif trans1=='crA2':    
                #second_ch_xf
                xtrans=data1.crossA2_x
                ytrans=data1.crossA2_y
                cdfname='trans_crossA2.nc'
            elif trans1=='crA3':    
                #second_ch_xf
                xtrans=data1.crossA3_x
                ytrans=data1.crossA3_y
                cdfname='trans_crossA3.nc'
        
            y1=ytrans
            indy = zip(*sorted([(val, i) for i, val in enumerate(y1)]))[1]
            
            indy=np.array(indy)
            y11=y1[indy]
            x11=xtrans[indy]
            
            yin=np.linspace(y11.min(),y11.max(), 100)
            xin=np.interp( yin, y11,x11)
        
            [xin_shape]=xin.shape
            u_all=np.zeros((xin_shape),dtype=float)
            v_all=np.zeros((xin_shape),dtype=float)
            z_all=np.zeros((xin_shape),dtype=float)
            h_all=np.zeros((xin_shape),dtype=float)
            u_dta=np.zeros((xin_shape),dtype=float)
            v_dta=np.zeros((xin_shape),dtype=float)
        
        
            if rotate2nri:
                import okean
                xin,yin = okean.calc.rot2d(xin,yin,ang=-np.pi/2.0)
                #u, v    = okean.calc.rot2d(u  ,v  ,ang=-np.pi/2.0)

            h_all[:]=interp3(xr,yr,h,xin,yin,'csa')
            u_all[:]=interp3(xr,yr,u,xin,yin,'csa')
            v_all[:]=interp3(xr,yr,v,xin,yin,'csa')
            z_all[:]=interp3(xr,yr,z,xin,yin,'csa')
            
            if real_data:
                u_dta[:] = interp3(x_sar[~m_sar],y_sar[~m_sar],u_sar[~m_sar],xin,yin,'csa')
                v_dta[:] = interp3(x_sar[~m_sar],y_sar[~m_sar],v_sar[~m_sar],xin,yin,'csa')
            
            all_trans.update({case+'_'+trans1:np.array([xin,yin,h_all,\
                            u_all,v_all,np.sqrt(u_all*u_all+v_all*v_all),\
                            u_dta,v_dta,np.sqrt(u_dta*u_dta+v_dta*v_dta),\
                            ])})   
    
    
    ### Part one
    ### start to plot pcolor 2D maps
    #var='h'
    for var in vars:
        print '   > ', var
        ## Aranging the axes
        if False:
            #fig, axe = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True,figsize=(7,9)) 
            fig, axe = plt.subplots(nrows=2, ncols=3,figsize=(7,9)) 
            axmap   = axe[0,1]
            axleft  = axe[0,0]
            axrigt  = axe[0,2]
            axbott  = axe[1,1]
            
            #make extra axes invisible
            axe[1,0].set_visible(False)
            axe[1,2].set_visible(False)
            
            # make some labels invisible
            #plt.setp(axmap.get_yticklabels() + axrigt.get_yticklabels() + axmap.get_xticklabels(),
            #         visible=False)
            #plt.setp(axleft.get_xticklabels() + axrigt.get_xticklabels() + axbott.get_yticklabels(),
            #         visible=True)
        else:    
            fig = plt.figure(figsize=(9,7))#, dpi=600)
            padx=0.1
            pady=0.1
            hm=0.5                       #height map
            wm=0.4                       #width map
            wl=(1-4*padx-wm)/2           #width left and right
            hb=1-hm-2*pady               #height bottom
            wb=1-2*padx                  #width bottom
            rect_left   = [padx             , 2*pady+hb  , wl   , hm-pady/2]
            rect_map    = [2*padx + wl      , 2*pady+hb  , wm   , hm-pady/2     ]
            rect_right  = [3*padx + wl + wm , 2*pady+hb  , wl   , hm-pady/2]
            rect_bottom = [padx             , pady/3     , wb   , hb     ]
        
            axleft  = plt.axes(rect_left)
            axmap   = plt.axes(rect_map)
            axrigt  = plt.axes(rect_right)
            axbott  = plt.axes(rect_bottom)
                
        axmap.set_aspect(1.)
        ##### Which Var to plot
        if var=='h':
            varu='[m]'
            color_lab='Depth'+varu
            val= h_cut
            valt=h_tru
            lim1= 0
            lim2= 9
            cmap=cmapj
            vec  = False
        
        if var=='h_dif':
            varu='[m]'
            color_lab='Depth - Depth_true'+varu
            #color_lab='Depth perturbations'+varu
            valt=h_tru
            val= h_cut-h_tru
            #lim1= -0.2
            #lim2=  0.2
            lim1= -1.5
            lim2=  1.5            
            cmap=cmapo
            vec  = False

        if var=='uv':
            varu='[ms$^{-1}$]'
            color_lab='Vel. '+varu
            if normalize: color_lab='normalized velocity [-] '
            val= np.sqrt(u_cut*u_cut+v_cut*v_cut)
            valt=np.sqrt(u_tru*u_tru+v_tru*v_tru)
            lim1= 0.2
            lim2= 0.6               
            if normalize:
                lim2= 1.0
            
            if True:
                #to mask out velocities smaller than
                min_mask=lim1
                val= np.ma.masked_where(val<min_mask,val)
                plt.cm.jet.set_bad(color='white', alpha=None)                
            
            cmap = cmapj
            vec  = True
                
        if var=='u':
            varu='[ms$^{-1}$]'
            color_lab='U '+varu
            val= u_cut
            valt=u_tru
            lim1= -1
            lim2=  1
            cmap=cmapo    
        
        if var=='v':
            varu='[ms$^{-1}$]'
            color_lab='V '+varu
            val= v_cut
            valt= v_tru
        
            lim1=  -1.0
            lim2=   0
            cmap=cmapj
        
        if var=='uv_dif':
            #prior stuf
            ma0b=h_tru.mask
            ma0=np.ones_like(h_tru)
            ma0[ma0b]=0
            #curent vel stuff
            ma1b=h_cut.mask
            ma1=np.ones_like(h_cut)
            ma1[ma1b]=0
            
            #general mask 
            mask01=ma0*ma1
            
            #apply mask to all
            u_tru=u_tru*mask01
            v_tru=v_tru*mask01
            u_cut =u_cut *mask01
            v_cut =v_cut *mask01                
            varu='[ms$^{-1}$]'
            color_lab='Vel. dif. '+varu
            
            dif_method='alex'
            if dif_method=='alex1':
                val= np.sqrt((u_cut-u_tru)**2+ (v_cut-v_tru)**2)
                #only for prior and true
                lim1=   0.0
                lim2=   0.25
                #lim2=   0.04
                min_mask=lim1                    
                
                #if itr!=true_itr:
                #    val= np.ma.masked_where(val<min_mask,val)
                #    np.cm.jet.set_bad(color='white', alpha=None)
                cmap=cmapj 
            else:
                val= np.sqrt(u_cut**2+v_cut**2) - np.sqrt(u_tru**2+v_tru**2)
                #lim1=  -0.05
                #lim2=   0.05
                lim1=  -0.25
                lim2=   0.25
                cmap=cmapo
        
        if normalize and var !='h':
            print '         val.max()  > ',val.max()
            print '         val.min()  > ',val.min()
            
            val=val/abs(val).max()
    
        xmin= x_cut.min()   
        xmax= x_cut.max() 
        ymin= y_cut.min()   
        ymax= y_cut.max()
        
        if vec:
            alpha=0.85
        else:
            alpha=1.0
            
        pcolor_plot=True
        if pcolor_plot:
            print '     >','create pcolor ...'
            img1=axmap.pcolor (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),cmap=cmap,alpha=alpha)
            img1.set_clim(lim1,lim2)    # evry thing more than 2 will be dark blue
        else:
            levels=pl.linspace(lim1,lim2,25)
            plt.contourf (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),levels,cmap=cmap)
        
        #matplotlib.rcParams['axes.formatter.limits'] = (-1,1)    #to cjange axes to show 1e-1
        descrip=ps.ordinal[1]+' '+varn[itr]
        axmap.text(xmin+0.0*(xmax-xmin),ymax-0.05*(ymax-ymin) ,
                 descrip ,size=9.0,color='k',backgroundcolor='w')      
        
        if var != 'h_dif' and var != 'h':   
            axmap.set_title(date_sar.isoformat())
            axmap.text(xmin+0.1*(xmax-xmin),ymin- 0.15*(ymax-ymin),base_dir.split('/')[-1] ,size=3.0,color='k',alpha=0.5)
        
        print '     >','create colorbar ...'
        ##### colorbar stuff ################################### 
        pos_ax=np.array (axmap.get_position())
        xsub1=pos_ax[0][0]
        xsub2=pos_ax[1][0]
        ysub1=pos_ax[0][1]            
        ysub2=pos_ax[1][1]
        #print xsub1,xsub2,ysub1,ysub2
        cbticks=np.linspace(lim1, lim2, 3, endpoint=True)
        cb=plt.colorbar(img1,ticks=cbticks,shrink=0.6,orientation='horizontal') #
        cb.ax.set_position([xsub1- 0.35*(xsub2-xsub1), ysub1 -0.2 *(ysub2-ysub1), 0.7, 0.015])  # col=1
        
        cb.set_label(color_lab)
        [i.set_linewidth(0.05) for i in cb.ax.spines.itervalues()]
    
           
        if plot_mask:
            print '     >','create mask ...'
            mask  = z_cut.mask
            maskp = np.ones_like(~mask)
            maskp = np.ma.masked_array(maskp,~mask)
            axmap.pcolor (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(maskp),cmap=plt.cm.bone,alpha=1.0)
           
        #############################################
        if vec:
             print '     >','create vec ...'
             acol=(0.6,0.7,0.3)
             #vectors
             units='inches'
             scale=10
             #paper               
             kx=5
             ky=5               
    
             axmap.quiver( x_cut[::ky, ::kx], y_cut[::ky, ::kx], u_cut[::ky, ::kx], v_cut[::ky, ::kx], \
                          pivot='tail', color='k', scale=scale)
             # vel reference
             xtext = xmin + 0.35 * (xmax-xmin)
             ytext = ymax - 0.1 * (ymax-ymin)
             axmap.quiver([xtext],[ytext],[0.0,],[0.5,],pivot='tail',scale=scale,color='w',width=0.003)
             axmap.quiver([xtext],[ytext],[0.5,],[0.0,],pivot='tail',scale=scale,color='w',width=0.003)
             axmap.text  ( xtext-50 ,  ytext-0.04*(ymax-ymin) ,'0.5 [ms$^{-1}$]',size=7,color='w')  
        
        if plot_contour:
             # contour stuff
             print '     >','create contour ...'
             lim1=0
             lim2=9
             levels=np.linspace(lim1,lim2,11)   #7)
             con1 = axmap.contour (np.squeeze(x_cut),np.squeeze(y_cut),h_cut,levels,colors='k',linewidths=0.4,alpha=0.75)
             #levels=pl.linspace(lim1,lim2,5)   #7)
             #con1 = pl.contour (np.squeeze(x_cut),np.squeeze(y_cut),h_tru,levels,colors='k',linewidths=0.4,alpha=1.0)
             plt.clabel(con1, inline=1, fontsize=3,fmt='%1.1f')        
        
        
        
        ##################################################          
        # the linewidth of the rectangular axis frame
        fr_linewidth=0.3
        [i.set_linewidth(fr_linewidth) for i in axmap.spines.itervalues()]
        
        axmap.set_ylim(ymin,ymax)
        axmap.set_xlim(xmin,xmax)
        
        dev1 = 500
        xmin_tic = ceil(xmin/dev1) * dev1
        xmax_tic = ceil(xmax/dev1) * dev1

        #xxticks = np.linspace(xmin_tic, xmax_tic, 5 , endpoint=False)
        xxticks = np.arange(xmin_tic, xmax_tic,dev1 )
        axmap.xaxis.set_ticks(ticks=xxticks)
        for tick in axmap.xaxis.get_major_ticks():
            tick.label.set_rotation(10) 
        
        #xxticks=np.linspace(xmin//100 *100, xmax//100*100, 4 , endpoint=False)+200
        #axmap.xaxis.set_ticks(ticks=xxticks) 
        
        axmap.set_ylabel('Y [m]')     
        axmap.set_xlabel('X [m]')
        
        #########################################
        #########################################
        ##### plot left and right and bottom
        trans_select = 'crA1'
        #trans_select = 'crA2'
        #trans_select = 'crA3'
        
        tru_key = 'tru_' + trans_select
        asi_key = 'asi_' + trans_select
        pri_key = 'pri_' + trans_select
        
        lwmap=1
        lw=1
    
        if   'h' in var:
            nnum=2
        elif 'uv'in var:
            nnum=5
        print '     >','create left ...'
        #plot left transect
        #plot A1-A! line on map
        xa1=all_trans[tru_key][0]
        ya1=all_trans[tru_key][1]
        axmap.plot(xa1,ya1,'k-D',lw=lwmap,label='Sec. A1-A1', markevery=20 , ms =5,mfc='none',markeredgewidth=lwmap)
        
        if not 'dif' in var:
            axleft.plot(all_trans[tru_key][nnum],ya1,'r-',lw=lw,label='true')
            axleft.plot(all_trans[asi_key][nnum],ya1,'b-',lw=lw,label='Assim. '+varn[itr])
            axleft.plot(all_trans[pri_key][nnum],ya1,'g-',lw=lw,label='prior')
        else:
            axleft.plot(-all_trans[tru_key][nnum]+all_trans[asi_key][nnum],ya1,'b-',lw=lw,label='assim-true')
            axleft.plot(-all_trans[tru_key][nnum]+all_trans[pri_key][nnum],ya1,'g-',lw=lw,label='prior-true')
        
        if real_data and var == 'uv':
            axleft.plot(all_trans[tru_key][8],ya1,'k-',lw=lw,label='Data')
        
        ############################################################################################################
        
        
        #axleft.legend()
        if   var=='h':
             lim1=0
             lim2=12.0
        elif var=='uv':
             lim1=0
             lim2=2.0  
        
        axleft.set_ylabel('Y [m]')
        axleft.set_xlabel(color_lab)
        #if 'dif' in var:
        axleft.set_xlim(lim1,lim2)
        xxticks=np.linspace(lim1, lim2, 3 , endpoint=True)
        axleft.xaxis.set_ticks(ticks=xxticks)
        #else:
        #    for tick in axleft.xaxis.get_major_ticks():
        #        tick.label.set_rotation(45)
        
        #add text
        descrip=ps.ordinal[0]+' Sec.\n    A1-A1'
        axleft.text(lim1+0.2*(lim2-lim1),340, descrip ,size=11.0,color='k')#,backgroundcolor='w')  
        
        print '     >','create right ...'
        #plot right transect
        #plot B-B line on map
        xa1=all_trans['tru_crBB'][0]
        ya1=all_trans['tru_crBB'][1]
        axmap.plot(xa1,ya1,'k-<',lw=lwmap,label='Sec. B-B' , markevery=5 , ms =5,mfc='none',markeredgewidth=lwmap)
        if not 'dif' in var:
            axrigt.plot(all_trans['tru_crBB'][nnum],ya1,'r-',lw=lw,label='true')
            axrigt.plot(all_trans['asi_crBB'][nnum],ya1,'b-',lw=lw,label='Assim. '+varn[itr])
            axrigt.plot(all_trans['pri_crBB'][nnum],ya1,'g-',lw=lw,label='prior')
        else:
            axrigt.plot(-all_trans['tru_crBB'][nnum]+all_trans['asi_crBB'][nnum],ya1,'b-',lw=lw,label='assim-true')
            axrigt.plot(-all_trans['tru_crBB'][nnum]+all_trans['pri_crBB'][nnum],ya1,'g-',lw=lw,label='prior-true')
        
        if real_data and var == 'uv':
            axrigt.plot(all_trans['tru_crBB'][8],ya1,'k-',lw=lw,label='Data')
        #axrigt.legend()
#         
        if   var=='h':
            lim1=0
            lim2=6.0
        elif var=='uv':
            lim1=0
            lim2=2.0   


        axrigt.set_ylabel('Y [m]')
        axrigt.set_xlabel(color_lab)
        #add text
        descrip=ps.ordinal[2]+' Sec.\n    B-B'
        axrigt.text(lim1+0.2*(lim2-lim1),-1410, descrip ,size=11.0,color='k')#,backgroundcolor='w')  
        
        #if 'dif' in var:
        axrigt.set_xlim(lim1,lim2)
        xxticks=np.linspace(lim1, lim2, 3 , endpoint=True)
        axrigt.xaxis.set_ticks(ticks=xxticks)
        #else:
        #    for tick in axrigt.xaxis.get_major_ticks():
        #        tick.label.set_rotation(45)
        
        print '     >','create bottom ...'
        #plot bottom transect
        #plot main line on map
        xa1=all_trans['tru_main'][0]
        ya1=all_trans['tru_main'][1]
        axmap.plot(xa1,ya1,'k-o',lw=lwmap,label='Main channel', markevery=5 , ms =5,mfc='none',markeredgewidth=lwmap)
        if not 'dif' in var:
            axbott.plot(xa1,all_trans['tru_main'][nnum],'r-',lw=lw,label='true')
            axbott.plot(xa1,all_trans['asi_main'][nnum],'b-',lw=lw,label='Assim. '+varn[itr])
            axbott.plot(xa1,all_trans['pri_main'][nnum],'g-',lw=lw,label='prior')
        else:
            axbott.plot(xa1,-all_trans['tru_main'][nnum]+all_trans['asi_main'][nnum],'b-',lw=lw,label='assim-true')
            axbott.plot(xa1,-all_trans['tru_main'][nnum]+all_trans['pri_main'][nnum],'g-',lw=lw,label='prior-true') 
        
        if real_data and var == 'uv':
            axbott.plot(xa1,all_trans['tru_main'][8],'k-',lw=lw,label='Data')
        #axrigt.legend()            
            
        if   var=='h':
            lim1=0
            lim2=12.0
        elif var=='uv':
            lim1=0
            lim2=2.0  

        #if not 'dif' in var:
        
        axbott.set_ylim(lim1,lim2)
        axbott.set_xlim(-1400,1000)
        axbott.set_xlabel('X [m]')
        axbott.set_ylabel(color_lab)
        axbott.legend(ncol=3,frameon=False,loc='lower center',prop=dict(size='medium'))
        
        #add text
        descrip=ps.ordinal[3]+' Main channel'
        if  var=='h':
            axbott.invert_yaxis()
            ytext2=1
        elif var=='uv':
            ytext2=0.7
        else:
            ytext2=lim2* 0.8
            
      
        axbott.text(500,ytext2, descrip ,size=11.0,color='k')#,backgroundcolor='w')  
        axmap.legend(ncol=1,frameon=True,loc='lower left',prop=dict(size='small'))
        
        
        outdir=string.join(files[ntr[0]].split('/')[0:-3],'/')+'/post/transects/out_transect/'
        os.system('mkdir -p '+outdir)
        for ftype in ftypes:
            outfile=outdir+'/'+str(itr+1000)+'maps_'+varn[itr]+'_'+var+'.'+ftype
            
            if ftype=='png':
               plt.savefig(outfile,dpi=600)
            else:
               plt.savefig(outfile)    
         
        plt.close()
    
    
    #Back_up scr
    args=sys.argv
    scr_name=args[0]    
    command='mkdir -p '+ outdir+'/scr/'
    os.system(command)
    scr_dir1=os.getcwd()
    os.system('cp -fr  '+scr_dir1+'/'+scr_name +'    '+outdir+'/scr')
    os.system('cp -fr  '+scr_dir1+'/../transect/data1.py'+'    '+outdir+'/scr')

##############################################################################
##############################################################################
##############################################################################
