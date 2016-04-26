#!/usr/bin/env python
# -*- coding: utf-8 -*-

#### Map plot for assimilation  outputs ####
#############################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0 03/25/2013 02:14:41 PM    prepared for general use
#
#

#
# <nbformat>2</nbformat>

# <markdowncell>

# <script type="text/javascript">
#     $('div.input').show();  
# </script>
# 
# 
# Hide Code  >> to show >> $('div.input').show(); 
# 
# 5 Mar 2013

# <markdowncell>

# 1. Experiment
# ==============
# We run 200 ROMS members for 4 days and the same members for SWAN at one stationary time step. The results for assimilation of diffrent combination of data will be presented.
# 
# ## Data sources:
# 
# *    3 days of swift data 08-10 may 2012
# *    2 velocity field derived from SAR images (10th may 2012)
# *    1 set of wave f-k data derived from Radar images by cBathy
# 
# 
# 
# 2. Real data vs. twin test assimilation 
# ===============
# It seems that bottom roughness could have some effects on bathymetry inversion. 
# To test this idea range of diffrent Quadratic bottom drag coefficient (rdrg2) have been tested.
# The results for min and max of this range for rdrg2=0.005 and rdrg2=0.04 presnted hereafter:
# 
# ##Assimilation constants:
# ###  Error range is very important and sensitive parameter. 
# 
# For this experiment:
# 
# * SAR velocity error       = 0.1 ms$^{-1}$   
# * Radar wave-number error    from cBathy out put
# 
# * Twin velocity error      = 0.05 ms$^{-1}$
# * Twin wave-number error    = 0.02 m$^{-1}$
# 
# 
# My guess of SAR velocity error for the very first tests not used in this experiment:  $$\quad ERR_{ij}=0.5  abs \left(1-\frac{Max(SWIFT)}{Max(SAR)}\right ) * SAR_{ij}$$ 

# <markdowncell>

# 2. SWIFT data test
# ===============
# 
# Here we like to test how swift could help the assimilation process.
# 
# 
# 
# ### For now:
# * SAR very rough error estimatimation:  $$\quad ERR_{ij}=0.5  abs \left(1-\frac{Max(SWIFT)}{Max(SAR)}\right ) * SAR_{ij}$$
# * f-k provided by cBathy 
# * SWIFT ERRxyt=0.2 m/s

# <codecell>

# Where to create figures
#cur_dir='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/set4_all/post/02-maps/5mar2012_assimilation_test'
#%cd $cur_dir
#!mkdir -p 5mar2012_assimilation_test/pic
#%cd 5mar2012_assimilation_test
#tmp=!pwd
#cur_dir=tmp.pop()

# <markdowncell>

# Functions

# <codecell>

global rotate2nri
#### Functions ####
def read_nc_file(itr,date_assim,surface_vel):
    #Read in his_file
    base=base_dir+'/run_'+str(1000+itr)+'/08_forward/'
    #ncf=base+'nri_his.nc'
    ncf=base+'nri_his_0001.nc'

    print itr, ncf
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
    n=pl.array(ind).item()
    #print '     >','chosen time step from netCDF file > ',n
    
    #n=142  # ebb
    #for one day run
    #n=42  # ebb
    #n=29   # flood
    
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
    u=pl.zeros_like(z)
    v=pl.zeros_like(z)
    
    #interpolate to rho points (aprox.)
    u[:,1:-1]=(ubar[:,:-1]+ubar[:,1:])/2.
    u[:,0]   = ubar[:,0]
    u[:,-1]  = ubar[:,-1]
    
    v[1:-1,:]=(vbar[:-1,:]+vbar[1:,:])/2.
    v[0   ,:]= vbar[0  ,:]
    v[-1  ,:]= vbar[-1 ,:]   
    
    m_cut=maskr [j1:j2:k,i1:i2:k]
    x_cut=xr    [j1:j2:k,i1:i2:k]
    y_cut=yr    [j1:j2:k,i1:i2:k]
    h_cut=h     [j1:j2:k,i1:i2:k]
    u_cut=u     [j1:j2:k,i1:i2:k]
    v_cut=v     [j1:j2:k,i1:i2:k]
    z_cut=z     [j1:j2:k,i1:i2:k]
    
    #x_cut=pl.ma.masked_where(m_cut==0,x_cut)
    #y_cut=pl.ma.masked_where(m_cut==0,y_cut)
    h_cut=pl.ma.masked_where(m_cut==0,h_cut)
    u_cut=pl.ma.masked_where(m_cut==0,u_cut)
    v_cut=pl.ma.masked_where(m_cut==0,v_cut)
    z_cut=pl.ma.masked_where(m_cut==0,z_cut)
    nc.close()
    
    if rotate2nri:
        import okean
        x_cut,y_cut=okean.calc.rot2d(x_cut,y_cut,ang=-np.pi/2.0)
        u_cut,v_cut=okean.calc.rot2d(u_cut,v_cut,ang=-np.pi/2.0)

    return x_cut,y_cut,h_cut,u_cut,v_cut,z_cut


def read_sar_file(time_step):
    #slice information        
    #Read in his_file
    #ncf=inp_dir+'/obs/sar/20120510_182231_stitched_ks.nc'
    #ncf='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run7_real_paper/02_ebb_1-iterate_01/inp/obs/sar/20120510_182231_stitched_ks.nc'
    #ncf=inp_dir+'/obs/sar/uASAR_05102012_smoothed4test_21feb2013.nc'
    #ncf='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run7_real_paper/02_ebb_1-iterate_01/inp/obs/sar/back/20120510_182231_stitched.nc'
    #ncf='/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run7_real_paper/02_ebb_1-iterate_01/inp/obs/sar/uASAR.nc'
    #print ncf
    #if copy1: os.system('cp '+ncf+'   nri_his_'+str(1000+itr)+'.nc')
    ncf   = sar_file
    nc    = netCDF4.Dataset(ncf)
    ncvar = nc.variables
    timename_sim='time'
    utime_sim  = netcdftime.utime(ncvar[timename_sim].units)
    sdates_sim = utime_sim.num2date(ncvar[timename_sim][:]) 
    nt    = len(sdates_sim)
    date  = sdates_sim[time_step]
    n     = time_step   # ebb
     
    xr = ncvar['x'][:]
    yr = ncvar['y'][:]
    u  = ncvar['u'][n,:]
    v  = ncvar['v'][n,:]
    
    x_cut = xr    [j1:j2:k,i1:i2:k]
    y_cut = yr    [j1:j2:k,i1:i2:k]
    u_cut = u     [j1:j2:k,i1:i2:k]
    v_cut = v     [j1:j2:k,i1:i2:k]
   
    u_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)| np.ma.is_masked(u_cut)),u_cut )
    v_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)| np.ma.is_masked(v_cut)),v_cut )
    m_cut = u_cut.mask

    #x_cut=pl.ma.masked_where(m_cut,x_cut)
    #y_cut=pl.ma.masked_where(m_cut,y_cut)
    u_cut=pl.ma.masked_where(m_cut,u_cut)
    v_cut=pl.ma.masked_where(m_cut,v_cut)
    nc.close()
    
    if rotate2nri:
        import okean
        x_cut,y_cut=okean.calc.rot2d(x_cut,y_cut,ang=-np.pi/2.0)
        u_cut,v_cut=okean.calc.rot2d(u_cut,v_cut,ang=-np.pi/2.0)
    
    
    return x_cut,y_cut,u_cut,v_cut,date

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


# <markdowncell>

# Plot Maps code

# <codecell>

import numpy as np
import pylab as pl
import netCDF4
import netcdftime
import datetime as datetime
import os,sys
import string
from matplotlib.ticker import OldScalarFormatter
sys.path.append('/home/server/pi/homes/moghimi/Desktop/cowast/projects/my_py_libs/plot')
import plot_settings as ps
pl.close('all')

import math

########  read directories  ##################
sys.path.append('../../pysim')
import base_info
global base_dir
global inp_dir
global sar_file

time_now=datetime.datetime.now().isoformat()[5:16]

# base directory for this set or iterations
base_dir=base_info.base_dir
inp_dir=base_info.inp_dir
prior=base_info.prior
   
global i1,i2,j1,j2,k
#slice information        

#for paper_whole
if True:
    i1=20
    i2=180
    j1=20
    j2=180

#for paper_whole #2
if False:
    i1=60
    i2=140
    j1=110
    j2=190

#just for uv_dif
#i1=5
#i2=195
#j1=5
#j2=195

#for subset
# i1=70
# i2=120
# j1=120
# j2=170

xlen=(i2-i1+1)*20
ylen=(j2-j1+1)*20
print 'xlen plot= ', xlen
print 'ylen plot= ', ylen
print ' ylen/xlen=', ylen *1.0 /xlen

k  = 1
mm = 1
#mm =  2
#########
#var to plot
#all_var=['uv']  #,'uv_dif','v','u','h_dif']
#all_var=['h_dif']
#all_var=['h','uv']
#all_var=['h','uv']


all_var=['h']
#all_var=['h_dif','uv_dif']
#all_var=['h','uv','h_dif']
#all_var=['uv']

##plot SAR velocity field
rotate2nri        = True
animation         = True
plot_sar          = False
plot_vel          = True
plot_data_points  = False
normalize         = False
plot_mask         = True
plot_contour      = True
plot_rect         = False
plot_true         = False
calc_stat         = True
ebb_data          = False
xtra_labl         = True
surface_vel       = True
ftypes=['png','pdf']
#ftypes=['png']

if True:
    #SYN paper
    #please give this data
    # number of iterations
    #ntr=[0,1,2,3,4,5]    # itration to be ploted
    #ntr=[0,1,2,4,6,7]    # itration to be ploted
    ntr=[3,4,5]          # iwaves
    #ntr=[0]    # itration to be ploted
    true_itr=8   #  itr number for true state
    base_dir='/home/shusin5/users/moghimi/assimilation/assim_local/wave_current/run6_syn_using_3D_runs/Compare_3d_4_paper/'
    lab1=['Prior',
          'CUR-ebb',
          'CUR-flood',
          'WAV-Tm01=4s',
          'WAV-Tm01=7s',
          'WAV-Tm01=14s',
          'CW-ebb',
          'CW-flood',
          'True']



#prepare output directory

vars_join=string.join(all_var)
vars_join=vars_join.replace(' ','-')
xlab1=[]
for i in range(len(ntr)): 
    xlab1.append(str(ntr[i]))
xlab1.append(str(true_itr))
xlab1=string.join(xlab1)
xlab1=xlab1.replace(' ','-')
outdir=base_dir+'/post/maps/'+vars_join+'_'+xlab1
command='mkdir -p '+ outdir+'/scr/'
os.system(command)
#from octant import plotting as ocpl
#cmap=ocpl.cmap_brightened(ocpl.cmap_discretize(pl.cm.jet,5),0.6)
##cmap=ocpl.cmap_brightened(ocpl.jetWoGn(),0.6)
cmapj  = pl.cm.jet_r
#cmapj = pl.cm.jet
cmapo=ps.jetWoGn(reverse=True)
#cmapo = ps.my_cmap
#cmap = pl.cm.spectral
#cmap = pl.cm.gist_earth

#cmap=pl.cm.Paired
##############
if plot_true:
    ntr.append(true_itr)

# to take into account plotting SAR velocity field

arg=sys.argv
#SAR DATAinp_dir
if plot_sar:
    sar_file=base_dir+'/inp/obs/sar/uASAR_full.nc'
    xs,ys,us,vs,date_sar=read_sar_file(time_step=0)
elif len(arg)>1:
    print arg[1] 
    import dateutil.parser
    date_str=arg[1]
    date_sar = dateutil.parser.parse(date_str)
else:
    #assim time
    if ebb_data:
        date_sar  = datetime.datetime(2012,05,10,20)     #ebb
    else:
        #date_sar = datetime.datetime(2012,05,11,1,30)   #flood
        date_sar = datetime.datetime(2012,05,11,2,30)   #flood

print '** > > plot date > ', date_sar

#plot mask
xm,ym,hm,um,vm,zm=read_nc_file(1,date_sar,surface_vel)
maskplot=hm.mask
maskpl=np.zeros(maskplot.shape)
maskpl=np.ma.masked_array(maskpl,~maskplot)


################ figure setting stuff >>>>>
ps.set_mode("m")
fig_size1=ps.get_figsize(500)

linewidth=0.2
#====== subplot adjustments ===============
left1  = 0.15 # the left side of the subplots of the figure
right1 = 0.95    # the right side of the subplots of the figure
bottom1= 0.3  # the bottom of the subplots of the figure   (ntr==16   bottom=0.05)
top1   = 0.95      # the top of the subplots of the figure
wspace1= 0.01   # the amount of width reserved for blank space between subplots
hspace1= 0.01   # the amount of height reserved for white space between subplots
###############################################<<<<<<

#subplot params
icount=1
rrow=2
# figs param
dpi=600
fwidth1,fheight1= fig_size1
my_fs=3
fontsize=6

if True:
    for var in all_var:
        #to plot vectors
        vec=False
        if var=='uv':
            vec=True
            #vec=False
        print '***************************'
        print '** > > outputs options >>'
        print 'Var              >', var
        print 'vector plot      >',vec
        print 'Rotae to NRI     >', rotate2nri
        print 'SAR plot         >', plot_sar
        print 'plot_data_points >',plot_data_points
        print 'normalize        >',normalize
        print 'plot_mask        >',plot_mask
        print 'plot_contour     >',plot_contour
        print 'plot_rect        >',plot_rect
        print 'plot_true        >',plot_true
        print 'calc_stat        >',calc_stat
        print 'outputs          >',ftypes
        print 'number of panels >',len(ntr)
        print '*****************************'
       
        
        if plot_sar and not var=='h': 
            ntr.append(true_itr+1)        
        
        nall=len(ntr)
        icol= nall

        #only 1 graph
        #fheight=fheight * nall  * 1.5  
        #fwidth=fwidth   * icol  * 1.5  
        #####################
        #icol = 1
        #fheight=fheight1 * 1.2 
        #fwidth =fheight1 * 1.2 
        ##########################
        #icol=2
        #fheight=fheight1    
        #fwidth= fheight1    * nall * 0.92
        #####################################
        #icol=3
        fheight=fheight1    
        fwidth= fheight1    * nall *0.92
        ##############
        # 4 col  and 2 row
        #icol= nall/2.0
        #fheight=fheight1    
        #fwidth= fheight1  * nall * 0.9 
        ###############
        # 4 col 
        #fheight=fheight1    * 1.0
        #fwidth= fheight1    * nall *0.84
        #bottom1= 0.3   # the bottom of the subplots of the figure   (ntr==16   bottom=0.05)
        ##############
        # 5 col 
        #fheight=fheight1  *  1.0   * 1.2 *0.8
        #fwidth =fheight1  * nall   * 1.0 *0.8
        ###############

        nn=0
        f=pl.figure(figsize= (fwidth,fheight))   # for 2 x 2 plot
        for itr in ntr:
            ######## SUBPLOT ARRANGMENTS    #########################
            nn+=1
            irow= int(nall/icol)
               
            print  '     >',irow,icol,nn
            axp=f.add_subplot(irow,icol,nn) # (#rows,#columns,subplot#)
            pl.subplots_adjust(left=left1, bottom=bottom1, right=right1, top=top1,
                 wspace=wspace1, hspace=hspace1)
            axp.set_aspect(1.0)
            
            pos_ax1=axp.get_position ()
            pos_ax=np.array (axp.get_position ())
            #print 'y',pos_ax[0][1], pos_ax[1][1]
            #print 'x',pos_ax[0][0], pos_ax[1][0]
            #print '-----------------------------'
            aheight=pos_ax[1][1] -pos_ax[0][1]
            awidth= pos_ax[1][0] -pos_ax[0][0]
            ##########################################################
            if xtra_labl:
                varn=lab1[itr]
            else:
                varn= 'iteration '+str(itr)  
            
            #if itr==0:         varn= 'Prior' 
            #if itr==true_itr:  varn= 'True' 
            
            #read netcdf file
            if itr!=true_itr+1:
                x_cut,y_cut,h_cut,u_cut,v_cut,z_cut=read_nc_file(itr,date_sar,surface_vel)
                u_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)) ,u_cut )
                v_cut = np.ma.masked_where(((np.abs(u_cut) > 5)|(np.abs(v_cut) > 5)) ,v_cut )
            else:
                print '     >','SAR >>>>>>>'
                x_cut,y_cut,u_cut,v_cut,date_sar = read_sar_file(time_step=0)
                varn = 'SAR data'
        
            rmaskp=u_cut.mask
            uv=pl.sqrt(u_cut**2+v_cut**2)
         
			#True value
            x_tru,y_tru,h_tru,u_tru,v_tru,z_tru=read_nc_file(true_itr,date_sar,surface_vel)
            
            if var=='h':
                varu='[m]'
                color_lab='Depth'+varu
                val  = h_cut
                valt = h_tru
                lim1= 0
                lim2= 9
                cmap=cmapj

            if var=='h_dif':
                varu='[m]'
                
                color_lab='Depth - Depth_true'+varu
                #color_lab='Depth perturbations'+varu
                valt = h_tru
                val  = h_cut-h_tru
                lim1= -1
                lim2=  1
                cmap=cmapo
        
            if var   =='uv':
                varu ='[ms$^{-1}$]'
                color_lab='Vel. '+varu
                if normalize: color_lab='normalized velocity [-] '
                val  = pl.sqrt(u_cut*u_cut+v_cut*v_cut)
                valt = pl.sqrt(u_tru*u_tru+v_tru*v_tru)
                if plot_sar:
                    valt = pl.sqrt(us*us+vs*vs)
                else:
                    valt = pl.sqrt(u_tru*u_tru+v_tru*v_tru)
                
                lim1 = 0.2
                #lim2 = 0.6               
                lim2 = 1.0               

                if normalize:
                    lim2= 1.0

                    u_cut = u_cut / abs(val).max()
                    v_cut = v_cut / abs(val).max()
                    u_tru = u_tru / abs(valt).max()
                    v_tru = v_tru / abs(valt).max()
                    val  = val  /abs(val).max()
                    valt = valt /abs(valt).max()                    
                    
                    
                if True:
                    #to mask out velocities smaller than
                    min_mask=lim1
                    val= np.ma.masked_where(val<min_mask,val)
                    pl.cm.jet.set_bad(color='white', alpha=None)                
                
                cmap=cmapj
            if var=='u':
                varu='[ms$^{-1}$]'
                color_lab='U '+varu
                val  = u_cut
                if plot_sar:
                    valt = us
                else:
                    valt = u_tru
                    
                lim1= -1
                lim2=  1
                cmap=cmapo    
        
            if var=='v':
                varu='[ms$^{-1}$]'
                color_lab='V '+varu
                val= v_cut
                if plot_sar:
                    valt = vs
                else:
                    valt = v_tru
                lim1=  -1.0
                lim2=   0
                cmap=cmapj
            
            if var=='uv_dif':
                #prior stuf
                x_tru,y_tru,h_tru,u_tru,v_tru,z_tru=read_nc_file(true_itr,date_sar,surface_vel)
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
                varn= lab1[itr]+ '-True'
                if ebb_data:
                    varn = varn+' (Ebb)'
                else:
                    varn = varn+' (Flood)'
                varu='[ms$^{-1}$]'
                color_lab='Vel. dif. '+varu
                
                dif_method='alex'
                if dif_method=='alex':
                    val= pl.sqrt((u_cut-u_tru)**2+ (v_cut-v_tru)**2)
                    #orig
                    lim1=   0.04
                    lim2=   0.1
                    min_mask=0.04
                    #only for prior and true
                    lim1=   0.1
                    lim2=   0.5
                    min_mask=lim1                    
                    
                    if itr!=true_itr:
                        val= np.ma.masked_where(val<min_mask,val)
                        pl.cm.jet.set_bad(color='white', alpha=None)
                    cmap=cmapj 
                else:
                    val= pl.sqrt(u_cut**2+v_cut**2) - pl.sqrt(u_tru**2+v_tru**2)
                    lim1=  -0.08
                    lim2=   0.08
                    cmap=cmapo

            
            val = np.ma.masked_where(np.isnan(val),val)
            if normalize and var !='h' and var !='uv':
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
                img1=pl.pcolor (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),cmap=cmap,alpha=alpha)
                pl.clim(lim1,lim2)    # evry thing more than 2 will be dark blue
            else:
                levels=pl.linspace(lim1,lim2,25)
                pl.contourf (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),levels,cmap=cmap)
            #color_lab=varn + varu
            
            
            if plot_data_points and itr>1:
                kk=2
                xdata=x_cut  [::kk,::kk]
                ydata=y_cut  [::kk,::kk]
                pl.plot(xdata.flatten(),ydata.flatten(),'k+',markersize=1, alpha=0.4,markeredgecolor='None')#,markerfacecolor='None')
            
                       #width=width,headwidth=headwidth,headlength=headlength ,scale=scale)
         
            #matplotlib.rcParams['axes.formatter.limits'] = (-1,1)    #to cjange axes to show 1e-1
            descrip=ps.ordinal[mm+nn-1]+' '+varn
            #descrip=str(itr)+') '+varn
            #pl.text(xmin+2.5,ymax-2 , ps.ordinal[itr] ,size=12.0)      
            pl.text(xmin+0.0*(xmax-xmin),ymax-0.05*(ymax-ymin) ,
                     descrip ,size=9.0,color='k',backgroundcolor='w')      

            if nn==1 and var != 'h_dif' and var != 'h':   
                pl.title(date_sar.isoformat())
                pl.text(xmin+0.1*(xmax-xmin),ymin- 0.25*(ymax-ymin),base_dir.split('/')[-1] ,size=3.0,color='k',alpha=0.5)

            if plot_mask:
                print '     >','create mask ...'
                pl.pcolor (np.squeeze(xm),np.squeeze(ym),np.squeeze(maskpl),cmap=pl.cm.bone,alpha=1.0)
            
            #title1= dirn+'   '+str(utime.num2date(ncv['time'][ntime]))
            #title(title1,fontsize=11)
            ##################################################          
            # the linewidth of the rectangular axis frame
            ax=pl.gca()
            fr_linewidth=0.3
            [i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
            
            pl.ylim(ymin,ymax)
            pl.xlim(xmin,xmax)
            
            dev1 = 500
            xmin_tic = math.ceil(xmin/dev1) * dev1
            xmax_tic = math.ceil(xmax/dev1) * dev1

            #xxticks = np.linspace(xmin_tic, xmax_tic, 5 , endpoint=False)
            xxticks = np.arange(xmin_tic, xmax_tic,dev1 )
            ax.xaxis.set_ticks(ticks=xxticks)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(10) 
            
            pl.ylabel('Y [m]')     
            pl.xlabel('X [m]')
            
            if not ax.is_last_row():
               pl.setp( ax, 'xticklabels', [] )
               pl.xlabel('')   
            
            if not ax.is_first_col():
               pl.setp( ax, 'yticklabels', [] )
               pl.ylabel('')   
            
            
            #subplots_adjust(bottom=bottom1)
            
            
            pos_ax=np.array (axp.get_position ())
        #    print '>>>>-----------------------------'
        #    print 'y',pos_ax[0][1], pos_ax[1][1]
        #    print 'x',pos_ax[0][0], pos_ax[1][0]
        #    print '>>>>-----------------------------'
            
            pl.subplots_adjust(left=left1, bottom=bottom1, right=right1, top=top1,
                 wspace=wspace1, hspace=hspace1)
            
            axp.set_position (pos_ax1)
        
        
            
            acol=(0.6,0.7,0.3)
            #vectors
            units='inches'
            scale=10
            if vec:
                print '     >','create vec ...'
                #paper               
                kx=5
                ky=5               
                #subset
                #kx=1
                #ky=1                       
                
                #if not itr==true_itr+1:
                pl.quiver( x_cut[::ky, ::kx], y_cut[::ky, ::kx], u_cut[::ky, ::kx], v_cut[::ky, ::kx], pivot='tail', color='k', scale=scale)
                # vel reference
                xtext= xmin+0.05*(xmax-xmin)
                ytext= ymax-0.2*(ymax-ymin)
                if rotate2nri:
                    ytext= ymax-0.8*(ymax-ymin)
                
                if not normalize:
                    pl.quiver([xtext],[ytext],[0.0,],[0.5,],pivot='tail',scale=scale,color='w',width=0.003)
                    pl.quiver([xtext],[ytext],[0.5,],[0.0,],pivot='tail',scale=scale,color='w',width=0.003)
                    pl.text  ( xtext-50 ,  ytext-0.04*(ymax-ymin) ,'0.5 [ms$^{-1}$]',size=7,color='w')  

            if plot_contour:
                # contour stuff
                print '     >','create contour ...'
                lim11=0
                lim21=9
                levels=pl.linspace(lim11,lim21,11)   #7)
                con1 = pl.contour (np.squeeze(x_cut),np.squeeze(y_cut),h_cut,levels,colors='k',linewidths=0.4,alpha=0.75)
                
                #levels=pl.linspace(lim1,lim2,5)   #7)
                #con1 = pl.contour (np.squeeze(x_cut),np.squeeze(y_cut),h_tru,levels,colors='k',linewidths=0.4,alpha=1.0)

                pl.clabel(con1, inline=1, fontsize=3,fmt='%1.1f')
            
            if plot_rect:
                
                #Wave_window
                x1=-1200.0
                x2= 1200.0
                y1=-1200.0
                y2= 540.0
                
                px=np.array([x1,x2,x2,x1,x1])
                py=np.array([y1,y1,y2,y2,y1])
                if rotate2nri:
                    import okean
                    px,py=okean.calc.rot2d(px,py,ang=-np.pi/2.0)
                
                #pl.plot(px,py,'r-',lw=2.5)                               
                
                #cur_window
                x1=-500.0
                x2=750.0
                y1=-600.0
                y2=1700.0
                px=np.array([x1,x2,x2,x1,x1])
                py=np.array([y1,y1,y2,y2,y1])
                if rotate2nri:
                    import okean
                    px,py=okean.calc.rot2d(px,py,ang=-np.pi/2.0)
                
                pl.plot(px,py,'g-',lw=2.5)

                #assim_window
                ix=20,180  # from i=60 to i=150
                jy=10,190
                x1=-1450.0
                x2= 1750.0
                y1=-1740.0
                y2= 1850.0                

                px=np.array([x1,x2,x2,x1,x1])
                py=np.array([y1,y1,y2,y2,y1])
                if rotate2nri:
                    import okean
                    px,py=okean.calc.rot2d(px,py,ang=-np.pi/2.0)
                #pl.plot(px,py,'b-',lw=2.5)
            
            vstatus= True
            try:
                valt
            except NameError:
                vstatus= False
            
            if calc_stat and vstatus and var != 'uv_dif' and var != 'h_dif' and lab1[itr] !='SAR':
                #if (var == 'h' and  not 'True' in lab1[itr]) or var != 'h':
                if (not 'True' in lab1[itr]):

                    print '     > calc statistics ...'
                    mask = (h_tru.mask  | val.mask | valt.mask | np.isnan(valt+val))
    
                    da = valt[~mask].data
                    mo = val [~mask].data
                   
                    bi,rm=statatistics(da,mo)	
                    
                    bias_txt='Bias = %.3f %s' % (bi,varu)
                    rmse_txt='RMSE = %.3f %s' % (rm,varu)
                    
                    xtext= xmin+0.0*(xmax-xmin)
                    ytext= ymax-0.1*(ymax-ymin)
                    pl.text  ( xtext ,   ytext                  ,bias_txt,size=6,color='k',backgroundcolor='w')    
                    pl.text  ( xtext ,   ytext-0.05*(ymax-ymin) ,rmse_txt,size=6,color='k',backgroundcolor='w')    

   
        
            if nn == nall//2+1:   
               print '     >','create colorbar ...',nn
               ##### colorbar stuff ################################### 
               xsub1=pos_ax[0][0]
               xsub2=pos_ax[1][0]
               ysub1=pos_ax[0][1]            
               ysub2=pos_ax[1][1]
               #print xsub1,xsub2,ysub1,ysub2
               cbticks=np.linspace(lim1, lim2, 3, endpoint=True)
               #cbax = f.add_axes([xsub1-0.3*(xsub2-xsub1), ysub1 - 0.3 *(ysub2-ysub1), 0.25, 0.03]) #2panel 
               #cbax = f.add_axes([xsub1+0.25*(xsub2-xsub1), ysub1 - 0.3 *(ysub2-ysub1), 0.15, 0.03]) #multi 
               #cbax  = f.add_axes([xsub1+0.35*(xsub2-xsub1), ysub1 - 0.3 *(ysub2-ysub1), 0.3, 0.02]) #1panel 
               cbax  = f.add_axes([xsub1+0.12*(xsub2-xsub1), ysub1 - 0.3 *(ysub2-ysub1), 0.2, 0.03]) #3panel 

               cb=pl.colorbar(img1,cax=cbax,ticks=cbticks,format='%1.3g',orientation='horizontal') #
               
               cb.set_label(color_lab)
               [i.set_linewidth(0.05) for i in cb.ax.spines.itervalues()]
               #############################################
        
        #        if plot_sar:
        #            x_cut,y_cut,u_cut,v_cut,=read_sar_file(time_step=0)
        #            pl.quiver( x_cut[::ky, ::kx], y_cut[::ky, ::kx], u_cut[::ky, ::kx], v_cut[::ky, ::kx], pivot='tail', color='r'scale=scale)
        #        else:
        #            #read truth solution
        #            xe1,ye1,he1,ue1,ve1,ze1=read_nc_file(itr=true_itr)
        #            pl.quiver( xe1[::ky, ::kx], ye1[::ky, ::kx], ue1[::ky, ::kx], ve1[::ky, ::kx], pivot='tail', color='r',scale=scale)#,alpha=0.5)
        #            #plot current itration on top
            
        print '>> generate fig file'    
        #dir_in='.'    
        #outdir=dir_in+'/pic'
        #command='mkdir -p '+ outdir
        #os.system(command)
        #outfile=pl+'_cross_section_'+dirn+'_'+str(1000+ntime)+'.'+ftype
        #show()
        for ftype in ftypes:
            #outfile=outdir+'/maps_'+var+str(1000+len(ntr))+'_'+time_now+'.'+ftype
            #outfile=outdir+'/maps_'+var+str(1000+len(ntr))+'_'+date_str+'.'+ftype
            outfile=outdir+'/maps_'+var+str(1000+len(ntr))+'.'+ftype
            
            if ftype=='png':
               pl.savefig(outfile,dpi=dpi)
            else:
               pl.savefig(outfile)    
         
        pl.close()

#        
#        txtout= outdir+'/statistics'+str(itr+1000)+'.txt'
#        print txtout
#        fout=open(txtout,'w')
#        fout.write('statistics out for \n') 
#        st1='var   '+'bias  '+'rmse'+'\n'
#        st2='u  '  +str(biasu) +'  ' +str(rmseu) +'\n'
#        st3='v  '  +str(biasv) +'  ' +str(rmsev) +'\n'
#        st4='uv '  +str(biasuv)+'  ' +str(rmseuv)+'\n'
#        st5='h  '  +str(biash) +'  ' +str(rmseh) +'\n'
#        fout.write(st1)
#        fout.write(st2)
#        fout.write(st3)
#        fout.write(st4)
#        fout.write(st5)
#        fout.close()


scatter = False
if scatter:
    print '>>>   SCATTER SAR wirh TRUE and ASSIM  >>>'  
    ntr=[0,6,7]    # itration to be ploted
    nall = icol = len(ntr)
    
    lim1=1
    lim2=7
    #SAR
    print 'SAR >>>>>>>'
    xs,ys,us,vs,date_sar = read_sar_file(time_step=0)    
    uvs  = np.sqrt(us*us+vs*vs)
    mask = us.mask
    
    kj = 1
    
    fwidth,fheight= 10,4
    #f=pl.figure(figsize= (fwidth,fheight))   # for 2 x 2 plot
    #====== subplot adjustments ===============
    left1  = 0.1 # the left side of the subplots of the figure
    right1 = 0.9    # the right side of the subplots of the figure
    bottom1= 0.3   # the bottom of the subplots of the figure   (ntr==16   bottom=0.05)
    top1   = 0.9      # the top of the subplots of the figure
    wspace1= 0.05   # the amount of width reserved for blank space between subplots
    hspace1= 0.2   # the amount of height reserved for white space between subplots
    ###############################################<<<<<<
    
    print 'Truth >>>>>>>'
    xt,yt,ht,ut,vt,zt=read_nc_file(true_itr,date_sar,surface_vel)
    uvt=np.sqrt(ut*ut+vt*vt)
    mask = (mask & np.isnan(uvt) & uvt.mask)
    
    uvsf = uvs[~mask].flatten().data
    uvtf = uvt[~mask].flatten().data
    htf  = ht [~mask].flatten().data
    
    nn = -1
    
    irow= int(nall//icol)
    f,axgrid = pl.subplots(nrows=irow, ncols=icol, sharex=True, sharey=True,
                    figsize=(fwidth, fheight),facecolor='w', edgecolor='k',dpi=dpi )
    axgrid = np.array(axgrid)
    axgrid = axgrid.reshape(icol*irow)
    
    for itr in ntr:
        ######## SUBPLOT ARRANGMENTS    #########################
        nn+=1
        print  '  >>> ',nn
        axp=axgrid[nn] 
    
        pos_ax1=axp.get_position ()
        pos_ax=np.array (axp.get_position ())
        #print 'y',pos_ax[0][1], pos_ax[1][1]
        #print 'x',pos_ax[0][0], pos_ax[1][0]
        #print '-----------------------------'
        aheight=pos_ax[1][1] -pos_ax[0][1]
        awidth= pos_ax[1][0] -pos_ax[0][0]
        ##########################################################
        if xtra_labl:
            varn = lab1[itr]
        else:
            varn = 'iteration '+str(itr)  
        
        if itr==0:         varn= 'Prior' 
        if itr==true_itr:  varn= 'Truth' 
        
        #read netcdf file
        print 'before reading nc file', itr
        xr,yr,hr,ur,vr,zr=read_nc_file(itr,date_sar,surface_vel)
        
        ur=np.ma.masked_array(ur,mask)
        vr=np.ma.masked_array(vr,mask)
        
        uvr=np.sqrt(ur*ur+vr*vr)
        uvrf=uvr[~mask].flatten().data
        num=len(uvrf)
        
        hrf=hr[~mask].flatten().data
    
        #axp.scatter(ur,us,c='r',label='u',alpha=0.5)
        #axp.scatter(vr,vs,c='b',label='v',alpha=0.5)
    
        #if nn==nall:
        #     print '>>>>>>>> replace last one with truth'
        #     uvrf=uvtf
        #     lab1[itr]='Truth forward run'
        
        if True:
            cb_text = 'True depth [m]'
            cc      =  htf[::kj]
        else:
            cb_text = 'Depth [m]'
            cc      =  hrf[::kj]                
        
        
        scr=axp.scatter(uvrf[::kj],uvsf[::kj],c=cc,vmin=lim1,vmax=lim2,s=2,alpha=0.7  \
                        ,edgecolors='none',marker='o')
        axp.legend(loc=2,ncol=1)
                   
        axp.plot([-10,10],[-10,10])
        #axp.set_aspect(1.0)
        axp.set_xlim(0,1.5)
        axp.set_ylim(0,1.5)
        axp.set_xlabel('ROMS Vel [ms$^{-1}$]')
        axp.set_ylabel('SAR  Vel [ms$^{-1}$]')
        axp.set_title(lab1[itr],size='small')
        #axp.set_aspect(1.0)
        
        
    #             fit = np.polyfit(uvrf[::kj],uvsf[::kj],1)
    #             fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    #             fit_txt = str( fit_fn)
    #             fit_txt = fit_txt.replace('\n', '\n y = ')
    #             
    #             xx = np.linspace(-10,10,100) 
    #             plt.plot(xx, fit_fn[xx] , color = 'b', lw=1.25 )
    #             print 'fit info > ',  xlab[ica] , fit_fn          
        
        if nn==1:   
            print 'create colorbar ...'
            ##### colorbar stuff ################################### 
            xsub1=pos_ax[0][0]
            xsub2=pos_ax[1][0]
            ysub1=pos_ax[0][1]            
            ysub2=pos_ax[1][1]
            #print xsub1,xsub2,ysub1,ysub2
            cbax = f.add_axes([xsub1 + 0.2 * (xsub2-xsub1), ysub1 - 0.3 * (ysub2-ysub1), 0.2, 0.02]) 
            cbticks=np.linspace(lim1,lim2 , 3, endpoint=True)
            cb=pl.colorbar(scr,cax=cbax,ticks=cbticks,format='%1.3g',orientation='horizontal') #
            #cb.ax.set_position([xsub1+0.3*(xsub2-xsub1), ysub1+0.08*(ysub2-ysub1), 0.15, 0.05])
            #text(xmin,1.3 * ymax  , color_lab )
            
            #cb.ax.set_position([xsub1-0.5*(xsub2-xsub1), ysub1- 0.205*(ysub2-ysub1), 0.5, 0.03]) # col=3
            #cb.ax.set_position([xsub1-0.0*(xsub2-xsub1), ysub1- 0.0*(ysub2-ysub1), (xsub2-xsub1), 0.02]) # col=4
            
            cb.set_label(cb_text)
            
            [i.set_linewidth(0.05) for i in cb.ax.spines.itervalues()]
            #############################################
        
        if not axp.is_last_row():
           #pl.setp( axp, 'xticklabels', [] )
           axp.set_xlabel('')   
        
        if not axp.is_first_col():
           #pl.setp( axp, 'yticklabels', [] )
           axp.set_ylabel('')   
        
        pl.subplots_adjust(left=left1, bottom=bottom1, right=right1, top=top1,wspace=wspace1, hspace=hspace1)
        #axp.set_position (pos_ax1)
        
    
    print '>> generate fig file'    
    #outfile=pl+'_cross_section_'+dirn+'_'+str(1000+ntime)+'.'+ftype
    #show()
    for ftype in ftypes:
        print ftype
        outfile=outdir+'/scat_'+var+str(100+len(ntr))+'_'+cb_text.replace(' ','_')+'.'+ftype
        if ftype=='png':
           pl.savefig(outfile,dpi=dpi)
        else:
           pl.savefig(outfile)    
    
    #################################################
    #################################################
    #################################################
    print '>>> SCATTER TRUE with ASSIM  >>> '  
    ntr=[0,6]    # itration to be ploted
    nall = icol = len(ntr)
    
    lim1=1
    lim2=7
    #SAR
    print 'SAR >>>>>>>'
    xs,ys,us,vs,date_sar = read_sar_file(time_step=0)    
    uvs  = np.sqrt(us*us+vs*vs)
    mask = us.mask
    
    kj = 1
    
    fwidth,fheight= 6.5,4
    #f=pl.figure(figsize= (fwidth,fheight))   # for 2 x 2 plot
    #====== subplot adjustments ===============
    left1  = 0.1 # the left side of the subplots of the figure
    right1 = 0.9    # the right side of the subplots of the figure
    bottom1= 0.3   # the bottom of the subplots of the figure   (ntr==16   bottom=0.05)
    top1   = 0.9      # the top of the subplots of the figure
    wspace1= 0.05   # the amount of width reserved for blank space between subplots
    hspace1= 0.2   # the amount of height reserved for white space between subplots
    ###############################################<<<<<<
    
    print 'Truth >>>>>>>'
    xt,yt,ht,ut,vt,zt=read_nc_file(true_itr,date_sar,surface_vel)
    uvt=np.sqrt(ut*ut+vt*vt)
    mask = (mask & np.isnan(uvt) & uvt.mask)
    
    uvsf = uvs[~mask].flatten().data
    uvtf = uvt[~mask].flatten().data
    htf  = ht [~mask].flatten().data
    
    nn = -1
    
    irow= int(nall//icol)
    f,axgrid = pl.subplots(nrows=irow, ncols=icol, sharex=True, sharey=True,
                    figsize=(fwidth, fheight),facecolor='w', edgecolor='k',dpi=dpi )
    axgrid = np.array(axgrid)
    axgrid = axgrid.reshape(icol*irow)
    
    for itr in ntr:
        ######## SUBPLOT ARRANGMENTS    #########################
        nn+=1
        print  '  >>> ',nn
        axp=axgrid[nn] 
    
        pos_ax1=axp.get_position ()
        pos_ax=np.array (axp.get_position ())
        #print 'y',pos_ax[0][1], pos_ax[1][1]
        #print 'x',pos_ax[0][0], pos_ax[1][0]
        #print '-----------------------------'
        aheight=pos_ax[1][1] -pos_ax[0][1]
        awidth= pos_ax[1][0] -pos_ax[0][0]
        ##########################################################
        if xtra_labl:
            varn = lab1[itr]
        else:
            varn = 'iteration '+str(itr)  
        
        if itr==0:         varn= 'Prior' 
        if itr==true_itr:  varn= 'Truth' 
        
        #read netcdf file
        print 'before reading nc file', itr
        xr,yr,hr,ur,vr,zr=read_nc_file(itr,date_sar,surface_vel)
        
        ur=np.ma.masked_array(ur,mask)
        vr=np.ma.masked_array(vr,mask)
        
        uvr=np.sqrt(ur*ur+vr*vr)
        uvrf=uvr[~mask].flatten().data
        num=len(uvrf)
        
        hrf=hr[~mask].flatten().data
    
        #axp.scatter(ur,us,c='r',label='u',alpha=0.5)
        #axp.scatter(vr,vs,c='b',label='v',alpha=0.5)
    
        #if nn==nall:
        #     print '>>>>>>>> replace last one with truth'
        #     uvrf=uvtf
        #     lab1[itr]='Truth forward run'
        
        if True:
            cb_text = 'True depth [m]'
            cc      =  htf[::kj]
        else:
            cb_text = 'Depth [m]'
            cc      =  hrf[::kj]                
        
        
        scr=axp.scatter(uvrf[::kj],uvtf[::kj],c=cc,vmin=lim1,vmax=lim2,s=2,alpha=0.7  \
                        ,edgecolors='none',marker='o')
        axp.legend(loc=2,ncol=1)
                   
        axp.plot([-10,10],[-10,10])
        #axp.set_aspect(1.0)
        axp.set_xlim(0,1.5)
        axp.set_ylim(0,1.5)
        axp.set_xlabel('ASIM Vel [ms$^{-1}$]')
        axp.set_ylabel('True Vel [ms$^{-1}$]')
        axp.set_title(lab1[itr],size='small')
        #axp.set_aspect(1.0)
        
        
    #             fit = np.polyfit(uvrf[::kj],uvsf[::kj],1)
    #             fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    #             fit_txt = str( fit_fn)
    #             fit_txt = fit_txt.replace('\n', '\n y = ')
    #             
    #             xx = np.linspace(-10,10,100) 
    #             plt.plot(xx, fit_fn[xx] , color = 'b', lw=1.25 )
    #             print 'fit info > ',  xlab[ica] , fit_fn          
        
        if nn==1:   
            print 'create colorbar ...'
            ##### colorbar stuff ################################### 
            xsub1=pos_ax[0][0]
            xsub2=pos_ax[1][0]
            ysub1=pos_ax[0][1]            
            ysub2=pos_ax[1][1]
            #print xsub1,xsub2,ysub1,ysub2
            cbax = f.add_axes([xsub1 + 0.2 * (xsub2-xsub1), ysub1 - 0.3 * (ysub2-ysub1), 0.2, 0.02]) 
            cbticks=np.linspace(lim1,lim2 , 3, endpoint=True)
            cb=pl.colorbar(scr,cax=cbax,ticks=cbticks,format='%1.3g',orientation='horizontal') #
            #cb.ax.set_position([xsub1+0.3*(xsub2-xsub1), ysub1+0.08*(ysub2-ysub1), 0.15, 0.05])
            #text(xmin,1.3 * ymax  , color_lab )
            
            #cb.ax.set_position([xsub1-0.5*(xsub2-xsub1), ysub1- 0.205*(ysub2-ysub1), 0.5, 0.03]) # col=3
            #cb.ax.set_position([xsub1-0.0*(xsub2-xsub1), ysub1- 0.0*(ysub2-ysub1), (xsub2-xsub1), 0.02]) # col=4
            
            cb.set_label(cb_text)
            
            [i.set_linewidth(0.05) for i in cb.ax.spines.itervalues()]
            #############################################
        
        if not axp.is_last_row():
           #pl.setp( axp, 'xticklabels', [] )
           axp.set_xlabel('')   
        
        if not axp.is_first_col():
           #pl.setp( axp, 'yticklabels', [] )
           axp.set_ylabel('')   
        
        pl.subplots_adjust(left=left1, bottom=bottom1, right=right1, top=top1,wspace=wspace1, hspace=hspace1)
        #axp.set_position (pos_ax1)
        
    
    print '>> generate fig file'    
    #outfile=pl+'_cross_section_'+dirn+'_'+str(1000+ntime)+'.'+ftype
    #show()
    for ftype in ftypes:
        print ftype
        outfile=outdir+'/scat_true_'+var+str(100+len(ntr))+'_'+cb_text.replace(' ','_')+'.'+ftype
        if ftype=='png':
           pl.savefig(outfile,dpi=dpi)
        else:
           pl.savefig(outfile)    
    

        
#########################################################################################
print 'END and backing up script'
args=sys.argv
scr_name=args[0]    
scr_dir1=os.getcwd()
os.system('cp -fr  '+scr_dir1+'/'+scr_name +'    '+outdir+'/scr')

#########################################################################################
    
    # <codecell>
#    
#    imgfile=picdir+'/scat_1008.png'
#    Image(filename=imgfile)

#
# <markdowncell>

# # Bathymetry field for differenet cases:


# <codecell>
#
#
#if True:
#    from IPython.core.display import Image 
#    
#    picdir=cur_dir+'/pic'
#    #!ls -l $picdir
#    #!rm $picdir/*.png
#    imgfile=picdir+'/maps_h1008.png'
#    Image(filename=imgfile)
#    
#    # <markdowncell>
#    
#    # <img src="image.png">
#    
#    # <markdowncell>
#    
#    # # Velocity field based on assimilation of differenet data sources:
#    
#    # <codecell>
#    
#    imgfile=picdir+'/maps_uv1008.png'
#    Image(filename=imgfile)
    
    # <markdowncell>
    
    # # SAR and assimilated forward model velocity corelation:
    
    # <codecell>

