# -*- coding: utf-8 -*-
#from octant.grid import *
#from grid import *
import numpy as np
import pylab as pl
import netCDF4
import netcdftime
import glob
#import datetime
import os,sys

#from   okean import calc
#import scipy.io as sio
from matplotlib.ticker import OldScalarFormatter

sys.path.append('/home/server/pi/homes/moghimi/Desktop/cowast/projects/my_py_libs/plot')
import plot_settings as ps
pl.close('all')


#from octant import plotting as ocpl
#cmap=ocpl.cmap_brightened(ocpl.cmap_discretize(pl.cm.jet,5),0.6)
##cmap=ocpl.cmap_brightened(ocpl.jetWoGn(),0.6)
cmap=pl.cm.jet
#cmap = pl.cm.spectral
#cmap = pl.cm.gist_earth

#cmap=pl.cm.Paired
##############
i1=20
i2=160

j1=30
j2=175

k=1

plot_vel=False
plot_data_points=False

#ftypes=['png','pdf']
ftypes=['png']


# number of iterations
nstart=1
nend=48
nall=nend+1-nstart

#ntr=[0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
#ntr=[0,1,4,6,8,10,12,16]
#nall=len(ntr)


########  read directories  ##################
sys.path.append('../../')
#import base_info
#global base_dir
#global inp_dir
# base directory for this set or iterations
#base_dir=base_info.base_dir
#inp_dir=base_info.inp_dir
#prior=base_info.prior


#########  Read  input files info
itr=1
#base_dir= '/home/shusin3/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/04_wave_limit_test_final/99_wav_cur_true2true/'
#base_dir= '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/800_tru_tru_new_cov'
#base_dir= '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/801_tru_tru_simple_alpha01_hcrit0.5'

base_dir= '/home/shusin4/users/moghimi/assimilation/assim_local/wave_current/run5_syn_paper_cases/05_wave_limit_test_final_new_chh/801_tru_tru_simple_alpha01_hcrit0.5'

base=base_dir+'/run_'+str(1000+itr)+'/04_mem_adj/m*'
flist=glob.glob(base)
flist.sort()

################ figure setting stuff >>>>>
# Publishable quality image
if nall < 2:
   ps.set_mode("publish")
   fig_size1=ps.get_figsize(250)
else:
   ps.set_mode("medium")
   fig_size1=ps.get_figsize(500)

#subplot params
icount=1
rrow=1
# figs param
dpi=400
fwidth,fheight= fig_size1
my_fs=3
fontsize=6

icol=8

fheight=fheight * nall / 9 * 1.1 * 1.15
fwidth=fwidth   * icol / 6 * 1.5 * 5

f=pl.figure(figsize= (fwidth,fheight))   # for 2 x 2 plot
linewidth=0.2
#====== subplot adjustments ===============
left1  = 0.05 # the left side of the subplots of the figure
right1 = 0.95    # the right side of the subplots of the figure
bottom1= 0.1   # the bottom of the subplots of the figure   (ntr==16   bottom=0.05)
top1   = 0.95      # the top of the subplots of the figure
wspace1= 0.01   # the amount of width reserved for blank space between subplots
hspace1= 0.01   # the amount of height reserved for white space between subplots
#####################################################
ncf   = base_dir+'/inp/const/grd_local_20120510.nc'
nc    = netCDF4.Dataset(ncf)
ncvar = nc.variables
mt    = ncvar['mask_rho'] [:]
xt    = ncvar['x_rho'][:]
yt    = ncvar['y_rho'][:]
mtp   = mt==1
mt = np.ma.masked_array(mt,mtp)
#####################################################


nn=0
ntr=pl.arange(0,nall)
for itr in ntr:
    ######## SUBPLOT ARRANGMENTS    #########################
    nn+=1
    irow= int(nall/icol)
       
    print  irow,icol,nn
    axp=f.add_subplot(irow,icol,nn) # (#rows,#columns,subplot#)
    #pl.subplots_adjust(left=left1, bottom=bottom1, right=right1, top=top1,
    #     wspace=wspace1, hspace=hspace1)
    #axp.set_aspect(1.0)
    
    pos_ax1=axp.get_position ()
    pos_ax=np.array (axp.get_position ())
    #print 'y',pos_ax[0][1], pos_ax[1][1]
    #print 'x',pos_ax[0][0], pos_ax[1][0]
    #print '-----------------------------'
    aheight=pos_ax[1][1] -pos_ax[0][1]
    awidth= pos_ax[1][0] -pos_ax[0][0]
    ##########################################################
    #Read in his_file
    #base=base_dir+'/run_'+str(1000+itr)+'/08_forward/'
    #ncf=base+'nri_his.nc'
    ncf=flist[itr]+'/nri_his1000.nc'
    print itr, ncf
    #if copy1: os.system('cp '+ncf+'   nri_his_'+str(1000+itr)+'.nc')
    nc=netCDF4.Dataset(ncf)
    ncvar=nc.variables
#    timename_sim='ocean_time'
#    utime_sim=netcdftime.utime(ncvar[timename_sim].units)
#    sdates_sim=utime_sim.num2date(ncvar[timename_sim][:]) 
#    nt=len(sdates_sim)
#    n=-2
#    
#    ubar=ncvar['ubar']
#    vbar=ncvar['vbar']
    z      = np.squeeze(ncvar['zeta'][:])
    h      = ncvar['h']    [:]
    maskr  = ncvar['mask'] [:]
    xr     = ncvar['x_rho'][:]
    yr     = ncvar['y_rho'][:]
    
#    u=pl.zeros_like(z)
#    v=pl.zeros_like(z)
#    
#    #interpolate to rho points (aprox.)
#    u[:,1:-1]=(ubar[n,:,:-1]+ubar[n,:,1:])/2.
#    u[:,0]   = ubar[n,:,0]
#    u[:,-1]  = ubar[n,:,-1]
#    
#    v[1:-1,:]=(vbar[n,:-1,:]+vbar[n,1:,:])/2.
#    v[0   ,:]=vbar[n,0  ,:]
#    v[-1  ,:]=vbar[n,-1 ,:]   
#    
    m_cut = maskr [j1:j2:k,i1:i2:k]
    x_cut = xr    [j1:j2:k,i1:i2:k]
    y_cut = yr    [j1:j2:k,i1:i2:k]
    h_cut = h     [j1:j2:k,i1:i2:k]
#    u_cut=u     [j1:j2:k,i1:i2:k]
#    v_cut=v     [j1:j2:k,i1:i2:k]
    z_cut =z     [j1:j2:k,i1:i2:k]
#    
#    x_cut=pl.ma.masked_where(m_cut==0,x_cut)
#    y_cut=pl.ma.masked_where(m_cut==0,y_cut)
    h_cut = pl.ma.masked_where(m_cut==0,h_cut)
#    u_cut=pl.ma.masked_where(m_cut==0,u_cut)
#    v_cut=pl.ma.masked_where(m_cut==0,v_cut)
#    z_cut=pl.ma.masked_where(m_cut==0,z_cut)
#    
#    rmaskp=x_cut.mask
#    uv=pl.sqrt(u_cut**2+v_cut**2)
   
    varn= 'iteration '+str(itr)  
    
    if itr==0:  varn= 'Prior' 
    if itr==16: varn= 'Solution' 
    
    varu='[m]'
    val= h_cut + z_cut

    xmin= x_cut.min()   
    xmax= x_cut.max() 
    ymin= y_cut.min()   
    ymax= y_cut.max()

    lim1= 0.25
    lim2= 12
    pcolor_plot=False
    if pcolor_plot:
        pl.pcolor (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),cmap=cmap)
        pl.clim(lim1,lim2)    # evry thing more than 2 will be dark blue
    else:
        levels=pl.linspace(lim1,lim2,25)
        pl.contourf (np.squeeze(x_cut),np.squeeze(y_cut),np.squeeze(val),levels,cmap=cmap)
    #color_lab=varn + varu
    color_lab='depth'+varu
    
    if plot_data_points and itr>1:
        kk=2
        xdata=x_cut  [::kk,::kk]
        ydata=y_cut  [::kk,::kk]
        pl.plot(xdata.flatten(),ydata.flatten(),'k+',markersize=1, alpha=0.4,markeredgecolor='None')#,markerfacecolor='None')
    
    #matplotlib.rcParams['axes.formatter.limits'] = (-1,1)    #to cjange axes to show 1e-1
    #descrip=ps.ordinal[itr]+' '+varn
    #pl.text(xmin+2.5,ymax-2 , ps.ordinal[itr] ,size=12.0)      
    #pl.text(xmin+50,ymax-100 , descrip ,size=11.0)      
    
    if nn==nall- icol // 2:   
       print 'create colorbar ...'
       ##### colorbar stuff ################################### 
       xsub1=pos_ax[0][0]
       xsub2=pos_ax[1][0]
       ysub1=pos_ax[0][1]            
       ysub2=pos_ax[1][1]
       #print xsub1,xsub2,ysub1,ysub2
       cbticks=np.linspace(lim1, lim2, 6, endpoint=True)
       cb=pl.colorbar(ticks=cbticks,shrink=0.6,format=OldScalarFormatter() ,orientation='horizontal') #
       #cb.ax.set_position([xsub1+0.3*(xsub2-xsub1), ysub1+0.08*(ysub2-ysub1), 0.15, 0.05])
       #text(xmin,1.3 * ymax  , color_lab )
    
       #cb.ax.set_position([xsub1-2.3*(xsub2-xsub1), ysub1- 0.22*(ysub2-ysub1), 0.5, 0.03]) # col=4
       cb.ax.set_position([xsub1-0.5*(xsub2-xsub1), ysub1- 0.25*(ysub2-ysub1), 0.5, 0.02]) # col=4
       
       cb.set_label(color_lab)
    
       [i.set_linewidth(0.05) for i in cb.ax.spines.itervalues()]
       #############################################
    pl.pcolor(xt,yt,mt,cmap=pl.cm.bone)
    
    # contour stuff
    #matplotlib.rcParams['contour.negative_linestyle'] = 'solid' #all lines be solid
    levels=pl.linspace(lim1,lim2,7)
    #Cabc = pl.contour (np.squeeze(x_cut),np.squeeze(y_cut),val,levels,colors='k',linewidths=linewidth)
    Cabc = pl.contour (np.squeeze(x_cut),np.squeeze(y_cut),h_cut,levels,colors='k',linewidths=linewidth)
    
    #title1= dirn+'   '+str(utime.num2date(ncv['time'][ntime]))
    #title(title1,fontsize=11)
    ##################################################          
    # the linewidth of the rectangular axis frame
    ax=pl.gca()
    fr_linewidth=0.3
    [i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
    
    pl.ylim(ymin,ymax)
    pl.xlim(xmin,xmax)
    
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

    
#pl.title(title1)
outdir=base_dir+'/post/maps_mems/'
command='mkdir -p '+ outdir
os.system(command)
for ftype in ftypes:
    outfile=outdir+'/maps'+str(1000+len(ntr))+'.'+ftype
    if ftype=='png':
       pl.savefig(outfile,dpi=300)
    else:
       pl.savefig(outfile)    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#pl.close('all')
#
#
#
#
#
#
#np=len(uv.flatten())
#cl = [cmap(float(i)/(np)) for i in xrange(np)]
#
#
#fig1 = pylab.figure(1)
#fig2 = pylab.figure(2)
#ax1 = fig1.add_subplot(111)
#ax2 = fig2.add_subplot(111)
#ax1.scatter(ha_cut[~rmask],h_cut[~rmask],c=cl,alpha=0.5)
#im=ax2.scatter(ha_cut[~rmask],h_cut[~rmask],c=cl,alpha=0.5)
#fig1.colorbar(im , ax=ax1)
#fig1.show()
#
#
#import matplotlib
#import matplotlib.pyplot as plt
#cm = matplotlib.cm.get_cmap('RdYlBu')
#zc = range(np)
#ref=h_cut[~rmask].flatten()
#asm=ha_cut[~rmask].flatten()
#sc = plt.scatter(ref,asm, c=zc, vmin=0, vmax=15, s=5, cmap=cm)
#plt.colorbar(sc)
#plt.show()
#
#
#
#
#
#fig1 = pylab.figure(1)
#fig2 = pylab.figure(2)
#ax1 = fig1.add_subplot(111)
#ax2 = fig2.add_subplot(111)
#ax1.scatter(range(10), range(10), c=range(10), alpha=0.2)
#im = ax2.scatter(range(10), range(10), c=range(10), alpha=1.0)
#fig1.colorbar(im, ax=ax1)
#fig1.show()
#
#
#
#import matplotlib
#import matplotlib.pyplot as plt
#cm = matplotlib.cm.get_cmap('RdYlBu')
#xy = range(20)
#z = xy
#sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
#plt.colorbar(sc)
#plt.show()
#
#
#
#
#
#
#
#
#
#
#
#
#

