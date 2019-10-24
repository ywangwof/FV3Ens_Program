import matplotlib
matplotlib.use("Agg")
import skewx
import numpy as np
import math
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.watch_type as watch_type
from sharppy.sharptab.constants import *
import sharppy.viz.barbs as barbs
import matplotlib.cbook as cbook
import sharppy.databases.inset_data as inset_data
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime, timedelta
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_colormaps
from news_e_colormaps import *
import sharppy.sharptab as sharptab

#plt.rcParams['font.family'] = 'Liberation Sans'
def plot_wof(prof, members, figname, xlat, xlon, idateobj, vdateobj, **kwargs):
#    '''
#    Plots SHARPpy SPC window as .png
#
#    Parameters
#    ----------
#    prof : a Profile Object from sharppy.sharptab.profile
#    
#    kwargs
#    ------
#    parcel_type: Parcel choice for plotting. 'sfc','ml','mu','fcst' Default is 'ml'
#    filename: Filename as a string. Default is 'sounding.png'
#    logo: Logo for upper-left portion of the skew-t. Default is 'None' and does not plot a logo.
#    logo_dxdy: Size of logo. Actual dimensions are dT and dp as it is plotted on the skewT. Default is (20,13) This worked for a 489x132 pixel image.
#    '''
    #kwargs
   parcel_type = kwargs.get('parcel_type', 'ml')
   xpts = kwargs.get('x_pts')
   ypts = kwargs.get('y_pts')

#Figure User Input

   p_grid_labels = ['1000','','','850','','','700','','','','500','','','','300','','200','','100'] #labels for the pressure ticks. Standard.
   p_grid = [1000,850,700,500,300,200,100] #where horizontal lines go across the skew-T

   my_dpi = 55 #dots per inch for the plot. This is a pretty hi-res image.

   pmax = 1050 #lowest pressure on the skew-T
   pmin = 100 #highest pressure on the skew-T
   dp = -10 #pressure spacing for creating skew-T background lines 

   presvals = np.arange(int(pmax), int(pmin)+dp, dp) #pressure values used for creating skew-T background lines

# Colors for wind speed bars and hodograph
   hgt_list_bar = [0,1000, 3000,6000,9000,20000]
   hgt_list_hodo = [0,1000, 3000,6000,9000,10000]

   hodo_color = [cb_colors.orange6, cb_colors.green6, cb_colors.blue6, cb_colors.purple6,cb_colors.red6]
   hodo_label = ['0-1km', '1-3km','3-6km','6-9km','9-10km']

#convoluted mess to get the title to be aligned how I wanted. This should be changed for others...
   spaces=10#22
#   title_text_1 = '' #site + ' ' + dt.strftime('%Y/%m/%d %H:%M UTC ' + data_type)
#   title_text_3 = 'Sounding Powered by SHARPpy'
   sharptext = 'Sounding Powered by SHARPpy'
#   title_text_2 = title_text_3 = ''
   title_text_3 = 'WoFS Sounding {}N, {}W'.format(xlat,xlon) + (' '*spaces) + 'Init: {}     Valid: {}'.format(idateobj.strftime('%Y-%m-%d %H%M UTC'),vdateobj.strftime('%Y-%m-%d %H%M UTC'))
   #xlat, xlon, initdate, validdate 
   title_text = title_text_3 #title_text_1 + (' '*spaces) +title_text_2 + (' '*spaces) + title_text_3

#Figures out where at which height the sounding reached in the list of h_new
#Then interpolates pressure to height levels
   h_new = [0,1000,3000,6000,9000,12000,15000]
    
   for i in range(len(h_new)):
      if np.max(prof.hght)>h_new[i]:
         index=i
   h_new_labels = ['0 km','1 km','3 km','6 km','9 km','12 km','15 km']
   h_new_labels = h_new_labels[0:index+1]
    #p_interped_func = interpolate.interp1d(prof.hght, prof.pres)
    
   p_interped = sharptab.interp.pres(prof, sharptab.interp.to_msl(prof, h_new[0:index+1]))

#Thin out the winds for better plotting (significant level data points seem to bunch together to closely
   minimum_separation = 250. #minimum spacing between wind barbs (meters)
   h_barb = np.array(prof.hght).tolist()
   p_barb = np.array(prof.pres).tolist()
   spd_barb = np.array(prof.wspd).tolist()
   direc_barb = np.array(prof.wdir).tolist()
   
#adds units to our newly created pressure, speed, and direction arrays for wind barb plotting
   #p_barb = p_barb * units.mbar
   #spd_barb = spd_barb * units.knot
   #direc_barb = direc_barb * units.deg

# Convert wind speed and direction to components
   #u, v = get_wind_components(prof.wspd * units.knot, prof.wdir * units.deg)
   u, v = utils.vec2comp(prof.wdir, prof.wspd)
   u_barb, v_barb = utils.vec2comp(prof.wdir, prof.wspd)

#SELECT PARCEL AND GET PARCEL DATA FROM SPC_UTILS
   sfcpcl = prof.sfcpcl#params.parcelx( prof, flag=1 )
   fcstpcl = prof.fcstpcl#params.parcelx( prof, flag=2)
   mupcl = prof.mupcl#params.parcelx( prof, flag=3 )
   mlpcl = prof.mlpcl#params.parcelx( prof, flag=4 )
   if parcel_type=='sfc':
      pcl = sfcpcl
      pcl_box_level = -0.065
      pcl_type = 1
   elif parcel_type=='fcst':
      pcl = fcstpcl
      pcl_box_level =  -0.0875
      pcl_type = 2
   elif parcel_type=='mu':
      pcl = mupcl
      pcl_box_level = -0.1325
      pcl_type = 4
   elif parcel_type=='ml':
      pcl = mlpcl
      pcl_box_level = -0.11
      pcl_type = 3
   else:
      print( "ERROR! Select 'sfc', 'fcst', 'mu', or 'ml' for parcel type. (plot_spc(prof,parcel_type)")
      print( "Defaulting to surface parcel...")
      pcl = sfcpcl
      pcl_box_level = -0.065


#PLOTTING *************************************************************************************************************

#Create full figure
   fig = plt.figure(figsize=(1180/my_dpi, 800/my_dpi), dpi=my_dpi, frameon=False)

#SKEW T ***************************************************
   ax = fig.add_subplot(111, projection='skewx', facecolor="w") #skewed x-axis

# plot dashed temperature lines
   ax.xaxis.grid(color='k', linestyle='--', dashes=(3,3), alpha=0.5,zorder=0)

# plot the moist-adiabats
   for temp in np.arange(-10,45,5):
      tw = []
      for pres in presvals:
         tw.append(thermo.wetlift(1050., temp, pres))
      ax.semilogy(tw, presvals, color=cb_colors.purple6,linestyle='--', dashes=(5,2), alpha=.3)#cb_colors.purple6

# plot the dry adiabats
   for t in np.arange(-50,80,20):
      thetas = ((t + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK
      ax.semilogy(thetas, presvals, 'k', alpha=.3)

#plot mixing ratio lines
   mixing_ratio_list = range(6,36,4)
   for mr in mixing_ratio_list:
      plt.plot((thermo.temp_at_mixrat(mr, 1050)-273, thermo.temp_at_mixrat(mr, 600)-273), (1050,600), 'g-',lw=1.0, zorder=3, alpha=0.6)
      ax.annotate(str(mr), xy=((thermo.temp_at_mixrat(mr, 600)-273), (600-3)), xytext=((thermo.temp_at_mixrat(mr, 600)-273), (600-3)), ha='center', color='g', family='sans-serif', weight='bold', zorder=3, fontsize=10, alpha=0.6)
    
#plot horizontal lines at standard pressure levels
   for p_loc in p_grid:
      ax.axhline(y=p_loc,ls='-',c='k',alpha=0.5, linewidth=1.5, zorder=3)


# PLOT THE DATA ON THE SKEW-T
    
# Plot the data using normal plotting functions, in this case using log scaling in Y, as dicatated by the typical meteorological plot

   ax.semilogy(prof.wetbulb, prof.pres, c="c", linestyle='-', lw=1, alpha=1.0, zorder=3) # Plot the wetbulb profile
   ax.annotate(str(int(round(thermo.ctof(prof.wetbulb[prof.sfc])))), xy=(prof.wetbulb[prof.sfc], prof.pres[prof.sfc]+30), xytext=(prof.wetbulb[prof.sfc], prof.pres[prof.sfc]+30),ha='left', color="c", family='sans-serif', weight='normal', zorder=7, fontsize=12, alpha=1.0) # annotate the sfc wetbulb in F
   ax.semilogy(prof.dwpc, prof.pres, c=cb_colors.blue6, linestyle='-', lw=3, zorder=3) # plot the dewpoint profile
   ax.annotate(str(int(round(thermo.ctof(prof.dwpc[prof.sfc])))), xy=(prof.dwpc[prof.sfc], prof.pres[prof.sfc]+30), xytext=(prof.dwpc[prof.sfc], prof.pres[prof.sfc]+30), ha='left', color=cb_colors.blue6, family='sans-serif', weight='bold', zorder=7, fontsize=12, alpha=1.0) # annotate the sfc dewpoint in F
   ax.semilogy(prof.tmpc, prof.pres, c=cb_colors.orange6, linestyle='-', lw=3, zorder=3) # Plot the temperature profile
   ax.semilogy(prof.vtmp, prof.pres, c=cb_colors.orange6, linestyle='--', lw=3, zorder=3) # Plot the temperature profile
   ax.annotate(str(int(round(thermo.ctof(prof.tmpc[prof.sfc])))), xy=(prof.tmpc[prof.sfc], prof.pres[prof.sfc]+30), xytext=(prof.tmpc[prof.sfc], prof.pres[prof.sfc]+30), ha='left', color=cb_colors.orange6, family='sans-serif', weight='bold', zorder=7, fontsize=12, alpha=1.0) # annotate the sfc temp in F
   ax.semilogy(pcl.ttrace, pcl.ptrace, c=cb_colors.gray6, linestyle='--', dashes=(3,3), lw=1.5, alpha=1.0, zorder=3) # plot the parcel trace 

   #member_cape = []    
   if members is not None:
   #   print( "Plotting members...")
      for m_idx in range(len(members['tmpc'])):
         ax.semilogy(members['dwpc'][m_idx], members['pres'][m_idx], c=cb_colors.blue6, linestyle='-', lw=1., zorder=1, alpha=.6) # plot the dewpoint profile
         ax.semilogy(members['tmpc'][m_idx], members['pres'][m_idx], c=cb_colors.orange6, linestyle='-', lw=1., zorder=1, alpha=.6) # Plot the temperature profile
# set label and tick marks for pressure and temperature
   ax.xaxis.set_major_locator(plt.MultipleLocator(10))
   ax.set_xticks(np.arange(-100,60,10))
   ax.set_xticklabels([str(i) for i in np.arange(-100,60,10)],color=cb_colors.gray7,fontsize=12)
   ax.set_xlim(-50,50)
   ax.yaxis.set_major_formatter(plt.ScalarFormatter())
   ax.set_yticks(np.linspace(1000,100,19))
   ax.set_yticklabels(p_grid_labels,color=cb_colors.gray7,fontsize=12)
   ax.set_ylim(1050,100)

#plot the title text
   plt.text(0.05, 0.97, title_text, fontsize=15, color=cb_colors.gray8, weight='bold', ha='left', transform=fig.transFigure)
   #x_hodo.annotate(sharptext, xy=(0.95, 0.95), xytext=(0.95, 0.95),xycoords='axes fraction',textcoords='axes fraction',ha='center', va='bottom', color=cb_colors.gray7, family='sans-serif', weight='bold', zorder=3,fontsize=14)
   plt.text(0.8, 0.97, sharptext, fontsize=15, color=cb_colors.gray8, weight='bold', ha='left', transform=fig.transFigure)
#adjust the skew-T plot to make room for the rest of the figures. This was important to make everything line up.
   plt.subplots_adjust(left=0.05, right=0.55, top=0.96, bottom=0.15)

#Plot the height labels on the left axis
   ax2 = ax.twinx() #makes a twin of the skew-T subplot that's not skewed at 45 degrees
   plt.yscale('log', nonposy='clip')
   plt.yticks(p_interped,h_new_labels,color=cb_colors.green4, ha='left')
   ax2.yaxis.tick_left()
   ax2.tick_params(direction='in', pad=-15, axis='both', which='major', colors=cb_colors.green4,length=10,width=1.5)
   ax2.set_yticklabels(h_new_labels, fontsize=12, weight='bold',color=cb_colors.green4)

   x = np.random.uniform(0.0,10.0,15)
   y = np.random.uniform(0.0,10.0,15)

# Plot LCL and LFC levels
   plt.plot((38, 42), (pcl.lfcpres,pcl.lfcpres), c="darkgreen",lw=2.0, zorder=3)
   ax2.annotate('LFC', xy=(40, pcl.lfcpres), xytext=(40, pcl.lfcpres),ha='center', va='bottom', color=cb_colors.green5, family='sans-serif', weight='bold', zorder=3,fontsize=12)
   plt.plot((38, 42), (pcl.lclpres,pcl.lclpres), c="r", linestyle='-',lw=2.0, zorder=3)
   ax2.annotate('LCL', xy=(40, pcl.lclpres+5.), xytext=(40, pcl.lclpres+5.),ha='center', va='top', color=cb_colors.red5, family='sans-serif', weight='bold', zorder=3,fontsize=12)
   plt.plot((38, 42), (pcl.elpres,pcl.elpres), c="m", linestyle='-',lw=2.0, zorder=3)
   ax2.annotate('EL', xy=(40, pcl.elpres), xytext=(40, pcl.elpres),ha='center', va='bottom', color=cb_colors.purple5, family='sans-serif', weight='bold', zorder=3,fontsize=12)
    
# Plot Eff Inflow Layer
   eff_inflow = (prof.ebottom, prof.etop)
   eff_inflow_top = sharptab.interp.to_agl(prof, sharptab.interp.hght(prof, eff_inflow[1]))
   eff_inflow_bottom = sharptab.interp.to_agl(prof, sharptab.interp.hght(prof, eff_inflow[0]))
   bunkers = prof.srwind
   effective_srh = prof.right_esrh
   plt.plot((-25, -20), (eff_inflow[0],eff_inflow[0]), c=cb_colors.red5, linestyle='-',lw=1.75, zorder=3)
   plt.plot((-25, -20), (eff_inflow[1],eff_inflow[1]), c=cb_colors.red5, linestyle='-',lw=1.75, zorder=3)
   plt.plot((-22.5, -22.5), (eff_inflow[0],eff_inflow[1]), c=cb_colors.red5, linestyle='-',lw=1.75, zorder=3)
   try:
      plt.annotate(str(int(eff_inflow_bottom))+'m  ', xy=(-25, eff_inflow[0]), xytext=(-25, eff_inflow[0]), ha='right', va='center', color=cb_colors.red5, zorder=3, fontsize=12, weight='bold')
      plt.annotate(str(int(eff_inflow_top))+'m  ', xy=(-25, eff_inflow[1]), xytext=(-25, eff_inflow[1]), ha='right', va='center', color=cb_colors.red5, zorder=3, fontsize=12, weight='bold')
      plt.annotate(str(int(effective_srh[0]))+' m$^2$/s$^2$', xy=(-22.5, eff_inflow[1]-10), xytext=(-22.5, eff_inflow[1]-10), ha='center', va='bottom', color=cb_colors.red5, zorder=3, fontsize=12, weight='bold')
   except:
      print ("NO EFF INFLOW")


# PLOT WINDBARBS
   p_barb = np.asarray(p_barb)
   pidx = idx = np.where(np.asarray(p_barb) >= 100.)[0]
   wind_barbs = ax2.barbs(55*np.ones(len(p_barb[idx])), p_barb[idx], u_barb[idx], v_barb[idx], barbcolor=cb_colors.gray7, flagcolor='None', zorder=10, lw=1.0, length=7)
   wind_barbs.set_clip_on(False)

   ax2.invert_yaxis()
   ax2.set_xlim(-50,50)
   ax2.set_ylim(1050,100)

   spd_barb = np.asarray(spd_barb)

# Create hodograph ********************************************************************************************
   ax_hod = fig.add_axes([0.60, 0.45, 0.38, 0.475], frameon=False) #, facecolor='k')

# Set the characteristics of the tick marks
   for tick in ax_hod.xaxis.get_major_ticks():
      tick.label.set_fontsize(12)
      tick.label.set_color(cb_colors.gray7)
      tick.label.set_weight('bold')
   for tick in ax_hod.yaxis.get_major_ticks():
      tick.label.set_fontsize(12)
      tick.label.set_color(cb_colors.gray7)
      tick.label.set_weight('bold')
   for i in range(10,90,10):

# Draw the range rings around the hodograph.
      circle = plt.Circle((0,0),i,linestyle='--',color='k',alpha=.3, fill=False)
      ax_hod.add_artist(circle)

# Set the tick parameters to displace the tick labels from the hodograph axes
   ax_hod.tick_params(axis='x', which='major', labelsize=10,color=cb_colors.gray7, pad=-235, length=0)
   ax_hod.tick_params(axis='y', which='major', labelsize=10,color=cb_colors.gray7, pad=-315, length=0)

# Plot the hodograph axes
   ax_hod.axvline(0, color=cb_colors.gray7, linestyle='-',linewidth=2.)
   ax_hod.axhline(0, color=cb_colors.gray7, linestyle='-',linewidth=2.)
   ax_hod.set_yticks(range(-60,65,10))
   ax_hod.set_xticks(range(-70,75,10))
   hod_yticklabels = [str(abs(i)) for i in range(-60,65,10)]
   #hod_yticklabels[len(hod_yticklabels)/2] = ''
   hod_xticklabels = [str(abs(i)) for i in range(-70,75,10)]
   #hod_xticklabels[len(hod_xticklabels)/2] = ''
   ax_hod.set_yticklabels(hod_yticklabels, fontsize=12, weight='bold', color=cb_colors.gray7)
   ax_hod.set_xticklabels(hod_xticklabels, fontsize=12, weight='bold', color=cb_colors.gray7)
    
# Plot the hodograph using the color scheme for different layers (0-3, 3-6, etc.)
   bounds = [0, 1000, 3000, 6000, 9000, 12000]
   for i in range(1, len(bounds),1):
      subset_idxs = np.where( (prof.hght <= sharptab.interp.to_msl(prof, bounds[i])) & (prof.hght >= sharptab.interp.to_msl(prof, bounds[i-1])))
      subset_hghts = np.ma.concatenate(([sharptab.interp.to_msl(prof, bounds[i-1])], prof.hght[subset_idxs], [sharptab.interp.to_msl(prof, bounds[i])]))
      u, v = sharptab.interp.components(prof, sharptab.interp.pres(prof, subset_hghts))
      ax_hod.plot(u, v, c=hodo_color[i-1], linewidth=2.5,label=hodo_label[i-1],zorder=3)
    
   if members is not None:
      for mprof in members['member_profs']:
         for i in range(1, len(bounds),1):
            subset_idxs = np.where( (mprof.hght <= sharptab.interp.to_msl(mprof, bounds[i])) & (mprof.hght >= sharptab.interp.to_msl(mprof, bounds[i-1])))
            subset_hghts = np.ma.concatenate(([sharptab.interp.to_msl(mprof, bounds[i-1])], mprof.hght[subset_idxs], [sharptab.interp.to_msl(mprof, bounds[i])]))
            u, v = sharptab.interp.components(mprof, sharptab.interp.pres(mprof, subset_hghts))
            ax_hod.plot(u, v, c=hodo_color[i-1], linewidth=1.25, alpha=0.6,label=hodo_label[i-1], zorder=1)

# Get the Bunkers storm motions and convert them to strings to plot
   bunkers = srwind = prof.srwind
   bunkers_rt = utils.comp2vec(bunkers[0],bunkers[1])
   bunkers_lf = utils.comp2vec(bunkers[2],bunkers[3])
   bunkers_rt_str = str(int(np.ma.around(bunkers_rt[0],0)))+"/"+str(int(np.ma.around(bunkers_rt[1],0)))
   bunkers_lf_str = str(int(np.ma.around(bunkers_lf[0],0)))+"/"+str(int(np.ma.around(bunkers_lf[1],0)))
    
# Plot the effective inflow layer on the hodograph
   effubot, effvbot = sharptab.interp.components(prof, eff_inflow[0])
   effutop, effvtop = sharptab.interp.components(prof, eff_inflow[1])
   ax_hod.plot([effubot, srwind[0]], [effvbot, srwind[1]], c='c', linewidth=1.5)
   ax_hod.plot([effutop, srwind[0]], [effvtop, srwind[1]], c='c', linewidth=1.5)
    
# Annotate where the Bunkers storm motion vectors are on the hodograph
   ax_hod.plot(srwind[0], srwind[1], marker='o', fillstyle='none', markeredgecolor=cb_colors.blue8, markeredgewidth=1.5, markersize=11)
   ax_hod.annotate(bunkers_rt_str+' RM', (srwind[0]+1.5, srwind[1]-1.5), fontsize=12, va="top", ha="left", color='k', weight='bold',zorder=10)
   ax_hod.plot(srwind[2], srwind[3], marker='o', fillstyle='none', markeredgecolor=cb_colors.blue8, markeredgewidth=1.5, markersize=11)
   ax_hod.annotate(bunkers_lf_str+' LM', (srwind[2]+1.5, srwind[3]-1.5), fontsize=12, va="top", ha="left", color='k', weight='bold',zorder=10)

# Annotate where the Corfidi MBE vectors are on the hodograph
   corfidi = prof.upshear_downshear
   corfidi_up = utils.comp2vec(corfidi[0],corfidi[1])
   corfidi_dn = utils.comp2vec(corfidi[2],corfidi[3])
   c = 'k'#'#0A74C6'
   ax_hod.plot(corfidi[0], corfidi[1], marker='o', fillstyle='none', markeredgecolor=c, markeredgewidth=1.5, markersize=9)
   ax_hod.annotate(str(int(corfidi_up[0])) + '/' + str(int(corfidi_up[1])) +' UP', (corfidi[0]+1.5, corfidi[1]-1.5), fontsize=12, va="top", ha="left", color=cb_colors.purple8, weight='bold',zorder=10)
   ax_hod.plot(corfidi[2], corfidi[3], marker='o', fillstyle='none', markeredgecolor=c, markeredgewidth=1.5, markersize=9)
   ax_hod.annotate(str(int(corfidi_dn[0])) + '/' + str(int(corfidi_dn[1])) +' DN', (corfidi[2]+1.5, corfidi[3]-1.5), fontsize=12, va="top", ha="left", color=cb_colors.purple8, weight='bold',zorder=10)
 
# Get the cloud-layer mean wind
   mean_cloudlayer = winds.mean_wind(prof, pbot= pcl.lclpres, ptop=pcl.elpres)
   mean_cloudlayer_comp = utils.comp2vec(mean_cloudlayer[0],mean_cloudlayer[1])
   try:
      mean_cloudlayer_str = str(int(np.ma.around(mean_cloudlayer_comp[0],0)))+"/"+str(int(np.ma.around(mean_cloudlayer_comp[1],0)))
   except:
      mean_cloudlayer_str = 'M/M'

# Write the critical angle to the hodograph.
   if eff_inflow[0] == prof.pres[prof.sfc]:
      ax_hod.annotate('Critical Angle = '+str(int(prof.critical_angle)), (-65, -50), fontsize=12, va="bottom", ha="left", color=cb_colors.green6, weight='bold',zorder=10)

   ax_hod.set_xlim(-80,80)
   ax_hod.set_ylim(-70,70)

#BELOW IS STUFF FOR BOXES/BORDERS ******************************************************************************

   ax3 = ax2.twinx()
   ax3.axes.get_xaxis().set_visible(False)
   ax3.axes.get_yaxis().set_visible(False)
   ax3.set_yticks([])
   ax3.set_yticklabels([])

#Big Thick Box around Skew-T
   #box = ax3.add_patch(patches.Rectangle((-50, 0), 110.0, 1.0,fill=False,linewidth=2,edgecolor="w",zorder=3))
   #box.set_clip_on(False)

#box around hodograph
   #box = ax_hod.add_patch(patches.Rectangle((-80., -70.), 160., 140.,fill=False,linewidth=2,edgecolor="w",zorder=4))
   #box.set_clip_on(False)

   inset_color = cb_colors.gray7

#THICK TEXT BOX around Thermodynamics Text
   box = ax3.add_patch(patches.Rectangle((-58.0, -0.035), 54, -0.12,fill=False,linewidth=2,edgecolor=inset_color,zorder=4))
   box.set_clip_on(False)

#THICK TEXT BOX around Kinematics Text
   box = ax3.add_patch(patches.Rectangle((-3.0, -0.035), 60, -0.12,fill=False,linewidth=2,edgecolor=inset_color,zorder=4))
   box.set_clip_on(False)


#THICK TEXT BOX Around Dynamics Text
#   box = ax3.add_patch(patches.Rectangle((9.0, -0.04), 60, -0.38,fill=False,linewidth=2,edgecolor=inset_color,zorder=4))
#   box.set_clip_on(False)

#THICK TEXT BOX Around SARS Text
#   box = ax3.add_patch(patches.Rectangle((69.0, -0.04), 55, -0.38,fill=False,linewidth=2,edgecolor=inset_color,zorder=4))
#   box.set_clip_on(False)

#Thermodynamics
   #box = ax3.add_patch(patches.Rectangle((-55.0, -0.04), 64.0, -0.12,fill=False,linewidth=1,edgecolor=inset_color,zorder=4))
   #box.set_clip_on(False)
   #box = ax3.add_patch(patches.Rectangle((-55.0, -0.04), 64.0, -0.025,fill=False,linewidth=1,edgecolor=inset_color,zorder=4))
   #box.set_clip_on(False)
    
# Write the parcel properties to the inset.
   #x_list = [0, 0.08, 0.17, 0.25, 0.33, 0.39, 0.47]
   x_list = np.array([0, 0.08, 0.17, 0.25, 0.33, 0.39, 0.47])-0.075
   y_list = [-0.045, -0.07, -0.0925, -0.115, -0.1375]
   A = ["SFC", prof.sfcpcl.bplus, int(prof.sfcpcl.bminus), prof.sfcpcl.lclhght,prof.sfcpcl.li5, prof.sfcpcl.lfchght, prof.sfcpcl.elhght]
#   B = ["FCST", prof.fcstpcl.bplus, int(prof.fcstpcl.bminus), prof.fcstpcl.lclhght, prof.fcstpcl.li5, prof.fcstpcl.lfchght, prof.fcstpcl.elhght]
   C = ["ML", prof.mlpcl.bplus, int(prof.mlpcl.bminus), prof.mlpcl.lclhght, prof.mlpcl.li5, prof.mlpcl.lfchght, prof.mlpcl.elhght]
   D = ["MU", prof.mupcl.bplus, int(prof.mupcl.bminus), prof.mupcl.lclhght, prof.mupcl.li5, prof.mupcl.lfchght, prof.mupcl.elhght]
   #mlcape = C[1]
   #print('mlcape',mlcape)
   data = np.array([["PCL", "CAPE", "CINH", "LCL", "LI", "LFC", "EL"],
                      [ str(int(np.ma.around(A[i],0))) if (type(A[i])== np.float64) else str(A[i]) for i in range(len(A)) ],
                      [ str(int(np.ma.around(C[i],0))) if (type(C[i])== np.float64) else str(C[i]) for i in range(len(C)) ],
                      [ str(int(np.ma.around(D[i],0))) if (type(D[i])== np.float64) else str(D[i]) for i in range(len(D)) ]])

   for i in range(data.shape[0]):
      for j in range(data.shape[1]):
         ax2.annotate(data[i,j], (x_list[j], y_list[i]), xycoords="axes fraction", fontsize=12, va="top", ha="left", color=cb_colors.gray7, weight='bold')

# Draw a box around the selected parcel being shown in the Skew-T
   box = ax3.add_patch(patches.Rectangle((-57.7, pcl_box_level), 53., 0.0225,fill=False,linewidth=1,edgecolor=cb_colors.purple4,zorder=4))
   box.set_clip_on(False)

# Write the lapse rates to the inset.
   #box = ax3.add_patch(patches.Rectangle((-55.0, -0.305), 43.0, -0.115,fill=False,linewidth=1,edgecolor="w",zorder=4))
   #box.set_clip_on(False)
   x_list = [0.15,0.16]
   y_list = np.arange(-0.315, -.40, -0.0225)
   data = np.array([["0-3km AGL LR =",str(np.ma.around(prof.lapserate_3km,1))+" C/km"],
                       ["3-6km AGL LR =",str(np.ma.around(prof.lapserate_3_6km,1))+" C/km"],
                       ["850-500mb LR =",str(np.ma.around(prof.lapserate_850_500,1))+" C/km"],
                       ["700-500mb LR =",str(np.ma.around(prof.lapserate_700_500,1))+" C/km"]])
   for i in range(data.shape[0]):
      for j in range(data.shape[1]):
         if (j % 2 == 0):
            ax2.annotate(data[i,j], (x_list[j], y_list[i]), xycoords="axes fraction",
                        fontsize=12, va="top", ha="right", color='k', weight='bold')
         else:
            ax2.annotate(data[i,j], (x_list[j], y_list[i]), xycoords="axes fraction",
                        fontsize=12, va="top", ha="left", color='k', weight='bold')

#Severe Indices
   #box = ax3.add_patch(patches.Rectangle((-12.0, -0.305), 21.0, -0.115,fill=False,linewidth=1,edgecolor="w",zorder=4))
   #box.set_clip_on(False)
   x_list = [0.52,0.53]
   y_list = np.arange(-0.315, -.40, -0.0225)

# This looks lifted from the Profile class.  Don't need this.
   sfc = prof.pres[prof.sfc]
   p6km = sharptab.interp.pres(prof, sharptab.interp.to_msl(prof, 6000.))
   p8km = sharptab.interp.pres(prof, sharptab.interp.to_msl(prof, 8000.))
#   ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
#   etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
    
# Mean winds
   #mean_1km = winds.mean_wind(prof, pbot=sfc, ptop=p1km)
   mean_1km_comp = prof.mean_1km#utils.comp2vec(mean_1km[0],mean_1km[1])
   #mean_3km = winds.mean_wind(prof, pbot=sfc, ptop=p3km)
   mean_3km_comp = prof.mean_3km #utils.comp2vec(mean_3km[0],mean_3km[1])
   #mean_eff = winds.mean_wind(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
   mean_eff_comp = utils.comp2vec(*prof.mean_eff)#utils.comp2vec(mean_eff[0],mean_eff[1])
    
   if type(eff_inflow[0])!=np.float64:
      mean_eff_comp = ['---','--']
   mean_6km = winds.mean_wind(prof, pbot=sfc, ptop=p6km)
   mean_6km_comp = prof.mean_6km #utils.comp2vec(mean_6km[0],mean_6km[1])
   mean_8km = winds.mean_wind(prof, pbot=sfc, ptop=p8km)
   mean_8km_comp = prof.mean_8km #utils.comp2vec(mean_8km[0],mean_8km[1])
   mean_cloudlayer = winds.mean_wind(prof, pbot= pcl.lclpres, ptop=pcl.elpres)
   mean_cloudlayer_comp = prof.mean_lcl_el#utils.comp2vec(mean_cloudlayer[0],mean_cloudlayer[1])
   mean_ebwd = winds.mean_wind(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
   mean_ebwd_comp = utils.comp2vec(*prof.mean_ebw)#utils.comp2vec(mean_ebwd[0],mean_ebwd[1])
    
   if type(eff_inflow[0])!=np.float64:
      mean_ebwd_comp = ['---','--']
    
   bunkers_rt = utils.comp2vec(bunkers[0],bunkers[1])
   bunkers_lf = utils.comp2vec(bunkers[2],bunkers[3])
   corfidi_up = utils.comp2vec(corfidi[0],corfidi[1])
   corfidi_dn = utils.comp2vec(corfidi[2],corfidi[3])
   srw_1km = prof.srw_1km #utils.comp2vec(*winds.sr_wind(prof, pbot=sfc, ptop=p1km, stu=prof.srwind[0], stv=prof.srwind[1]))
   srw_3km = prof.srw_3km #utils.comp2vec(*winds.sr_wind(prof, pbot=sfc, ptop=p3km, stu=prof.srwind[0], stv=prof.srwind[1]))
   srw_eff = utils.comp2vec(*prof.srw_eff) #utils.comp2vec(*winds.sr_wind(prof, pbot=eff_inflow[0], ptop=eff_inflow[1], stu=prof.srwind[0], stv=prof.srwind[1]))
   if type(eff_inflow[0])!=np.float64:
      srw_eff = ['---','--']
   srw_6km = prof.srw_6km #utils.comp2vec(*winds.sr_wind(prof, pbot=sfc, ptop=p6km, stu=prof.srwind[0], stv=prof.srwind[1]))
   srw_8km = prof.srw_8km #utils.comp2vec(*winds.sr_wind(prof, pbot=sfc, ptop=p8km, stu=prof.srwind[0], stv=prof.srwind[1]))
   srw_cloudlayer = prof.srw_lcl_el#utils.comp2vec(*winds.sr_wind(prof, pbot=pcl.lclpres, ptop=pcl.elpres, stu=prof.srwind[0], stv=prof.srwind[1]))
   srw_ebwd = utils.comp2vec(*prof.srw_ebw) #utils.comp2vec(*winds.sr_wind(prof, pbot=eff_inflow[0], ptop=eff_inflow[1], stu=prof.srwind[0], stv=prof.srwind[1]))
   if type(eff_inflow[0])!=np.float64:
      srw_ebwd = ['---','--']
   srw_46km = utils.comp2vec(*prof.srw_4_6km)#utils.comp2vec(*winds.sr_wind(prof, pbot=p4km, ptop=p6km, stu=prof.srwind[0], stv=prof.srwind[1]))
   sfc_8km_shear = prof.sfc_8km_shear#winds.wind_shear(prof, pbot=sfc, ptop=p8km)
   sfc_6km_shear = prof.sfc_6km_shear#winds.wind_shear(prof, pbot=sfc, ptop=p6km)
   sfc_3km_shear = prof.sfc_3km_shear#winds.wind_shear(prof, pbot=sfc, ptop=p3km)
   sfc_1km_shear = prof.sfc_1km_shear#winds.wind_shear(prof, pbot=sfc, ptop=p1km)
   effective_shear = prof.eff_shear#winds.wind_shear(prof, pbot=eff_inflow[0], ptop=etop_hght)
   cloudlayer_shear = prof.lcl_el_shear#winds.wind_shear(prof,pbot= pcl.lclpres, ptop=pcl.elpres)
   srh3km = prof.srh3km#winds.helicity(prof, 0, 3000., stu = bunkers[0], stv = bunkers[1])
   srh1km = prof.srh1km#winds.helicity(prof, 0, 1000., stu = bunkers[0], stv = bunkers[1])
   stp_fixed = prof.stp_fixed#params.stp_fixed(pcl.bplus, pcl.lclhght, srh1km[0], utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])
   ship = prof.ship
   effective_srh = prof.right_esrh#winds.helicity(prof, ebot_hght, etop_hght, stu = bunkers[0], stv = bunkers[1])
   ebwd = prof.ebwd#winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
   ebwspd = prof.ebwspd
   scp = prof.right_scp
   stp_cin = prof.stp_cin#params.stp_cin(pcl.bplus, effective_srh[0], ebwspd, pcl.lclhght, pcl.bminus)
   brn_shear = pcl.brnshear

# Draw the SCP, STP, SHIP indices to the plot
# TODO: Include the color variations on this to denote intensity of the index.

#   if prof.stp_fixed is ma.masked: 
#      temp_stp_fixed = '0.0'
#   else: 
#      temp_stp_fixed = str(round(prof.stp_fixed,1))

#   if prof.stp_cin is ma.masked:
#      temp_stp_cin = '0.0'
#   else:
#      temp_stp_cin = str(round(prof.stp_cin,1))

#   if prof.right_scp is ma.masked:
#      temp_right_scp = '0.0'
#   else:
#      temp_right_scp = str(round(prof.right_scp,1))

#   data = np.array([["Supercell =",temp_right_scp],
#                       ["STP (cin) =",temp_stp_cin],
#                       ["STP (fix) =",temp_stp_fixed]])

   data = np.array([["Supercell =",str(np.ma.around(prof.right_scp,1))],
                       ["STP (cin) =",str(np.ma.around(prof.stp_cin,1))],
                       ["STP (fix) =",str(np.ma.around(prof.stp_fixed,1))]])

   data = np.array([["Supercell =",str(np.ma.around(prof.right_scp,1))],
                       ["STP (cin) =",str(np.ma.around(prof.stp_cin,1))],
                       ["STP (fix) =",str(np.ma.around(prof.stp_fixed,1))],
                       ["SHIP =",str(np.ma.around(prof.ship,1))]])
   '''
   for i in range(data.shape[0]):
      for j in range(data.shape[1]):
         d = float(data[i,1])
         if i == 0:
            if d >= 19.95:
               c = MAGENTA
            elif d >= 11.95:
               c = RED
            elif d >= 1.95:
               c = YELLOW
            elif d >= .45:
               c = WHITE
            elif d >= -.45:
               c = LBROWN
            elif d < -0.45:
               c = CYAN
         elif i == 1:
            if d >= 8:
               c = MAGENTA
            elif d >= 4:
               c = RED
            elif d >= 2:
               c = YELLOW
            elif d >= 1:
               c = WHITE
            elif d >= .5:
               c = LBROWN
            elif d < .5:
               c = DBROWN
               stpCinColor = c
         elif i == 2:
            if d >= 7:
               c = MAGENTA
            elif d >= 5:
               c = RED
            elif d >= 2:
               c = YELLOW
            elif d >= 1:
               c = WHITE
            elif d >= 0.5:
               c = LBROWN
            else:
               c = DBROWN
         elif i == 3:
            if d >= 5:
               c = MAGENTA
            elif d >= 2:
               c = RED
            elif d >= 1:
               c = YELLOW
            elif d >= .5:
               c = WHITE
            else:
               c = DBROWN
         if (j % 2 == 0):
            ax2.annotate(data[i,j], (x_list[j], y_list[i]), xycoords="axes fraction",
                         fontsize=12, va="top", ha="right", color=c, weight='bold')
         else:
            ax2.annotate(data[i,j], (x_list[j], y_list[i]), xycoords="axes fraction",
                       fontsize=12, va="top", ha="left", color=cb_colors.gray7, weight='bold')
   '''    
# Draw the kinematic inset on the plot
#   box = ax3.add_patch(patches.Rectangle((9.0, -0.04), 60.0, -0.025,fill=False,linewidth=1,edgecolor="w",zorder=4))
#   box.set_clip_on(False)
   x_list = np.array([0.60, 0.84, 0.97, 1.08, 1.18])-0.12
   y_list = [-0.045]
   y_list.extend(np.arange(-.07,-0.12,-.0225).tolist())
   y_list.extend(np.arange(-.145,-0.22,-.0225).tolist())
   y_list.extend(np.arange(-.2425,-0.27,-.0225).tolist())
   y_list.extend(np.arange(-0.295, -.4, -0.0225).tolist())

   A = ["SFC-1km", srh1km[0], utils.comp2vec(sfc_1km_shear[0],sfc_1km_shear[1])[1]]
   A2 = [mean_1km_comp, srw_1km]
   B = ["SFC-3km", srh3km[0], utils.comp2vec(sfc_3km_shear[0],sfc_3km_shear[1])[1]]
   B2 = [mean_3km_comp, srw_3km]
   C = ["Eff Inflow Layer", effective_srh[0], utils.comp2vec(effective_shear[0],effective_shear[1])[1]]
   C2 = [mean_eff_comp, srw_eff]
#   D = ["SFC-6km", "", utils.comp2vec(sfc_6km_shear[0],sfc_6km_shear[1])[1]]
#   D2 = [mean_6km_comp, srw_6km]
#   E = ["SFC-8km", "", utils.comp2vec(sfc_8km_shear[0],sfc_8km_shear[1])[1]]
#   E2 = [mean_8km_comp, srw_8km]
#   F = ["LCL-EL (CLoud Layer)", "", utils.comp2vec(cloudlayer_shear[0],cloudlayer_shear[1])[1]]
#   F2 = [mean_cloudlayer_comp, srw_cloudlayer]
#   G = ["Eff Shear (EBWD)", "", utils.comp2vec(ebwd[0],ebwd[1])[1]]
#   G2 = [mean_ebwd_comp, srw_ebwd]
#   H = ["BRN Shear (m2/s2)", "", brn_shear, "", ""]
#   I = ["4-6km SR Wind", ""]
#   I2 = [str(int(round(srw_46km[0],0)))+"/"+str(int(round(srw_46km[1],0)))]
#   I3 = ["", ""]
#   J = ["...Storm Motion Vectors...", "", "", "", ""]
#   K = ["Bunkers Right", ""]
#   K2 = [str(int(round(bunkers_rt[0],0)))+"/"+str(int(round(bunkers_rt[1],0)))]
#   K3 = ["", ""]
#   L = ["Bunkers Left", ""]
#   L2 = [str(int(round(bunkers_lf[0],0)))+"/"+str(int(round(bunkers_lf[1],0)))]
#   L3 = ["", ""]
#   M = ["Corfidi Downshear", ""]
#   M2 = [str(int(round(corfidi_dn[0],0)))+"/"+str(int(round(corfidi_dn[1],0)))]
#   M3 = ["", ""]
#   N = ["Corfidi Upshear", ""]
#   N2 = [str(int(round(corfidi_up[0],0)))+"/"+str(int(round(corfidi_up[1],0)))]
#   N3 = ["", ""]

   data = np.array([np.array(["", "SRH (m2/s2)", "Shear (kt)", "MnWind", "SRW"]),
                    np.array([ str(int(round(A[i],0))) if (type(A[i])== np.float64) else str(A[i]) for i in range(len(A)) ]+\
                       [ str(int(np.ma.around(A2[i][0],0)))+"/"+str(int(np.ma.around(A2[i][1],0))) if (type(A2[i][0])== np.ma.core.MaskedArray) else str(A2[i][0])+"/"+str(A2[i][1]) for i in range(len(A2)) ]),

                    np.array([ str(int(round(B[i],0))) if (type(B[i])== np.float64) else str(B[i]) for i in range(len(B)) ]+\
                        [ str(int(np.ma.around(B2[i][0],0)))+"/"+str(int(np.ma.around(B2[i][1],0))) if (type(B2[i][0])== np.ma.core.MaskedArray) else str(B2[i][0])+"/"+str(B2[i][1]) for i in range(len(B2)) ]),

                    np.array([ str(int(np.ma.around(C[i],0))) if (type(C[i])== np.float64) else str(C[i]) for i in range(len(C)) ]+\
                        [ str(int(np.ma.around(C2[i][0],0)))+"/"+str(int(np.ma.around(C2[i][1],0))) if (type(C2[i][0])== np.ma.core.MaskedArray) else str(C2[i][0])+"/"+str(C2[i][1]) for i in range(len(C2)) ])]) #,

#                    np.array([ str(int(round(D[i],0))) if (type(D[i])== np.float64) else str(D[i]) for i in range(len(D)) ]+\
#                        [ str(int(round(D2[i][0],0)))+"/"+str(int(round(D2[i][1],0))) if (type(D2[i][0])== np.ma.core.MaskedArray) else str(D2[i][0])+"/"+str(D2[i][1]) for i in range(len(D2)) ]),

#                    np.array([ str(int(round(E[i],0))) if (type(E[i])== np.float64) else str(E[i]) for i in range(len(E)) ]+\
#                        [ str(int(round(E2[i][0],0)))+"/"+str(int(round(E2[i][1],0))) if (type(E2[i][0])== np.ma.core.MaskedArray) else str(E2[i][0])+"/"+str(E2[i][1]) for i in range(len(E2)) ]),

#                    np.array([ str(int(round(F[i],0))) if (type(F[i])== np.float64) else str(F[i]) for i in range(len(F)) ]+\
#                        [ str(int(round(F2[i][0],0)))+"/"+str(int(round(F2[i][1],0))) if (type(F2[i][0])== np.ma.core.MaskedArray) else str(F2[i][0])+"/"+str(F2[i][1]) for i in range(len(F2)) ]),

#                    np.array([ str(int(round(G[i],0))) if (type(G[i])== np.float64) else str(G[i]) for i in range(len(G)) ]+\
#                        [ str(int(round(G2[i][0],0)))+"/"+str(int(round(G2[i][1],0))) if (type(G2[i][0])== np.ma.core.MaskedArray) else str(G2[i][0])+"/"+str(G2[i][1]) for i in range(len(G2)) ]),

#                    np.array([ str(int(round(H[i],0))) if (type(H[i])== np.float64) else str(H[i]) for i in range(len(H)) ]),
#                    np.array(I+I2+I3),
#                    np.array(J),
#                    np.array(K+K2+K3),
#                    np.array(L+L2+L3),
#                    np.array(M+M2+M3),
#                    np.array(N+N2+N3)])
   #x_list = np.array([0, 0.08, 0.17, 0.25, 0.33, 0.39, 0.47])-0.05
   y_list = [-0.045, -0.07, -0.0925, -0.115, -0.1375]
   for i in range(data.shape[0]):
      for j in range(data[0].shape[0]):
         if j>0:
            ax2.annotate(data[i][j], (x_list[j], y_list[i]), xycoords="axes fraction",
                         fontsize=12, va="top", ha="right", color=cb_colors.gray7, weight='bold')
         else:
            ax2.annotate(data[i][j], (x_list[j], y_list[i]), xycoords="axes fraction",
                         fontsize=12, va="top", ha="left", color=cb_colors.gray7, weight='bold')

#    wind_1km = utils.vec2comp(prof.wind1km[0], prof.wind1km[1])
#    wind_6km = utils.vec2comp(prof.wind6km[0], prof.wind6km[1])
#    wind_barbs = ax2.barbs(58, 2200, wind_1km[0], wind_1km[1], color='#AA0000', zorder=3, lw=1.25,length=9)
#    wind_barbs.set_clip_on(False)
#    wind_barbs = ax2.barbs(58, 2200, wind_6km[0], wind_6km[1], color='#0A74C6', zorder=3, lw=1.25,length=9)
#    wind_barbs.set_clip_on(False)
#    ax2.annotate("1km & 6km AGL\nWind Barbs", (1.08,-0.37), xycoords="axes fraction", fontsize=12, va="top", ha="center", color='#0A74C6', weight='bold')

# Draw the CAPE vs. SRH Scatter
   #ax_EFSTP = fig.add_axes([0.74625, 0.0229, 0.20375, 0.2507], frameon=False)
   #x_EFSTP = fig.add_axes([0.72625, 0.0229, 0.20375, 0.2507], frameon=False)
   ax_EFSTP = fig.add_axes([0.625, 0.05, 0.35, 0.3], frameon=False)
   
   for member_no,c in zip(np.arange(1,19,1), np.tile([cb_colors.orange6, cb_colors.orange6, cb_colors.green6, cb_colors.green6, cb_colors.purple6, cb_colors.purple6], (3,1))):
      memidx = [0,10,11,12,13,14,15,16,17,1,2,3,4,5,6,7,8,9]
      ax_EFSTP.scatter(xpts[memidx[member_no-1],:,:].ravel(),ypts[memidx[member_no-1],:,:].ravel(), color=c,marker='o',s=5,alpha=0.7)
   ax_EFSTP.set_xticks(np.arange(0,5001,1000))
   ax_EFSTP.set_yticks(np.arange(0,601,100))
   #ax_EFSTP.set_xticklabels(np.arange(0,5001,1000),color='k',fontsize=12)
   #ax_EFSTP.set_yticklabels(np.arange(0,601,100),color='k',fontsize=12)
   ax_EFSTP.set_xlabel('MLCAPE',weight='bold',fontsize=14)
   ax_EFSTP.set_ylabel('0-1km SRH',weight='bold',fontsize=14)
   ax_EFSTP.tick_params(axis='both',labelsize=12,labelcolor='k')
   ax_EFSTP.set_xlim(-200,5000)
   ax_EFSTP.set_ylim(-100,600)
   #ax_EFSTP.tick_params(direction='in', axis='x', which='major', colors=cb_colors.gray4,length=0,width=1.5,size=12)#pad=-10
   #ax_EFSTP.tick_params(direction='in', axis='y', which='major', colors=cb_colors.gray4,length=0,width=1.5,size=12)#pad=-23
   ax_EFSTP.grid(color=cb_colors.gray4, linestyle='--', dashes=(3,3), alpha=0.75, zorder=0, linewidth=1.25)
   ax_EFSTP.text(0.8, 0.95, 'YSU', color=cb_colors.orange6, transform=ax_EFSTP.transAxes, fontsize=16,weight='bold')
   ax_EFSTP.text(0.8, 0.89, 'MYJ', color=cb_colors.green6, transform=ax_EFSTP.transAxes, fontsize=16,weight='bold')
   ax_EFSTP.text(0.8, 0.83, 'MYNN', color=cb_colors.purple6, transform=ax_EFSTP.transAxes, fontsize=16,weight='bold')
   box = ax_EFSTP.add_patch(patches.Rectangle((0, -1), 13, 13,fill=False,linewidth=2,edgecolor=cb_colors.gray4,zorder=10))
   box.set_clip_on(False)
 
   ax_EFSTP.annotate("0-1 km SRH vs. 100-mb MLCAPE", (0.5, 1.075), xycoords="axes fraction",
                 fontsize=14, va="center", ha="center", color=cb_colors.gray7, weight='bold')
   ax_EFSTP.annotate("(9 km neighborhood)", (0.5, 1.025), xycoords="axes fraction",
                 fontsize=12, va="center", ha="center", color=cb_colors.gray7, weight='bold')


   plt.savefig(figname, facecolor=fig.get_facecolor())#, edgecolor=None)
   #stop
