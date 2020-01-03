#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
#from skimage.morphology import label
#from skimage.measure import regionprops
import math
from math import radians, tan, sin, cos, pi, atan, sqrt, pow, asin, acos
import pylab as P
import numpy as np
from numpy import NAN
import sys
import netCDF4
from optparse import OptionParser
import pickle as pl
from netcdftime import utime
import os
import time
import ctables
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input directory of summary files to plot")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
parser.add_option("-m", dest="mapname", type="string", help = "Path to Pickled Basemap instance")
parser.add_option("-n", dest="name", type="string", help = "Name of variable to process")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-f", dest="fcst_nt", type="int", help = "Total forecast timesteps")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.outdir == None) or (options.mapname == None) or (options.name == None) or (options.t == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   outdir = options.outdir
   mapname = options.mapname
   name = options.name
   t = options.t
   fcst_nt = options.fcst_nt

print 'PLOT SWATH RAINFALL: ', t, name
print 'TOTAL TIMESTEPS IN THIS FORECAST: ', fcst_nt

#################################### Get input file:  #####################################################

summary_files_temp = os.listdir(summary_dir)
timestep = str(t)
if (len(timestep) == 1):
   timestep = '0' + timestep

summary_files_temp.sort()
# Determine if its a 36 or 72 timestep forecast and find all PCP files #

if (fcst_nt == 36):
   pcp_found_file = 0
   for f, file in enumerate(summary_files_temp):
      if ((file[-28:-25] == 'PCP') and (file[-24:-22] == timestep)):
         infile = os.path.join(summary_dir, file)
         if ((file[-24:-22] == '36')):
            print 'All %s PCP files for this fcst are present.' % fcst_nt
            pcp_found_file = 1
            break
         else:
            print file[-24:-22]
            #print 'Not all %s PCP files are present for this fcst. Trying again in 30s' % fcst_nt
            #time.sleep(30)

elif (fcst_nt == 72):
    pcp_found_file = 0
    for f, file in enumerate(summary_files_temp):
       if ((file[-28:-25] == 'PCP') and (file[-24:-22] == timestep)):
          infile = os.path.join(summary_dir, file)
          if ((file[-24:-22] == '36')):
             print 'Have enough for 3-hr forecast, but not for 6-hr just yet.'
             continue
          if ((file[-24:-22] == '72')):
             print 'All %s PCP files for this fcst are present.' % fcst_nt
             pcp_found_file = 1
             break
          else:
             print 'Not all %s PCP files are present for this fcst. Trying again in 30s' % fcst_nt
             time.sleep(30)

#       while

#pcp_found_file = 0
#for i in range(0,100):
#   for f, file in enumerate(summary_files_temp):
#      if ((file[-28:-25] == 'SMR') and (file[-24:-22] == timestep)):
#         infile = os.path.join(summary_dir, file)
#         pcp_found_file = 1
#   if (pcp_found_file == 1):
#      break
#   else:
#      print 'PCP file total incomplete. Trying again in 30s'
#      time.sleep(30)


#   try:
#      if ((file[-28:-25] == 'PCP') and (file[-24:-22] == timestep)):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
#         infile = os.path.join(summary_dir, file)
#   except NameError:
#      print 'file is not there yet. Trying again in 30s'
#      time.sleep(30)

print 'Matched PCP File: ', timestep, '   ', infile
#################################### Basemap Variables:  #####################################################

resolution 	= 'h'
area_thresh 	= 1000.

damage_files = '' #['/scratch/skinnerp/2018_newse_post/damage_files/extractDamage_new/extractDamagePolys'] #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### User-Defined Variables:  #####################################################

domain          = 'full'        #vestigial variable that's still needed in the plotting subroutines ...
edge            = 7 		#number of grid points to remove from near domain boundaries
thin		= 6		#thinning factor for quiver values (e.g. 6 means slice every 6th grid point)

neighborhood 	= 15		#grid point radius of prob matched mean neighborhood

plot_alpha 	= 0.55		#transparency value for filled contour plots

comp_dz_thresh             = 45.                #40 dBZ
rain_1_thresh              = 0.5                #0.5 inches
rain_1_3hr_thresh	   = 1.
rain_1_6hr_thresh	   = 2.
rain_2_thresh              = 1.                 #1 inch
rain_2_3hr_thresh	   = 2.
rain_2_6hr_thresh	   = 3.
rain_3_thresh              = 2.                 #2 inches
rain_3_3hr_thresh	   = 3.
rain_3_6hr_thresh	   = 5.

rain_4_thresh        = 3.                       #3 inches

rain_5_thresh        = 5.                       #5 inches

#soil_moisture_thresh       = .95                #95% soil moisture

################################### Set plot object values for different variables:  #####################################################

if (name == 'rain'):
   var_label       = 'Accumulated Rainfall'
   var_units       = 'inches'
elif (name == 'rain_sat'):
   var_label       = 'Accumulated Rainfall on Sat. Soil'
   var_units       = 'inches'
elif (name == 'soil_moisture'):
   var_label       = 'Top Layer Soil Moisture'
   var_units       = '%'
#################################### Contour Levels:  #####################################################

prob_levels             = np.arange(0.1,1.1,0.1)                #(%)
rain_levels		= [0.01, 0.1, 0.25, 0.50, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]
#rain_levels		= [0, 0.1, 0.3, 0.5, 1., 2., 3., 4., 5., 7., 10.]
soil_moisture_levels    = np.arange(.35,1.,.05)
dz_levels_nws           = np.arange(20.0,80.,5.)                #(dBZ)

pmm_dz_levels           = [35., 50.]                            #(dBZ)
pmm_dz_colors_gray      = [cb_colors.gray8, cb_colors.gray8]    #gray contours

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

base_plot = web_plot('',									\
                     '', 									\
                     'Probability Matched Mean - Composite Reflectivity (dBZ)',                 \
                     cb_colors.gray6,    							\
                     '',           						         	\
                     pmm_dz_levels,       							\
                     '',                  							\
                     '',                  							\
                     pmm_dz_colors_gray,     							\
                     cb_colors.purple4,              						\
                     'none',              							\
                     cb_colors.wz_cmap_extend,            							\
                     'max',               							\
                   plot_alpha,               \
                   neighborhood)

rain_plot = web_plot('',									\
		     '',									\
		     'Probability Matched Mean - Composite Reflectivity (dBZ)',			\
		     cb_colors.gray6,								\
		     '',									\
                     pmm_dz_levels,								\
		     '',									\
		     '',									\
		     pmm_dz_colors_gray,							\
		     cb_colors.red9,								\
		     'none',									\
		     cb_colors.rain_cmap,							\
		     'max',									\
		     plot_alpha,		     \
		     neighborhood)

if (name == 'rain'):
   rain_plot.var1_levels = rain_levels
elif (name == 'rain_sat'):
   rain_plot.var1_levels = rain_levels
elif (name == 'soil_moisture'):
   base_plot.var1_levels = soil_moisture_levels
elif (name == 'comp_dz'):
   base_plot.var1_levels = dz_levels_nws
   base_plot.var2_tcolor = 'none'
   base_plot.var2_levels = [100., 110.]       #hack so nothing is plotted
   base_plot.var2_colors = ['none', 'none']
   base_plot.cmap = cb_colors.nws_dz_cmap
   base_plot.over_color = cb_colors.c75
   base_plot.alpha = 0.8

prob_plot = web_plot('',                                                                     \
                   '',		                                                             \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                \
                   cb_colors.gray6,                                                         \
                   prob_levels,                                                              \
                   pmm_dz_levels,                                                            \
                   '',                                                                       \
                   '',                                                                       \
                   pmm_dz_colors_gray,                                                          \
                   'none',	                                                             \
                   'none',		                                                     \
                   cb_colors.wz_cmap,                                                                 \
                   'neither',                                                                \
                   plot_alpha,               \
                   neighborhood)

######################################################################################################
#################################### Read Data:  #####################################################
######################################################################################################

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(infile, "r")
   print "Opening %s \n" % infile
except:
   print "%s does not exist! \n" %infile
   sys.exit(1)

######################### Read Attributes: ####################################

dx = fin.DX                                             #east-west grid spacing (m)
dy = fin.DY                                             #north-south grid spacing (m)
cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
true_lat_1 = fin.TRUE_LAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
true_lat_2 = fin.TRUE_LAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)
projection = fin.PROJECTION                             #domain projection
init_time_seconds = fin.INIT_TIME_SECONDS               #initialization time in seconds from 0000 UTC of case
valid_time_seconds = fin.VALID_TIME_SECONDS             #valid time in seconds from 0000 UTC of case
forecast_timestep = fin.FORECAST_TIME_STEP              #index of forecast time step

######################### Parse Init/Valid times from filename: ####################################

year = infile[-21:-17]
month = infile[-17:-15]
day = infile[-15:-13]
init_hour = infile[-12:-10]
init_min = infile[-10:-8]
valid_hour = infile[-7:-5]
valid_min = infile[-5:-3]

if ((int(valid_hour) < 12) and (int(init_hour) > 18)):
   temp_day = int(day) + 1
   valid_day = str(temp_day)
   if (len(valid_day) == 1):
      valid_day = '0' + valid_day
else:
   valid_day = day

init_label = 'Init: ' + year + '-' + month + '-' + day + ', ' + init_hour + init_min + ' UTC'
valid_label = 'Valid: ' + year + '-' + month + '-' + valid_day + ', ' + valid_hour + valid_min + ' UTC'
######################### Set domain: ####################################

xlat = fin.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][:,:]                     #longitude (dec deg; Lambert conformal)

xlat = xlat[edge:-edge,edge:-edge]
xlon = xlon[edge:-edge,edge:-edge]

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

######################### Read PMM composite reflectivity: ####################################

pmm_dz = fin.variables["comp_dz_pmm"][edge:-edge,edge:-edge]

################################# Make Figure Template: ###################################################

print 'basemap part'

#Load pickled basemap instance for faster plotting:

fig, ax1, ax2, ax3 = create_fig_nomap()

map_temp = pl.load(open(mapname, 'rb'))

P.sca(ax1)   #make sure boundaries are plotted in the right ax
map = mymap_boundaries(map_temp, damage_files)

#map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

###

if (name == 'rain'):
   var_90 = fin.variables["rain_90"][edge:-edge,edge:-edge]
   var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#   var_90 = np.ma.masked_less(var_90, 0.01)

   var_max = fin.variables["rain_max"][edge:-edge,edge:-edge]
   var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#   var_max = np.ma.masked_less(var_max, 0.01)

   var_med = fin.variables["rain_med"][edge:-edge,edge:-edge]
   var_med = np.where(var_med <= 0.01, -1., var_med)
#    var_med = np.ma.masked_less(var_med, 0.01)

   var_9km = fin.variables["rain_half_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_half_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_half_prob_27km"][edge:-edge,edge:-edge]

   rain_plot.name = name + '_90'
   rain_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
   env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')

   rain_plot.name = name + '_max'
   rain_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
   env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')
#
   rain_plot.name = name + '_med'
   rain_plot.var1_title = 'Ens. Median Value of ' + var_label +  ' (' + var_units + ')'
   env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_med[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')

   prob_plot.name = name + '_prob_3km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_one_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_one_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_one_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_two_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_two_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_two_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_three_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_three_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_three_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_three'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_4_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_three'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_4_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_three'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_4_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_five_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_five_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_five_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_five'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_5_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_five'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_5_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_five'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_5_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

#
   if ((t > 0) and ((t % 12) == 0)):
      var_90 = fin.variables["rain_90_hourly"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
##      var_90 = np.ma.masked_less(var_90, 0.01)
#
      var_max = fin.variables["rain_max_hourly"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
##      var_max = np.ma.masked_less(var_max, 0.01)
#
      var_med = fin.variables["rain_med_hourly"][edge:-edge,edge:-edge]
      var_med = np.where(var_med <= 0.01, -1., var_med)
#      var_med = np.ma.masked_less(var_med, 0.01)

      var_9km = fin.variables["rain_half_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_half_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_half_prob_27km_hourly"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_hourly'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_hourly'
      rain_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')
#
      rain_plot.name = name + '_med_hourly'
      rain_plot.var1_title = 'Ens. Median Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_med[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_one_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_one_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_one_prob_27km_hourly"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_two_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_two_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_two_prob_27km_hourly"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

   if ((t > 0) and ((t % 36) == 0)):
      var_90 = fin.variables["rain_90_3hr"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#      var_90 = np.ma.masked_less(var_90, 0.01)

      var_max = fin.variables["rain_max_3hr"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#      var_max = np.ma.masked_less(var_max, 0.01)

      var_med = fin.variables["rain_med_3hr"][edge:-edge,edge:-edge]
      var_med = np.where(var_med <= 0.01, -1., var_med)

      var_9km = fin.variables["rain_half_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_half_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_half_prob_27km_3hr"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_3hr'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of 3-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_3hr'
      rain_plot.var1_title = 'Ens. Max Value of 3-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_med_3hr'
      rain_plot.var1_title = 'Ens. Median Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_med[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_one_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_one_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_one_prob_27km_3hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_two_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_two_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_two_prob_27km_3hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')


   if ((t > 0) and ((t % 72) == 0)):
      var_90 = fin.variables["rain_90_6hr"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#      var_90 = np.ma.masked_less(var_90, 0.01)

      var_max = fin.variables["rain_max_6hr"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#      var_max = np.ma.masked_less(var_max, 0.01)

      var_med = fin.variables["rain_med_6hr"][edge:-edge,edge:-edge]
      var_med = np.where(var_med <= 0.01, -1., var_med)

      var_9km = fin.variables["rain_half_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_half_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_half_prob_27km_6hr"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_6hr'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of 6-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_6hr'
      rain_plot.var1_title = 'Ens. Max Value of 6-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_med_6hr'
      rain_plot.var1_title = 'Ens. Median Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_med[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_one_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_one_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_one_prob_27km_6hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_two_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_two_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_two_prob_27km_6hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')
###

###

if (name == 'rain_sat'):
   var_90 = fin.variables["rain_sat_90"][edge:-edge,edge:-edge]
   var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#   var_90 = np.ma.masked_less(var_90, 0.01)

   var_max = fin.variables["rain_sat_max"][edge:-edge,edge:-edge]
   var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#   var_max = np.ma.masked_less(var_max, 0.01)

   var_9km = fin.variables["rain_sat_half_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_sat_half_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_sat_half_prob_27km"][edge:-edge,edge:-edge]

   rain_plot.name = name + '_90'
   rain_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
   env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')

   rain_plot.name = name + '_max'
   rain_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
   env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')

   prob_plot.name = name + '_prob_3km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_half'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_sat_one_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_sat_one_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_sat_one_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_one'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   var_9km = fin.variables["rain_sat_two_prob_3km"][edge:-edge,edge:-edge]
   var_15km = fin.variables["rain_sat_two_prob_15km"][edge:-edge,edge:-edge]
   var_27km = fin.variables["rain_sat_two_prob_27km"][edge:-edge,edge:-edge]

   prob_plot.name = name + '_prob_3km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_15km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

   prob_plot.name = name + '_prob_27km_two'
   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
   prob_plot.var1_title = prob_plot.var1_title + var_units
   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')
#
   if ((t > 0) and ((t % 12) == 0)):
      var_90 = fin.variables["rain_sat_90_hourly"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#      var_90 = np.ma.masked_less(var_90, 0.01)
#
      var_max = fin.variables["rain_sat_max_hourly"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#      var_max = np.ma.masked_less(var_max, 0.01)
#
      var_9km = fin.variables["rain_sat_half_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_half_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_half_prob_27km_hourly"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_hourly'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_hourly'
      rain_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_1_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_one_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_one_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_one_prob_27km_hourly"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_two_prob_3km_hourly"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_two_prob_15km_hourly"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_two_prob_27km_hourly"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_hourly'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

   if ((t > 0) and ((t % 36) == 0)):
      var_90 = fin.variables["rain_sat_90_3hr"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#      var_90 = np.ma.masked_less(var_90, 0.01)

      var_max = fin.variables["rain_sat_max_3hr"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#      var_max = np.ma.masked_less(var_max, 0.01)

      var_9km = fin.variables["rain_sat_half_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_half_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_half_prob_27km_3hr"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_3hr'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of 3-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_3hr'
      rain_plot.var1_title = 'Ens. Max Value of 3-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3_hrthresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_1_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_one_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_one_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_one_prob_27km_3hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_2_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_two_prob_3km_3hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_two_prob_15km_3hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_two_prob_27km_3hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_3hr'
      prob_plot.var1_title = 'Probability of 3-hr ' + var_label + ' > %s ' % str(rain_3_3hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

   if ((t > 0) and ((t % 72) == 0)):
      var_90 = fin.variables["rain_sat_90_6hr"][edge:-edge,edge:-edge]
      var_90 = np.where(var_90 <= 0.01, -1., var_90)   #hack to not contour 0.00 inches of rain
#      var_90 = np.ma.masked_less(var_90, 0.01)

      var_max = fin.variables["rain_sat_max_6hr"][edge:-edge,edge:-edge]
      var_max = np.where(var_max <= 0.01, -1., var_max)   #hack to not contour 0.00 inches of rain
#      var_max = np.ma.masked_less(var_max, 0.01)

      var_9km = fin.variables["rain_sat_half_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_half_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_half_prob_27km_6hr"][edge:-edge,edge:-edge]

      rain_plot.name = name + '_90_6hr'
      rain_plot.var1_title = 'Ens. 90th Percentile Value of 6-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      rain_plot.name = name + '_max_6hr'
      rain_plot.var1_title = 'Ens. Max Value of 6-hr ' + var_label +  ' (' + var_units + ')'
      env_plot_rain(map, fig, ax1, ax2, ax3, x, y, rain_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')

      prob_plot.name = name + '_prob_3km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_half_6hr'
      prob_plot.var1_title = 'Probability of 6-hr ' + var_label + ' > %s ' % str(rain_1_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_one_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_one_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_one_prob_27km_6hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_one_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_one_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_one_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_2_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      var_9km = fin.variables["rain_sat_two_prob_3km_6hr"][edge:-edge,edge:-edge]
      var_15km = fin.variables["rain_sat_two_prob_15km_6hr"][edge:-edge,edge:-edge]
      var_27km = fin.variables["rain_sat_two_prob_27km_6hr"][edge:-edge,edge:-edge]

      prob_plot.name = name + '_prob_3km_two_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_15km_two_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

      prob_plot.name = name + '_prob_27km_two_6hr'
      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(rain_3_6hr_thresh)
      prob_plot.var1_title = prob_plot.var1_title + var_units
      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

###

#if (name == 'soil_moisture'):
#   var_90 = fin.variables["soil_moisture_90"][edge:-edge,edge:-edge]
#   var_max = fin.variables["soil_moisture_max"][edge:-edge,edge:-edge]
#   var_9km = fin.variables["soil_moisture_prob_3km"][edge:-edge,edge:-edge]
#   var_15km = fin.variables["soil_moisture_prob_15km"][edge:-edge,edge:-edge]
#   var_27km = fin.variables["soil_moisture_prob_27km"][edge:-edge,edge:-edge]
#
#   base_plot.name = name + '_90'
#   base_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
#   env_plot(map, fig, ax1, ax2, ax3, x, y, base_plot, var_90[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')
#
#   base_plot.name = name + '_max'
#   base_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
#   env_plot(map, fig, ax1, ax2, ax3, x, y, base_plot, var_max[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False', showmax='True')
#
#   prob_plot.name = name + '_prob_9km'
#   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#   prob_plot.var1_title = prob_plot.var1_title + var_units
#   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')
#
#   prob_plot.name = name + '_prob_15km'
#   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#   prob_plot.var1_title = prob_plot.var1_title + var_units
#   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')
#
#   prob_plot.name = name + '_prob_27km'
#   prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#   prob_plot.var1_title = prob_plot.var1_title + var_units
#   env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')
#
#   if ((t > 0) and ((t % 12) == 0)):
#      var_90 = fin.variables["soil_moisture_90_hourly"][edge:-edge,edge:-edge]
#      var_max = fin.variables["soil_moisture_max_hourly"][edge:-edge,edge:-edge]
#      var_9km = fin.variables["soil_moisture_prob_3km_hourly"][edge:-edge,edge:-edge]
#      var_15km = fin.variables["soil_moisture_prob_15km_hourly"][edge:-edge,edge:-edge]
#      var_27km = fin.variables["soil_moisture_prob_27km_hourly"][edge:-edge,edge:-edge]
#
#      base_plot.name = name + '_90_hourly'
#      base_plot.var1_title = 'Ens. 90th Percentile Value of ' + var_label +  ' (' + var_units + ')'
#      env_plot(map, fig, ax1, ax2, ax3, x, y, base_plot, var_90[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')
#
#      base_plot.name = name + '_max_hourly'
#      base_plot.var1_title = 'Ens. Max Value of ' + var_label +  ' (' + var_units + ')'
#      env_plot(map, fig, ax1, ax2, ax3, x, y, base_plot, var_max[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False', showmax='True')
#
#      prob_plot.name = name + '_prob_9km_hourly'
#      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#      prob_plot.var1_title = prob_plot.var1_title + var_units
#      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_9km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')
#
#      prob_plot.name = name + '_prob_15km_hourly'
#      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#      prob_plot.var1_title = prob_plot.var1_title + var_units
#      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_15km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')
#
#      prob_plot.name = name + '_prob_27km_hourly'
#      prob_plot.var1_title = 'Probability of ' + var_label + ' > %s ' % str(w_up_thresh)
#      prob_plot.var1_title = prob_plot.var1_title + var_units
#      env_plot(map, fig, ax1, ax2, ax3, x, y, prob_plot, var_27km[:,:], pmm_dz[:,:], (t/12), init_label, valid_label, domain, outdir, '', '', '', 1, 0, spec='False', quiv='False')

###

fin.close()
del fin

