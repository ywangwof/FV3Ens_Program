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
import time as timeit
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
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.outdir == None) or (options.mapname == None) or (options.t == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   outdir = options.outdir
   mapname = options.mapname
   t = options.t

#path to pickled map instance for plotting: 
#mapname = '/scratch2/patrick.skinner/images/map.pickle'

#################################### User-Defined Variables:  #####################################################

domain                     = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge                       = 7 		#number of grid points to remove from near domain boundaries

radius_max                 = 3                  #grid point radius for maximum value filter (3x3 square neighborhood)
radius_gauss               = 2                  #grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

plot_alpha          	   = 0.55		#transparency value for filled contour plots

#################################### Basemap Variables:  #####################################################

resolution 	= 'h'
area_thresh 	= 1000.

damage_files = '' #['/scratch/skinnerp/2018_newse_post/damage_files/extractDamage_new/extractDamagePolys'] #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Contour Levels:  #####################################################

dz_levels_nws       	= np.arange(20.0,80.,5.)		#(dBZ)

uh_2to5_levels 		= [100., 300., 7000.]				#(dBZ) 
pmm_dz_colors_gray	= [cb_colors.gray5, cb_colors.gray8, 'none']	#gray contours

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

dz_plot = web_plot('',                   \
                   '',			\
                   '',                   \
                   cb_colors.gray6,      \
                   dz_levels_nws,            \
                   uh_2to5_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.c75,		\
                   'none',		\
                   cb_colors.nws_dz_cmap,              \
                   'max',                \
                   0.35,               \
                   neighborhood)  

############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18

ens_files = []
summary_files_temp = os.listdir(summary_dir)

for f, file in enumerate(summary_files_temp):
   if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file)

ens_files.sort()

for tt in range(0, t+1):
   ens_file = ens_files[tt]

   infile = os.path.join(summary_dir, ens_file)

   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (tt == 0):

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

######################### Set domain: ####################################

      xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

      sw_lat_full = xlat[0,0]
      sw_lon_full = xlon[0,0]
      ne_lat_full = xlat[-1,-1]
      ne_lon_full = xlon[-1,-1]

######################### Read 2-5 km UH and composite reflectivity: ####################################

      uh_2to5                            = fin.variables["uh_2to5"][:,edge:-edge,edge:-edge]
      comp_dz                            = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]
   else:
      uh_2to5                            = np.where((fin.variables["uh_2to5"][:,edge:-edge,edge:-edge] > uh_2to5), fin.variables["uh_2to5"][:,edge:-edge,edge:-edge], uh_2to5)
      comp_dz                            = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]

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

   fin.close()
   del fin

#################################### Calc Swath:  #####################################################

print 'swath part'

########################################################################################################
### If updraft helicity, vertical vorticity, wind speed, graupel max or updraft plot, run maximum value and convolution filters
### over the raw data to spread and smooth the data
########################################################################################################

uh_2to5_convolve_temp = uh_2to5 * 0.
uh_2to5_convolve = uh_2to5 * 0.

kernel = gauss_kern(radius_gauss)

for n in range(0, uh_2to5.shape[0]):
   uh_2to5_convolve_temp[n,:,:] = get_local_maxima2d(uh_2to5[n,:,:], radius_max)
   uh_2to5_convolve[n,:,:] = signal.convolve2d(uh_2to5_convolve_temp[n,:,:], kernel, 'same')

################################# Make Figure Template: ###################################################

print 'basemap part'

#Load pickled basemap instance for faster plotting: 

fig, ax1, ax2, ax3 = create_fig_nomap()

map_temp = pl.load(open(mapname, 'rb'))

P.sca(ax1)
map = mymap_boundaries(map_temp, damage_files)

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])
xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

######################## Reflectivity Plots: #####################

for n in range(0, uh_2to5_convolve.shape[0]): 
   if ((n+1) < 10): 
      member = '0' + str(n+1)
   else: 
      member = str(n+1)

   dz_plot.name = 'member_' + member
   dz_plot.var1_title = 'Member %s Composite Reflectivity (dBZ)' % member
   dz_plot.var2_title = 'Member %s 2-5 km Updraft Helicity (m$^{2}$ s$^{-2}$)' % member
   mem_plot(map, fig, ax1, ax2, ax3, x, y, dz_plot, comp_dz[n,:,:], uh_2to5_convolve[n,:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

