#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
import skimage
from skimage.morphology import label
from skimage.measure import regionprops
import math
from math import radians, tan, sin, cos, pi, atan, sqrt, pow, asin, acos
import pylab as P
import numpy as np
from numpy import NAN
import sys
import netCDF4
from optparse import OptionParser
from netcdftime import utime
import os
import time as timeit
from optparse import OptionParser
import pickle as pl
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input directory of summary files to plot")
parser.add_option("-a", dest="date_dir", type="string", default= None, help="Input directory of forecasts for case")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for images)")
#parser.add_option("-m", dest="mapname", type="string", help = "Path to Pickled Basemap instance")
parser.add_option("-t", dest="t", type="int", help = "Timestep to process")
parser.add_option("-n", dest="nt", type="int", help = "Number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.date_dir == None) or (options.outdir == None) or (options.t == None) or (options.nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   date_dir = options.date_dir
   outdir = options.outdir
   t = options.t
   nt = options.nt

fcst_times = ['1900', '2000', '2030', '2100', '2130', '2200', '2230', '2300', '2330', '0000', '0030', '0100', '0130', '0200', '0230', '0300']
fcst_nt = [72, 60, 36, 48, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 35]


#path to pickled map instance for plotting: 
#mapname = '/scratch2/patrick.skinner/images/map.pickle'

#################################### User-Defined Variables:  #####################################################

newse_dz_thresh            = 45.002             

########## 99.95th percentile thresholds for NEWS-e 2017 cases: 

uh_0to2_thresh             = 14.217
uh_2to5_thresh             = 65.790

area_thresh_uh             = 10.            #Minimum area of rotation object
area_thresh_dz             = 12.           #Minimum area of reflectivity object
cont_thresh                = 2             #Must be greater than this number of forecast timesteps in a rotation object to be retained

swath_window               = 3             #+- 3 timesteps for swaths
domain                     = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge                       = 7          #number of grid points to remove from near domain boundaries

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

plot_alpha                 = 0.55                #transparency value for filled contour plots

#################################### Basemap Variables:  #####################################################

resolution      = 'h'
area_thresh     = 1000.

damage_files = '' #['/scratch/skinnerp/2018_newse_post/damage_files/extractDamage_new/extractDamagePolys'] #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

uh0to2_paint_plot = v_plot('uh0to2_paint_tl',                   \
                   'NEWS-e Member 0-2 km Updraft Helicity Objects',                        \
                   'Time-lagged Ensemble Paintball',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [900., 1000.],              \
                   [0.01, 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.tl_paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)

uh2to5_paint_plot = v_plot('uh2to5_paint_tl',                   \
                   'NEWS-e Member 2-5 km Updraft Helicity Objects',                        \
                   'Time-lagged Ensemble Paintball',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [900., 1000.],              \
                   [0.01, 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.tl_paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)

dz_paint_plot = v_plot('dz_paint_tl',                   \
                   'NEWS-e Member 45-dBZ Composite Reflectivity',                        \
                   'Time-lagged Ensemble Paintball',                  \
                   cb_colors.gray6,   \
                   cb_colors.gray6,     \
                   [900., 1000.],              \
                   [39., 1000.],     \
                   [cb_colors.gray8],                      \
                   cb_colors.tl_paintball_colors_list,                     \
                   '',  \
                   '',  \
                   'max',  \
                   0.6, \
                   neighborhood)

############################ Find WRFOUT files to process: #################################

### Find if enough forecast timesteps are available to create a rotation track ### 

if ((t < swath_window) or (t > (nt - (swath_window+1)))): 
   blank = 'True'
   print t, blank
else: 
   blank = 'False'

### Find ENS Summary files ### 

ne = 18

ens_files = []
summary_files_temp = os.listdir(summary_dir)

for f, file in enumerate(summary_files_temp):
   if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file)

ens_files.sort()

ens_file = ens_files[t]

infile = os.path.join(summary_dir, ens_file)

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

xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

fin.close()
del fin

######################### Initialize Object Variables: ####################################

uh_2to5_obj = np.zeros((len(fcst_times), ne, xlat.shape[0], xlat.shape[1]))
uh_0to2_obj = np.zeros((len(fcst_times), ne, xlat.shape[0], xlat.shape[1]))
comp_dz_obj = np.zeros((len(fcst_times), ne, xlat.shape[0], xlat.shape[1]))

######################### Find object from all forecasts valid at current time: ####################################

for d, fcst in enumerate(fcst_times): 

   temp_init_time = int(init_hour)
   if (temp_init_time < 12): 
      temp_init_time = temp_init_time + 24

   temp_fcst = int(fcst[0:2])
   if (temp_fcst < 12): 
      temp_fcst = temp_fcst + 24

   temp_init_min = int(init_min)
   temp_fcst_min = int(fcst[2:4])
      
   temp_init_time_seconds = temp_init_time * 3600. + temp_init_min * 60. 
   temp_fcst_time_seconds = temp_fcst * 3600. + temp_fcst_min * 60. 

   if (temp_fcst_time_seconds <= temp_init_time_seconds): 

### Find ENS Summary files ### 

      temp_nt = fcst_nt[d]
      ens_files = []
      temp_dir = os.path.join(date_dir, fcst)
      summary_files_temp = os.listdir(temp_dir)

      for f, file in enumerate(summary_files_temp):
         if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
            ens_files.append(file)

      ens_files.sort()

### Look for ENS file valid at current timei ### 

      for f, file in enumerate(ens_files): 
         temp_hour = file[-7:-5]
         temp_min = file[-5:-3]

         if ((temp_hour == valid_hour) and (temp_min == valid_min)):  
            temp_t = int(file[-24:-22])
 
            if ((temp_t < swath_window) or (temp_t > (temp_nt - (swath_window+1)))):
               temp_blank = 'True'
            else:
               temp_blank = 'False'

            ens_file = file

            infile = os.path.join(temp_dir, ens_file)
   
            try:                                                 #open WRFOUT file
               fin = netCDF4.Dataset(infile, "r")
               print "Opening %s \n" % infile
            except:
               print "%s does not exist! \n" %infile
               sys.exit(1)

########## Read comp_dz variable: #############

            comp_dz                            = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]
 
            fin.close()
            del fin

            if (temp_blank == 'False'): 
               for tt in range(temp_t-swath_window, temp_t+swath_window+1):
                  ens_file = ens_files[tt]

                  infile = os.path.join(temp_dir, ens_file)

                  try:                                                 #open WRFOUT file
                     fin = netCDF4.Dataset(infile, "r")
                     print "Opening %s \n" % infile
                  except:
                     print "%s does not exist! \n" %infile
                     sys.exit(1)

                  if (tt == (temp_t-swath_window)):
                     uh_2to5                            = fin.variables["uh_2to5"][:,edge:-edge,edge:-edge]
                     uh_0to2                            = fin.variables["uh_0to2"][:,edge:-edge,edge:-edge]

                     uh_2to5_indices                    = uh_2to5 * 0. + tt
                     uh_0to2_indices                    = uh_0to2 * 0. + tt

                  else:
                     uh_2to5_indices                    = np.where((fin.variables["uh_2to5"][:,edge:-edge,edge:-edge] > uh_2to5), tt, uh_2to5_indices)
                     uh_0to2_indices                    = np.where((fin.variables["uh_0to2"][:,edge:-edge,edge:-edge] > uh_0to2), tt, uh_0to2_indices)

                     uh_2to5                            = np.where((fin.variables["uh_2to5"][:,edge:-edge,edge:-edge] > uh_2to5), fin.variables["uh_2to5"][:,edge:-edge,edge:-edge], uh_2to5)
                     uh_0to2                            = np.where((fin.variables["uh_0to2"][:,edge:-edge,edge:-edge] > uh_0to2), fin.variables["uh_0to2"][:,edge:-edge,edge:-edge], uh_0to2)

######## Close file: #######

                  fin.close()
                  del fin

###

            radmask = comp_dz[0,:,:] * 0. #dummy variable (will be used when MRMS obs are included) 
 
            if (temp_blank == 'False'): 
               uh_2to5_obj[d,:,:,:] = find_objects_swath(uh_2to5, uh_2to5_indices, radmask, uh_2to5_thresh, area_thresh_uh, cont_thresh)
               uh_0to2_obj[d,:,:,:] = find_objects_swath(uh_0to2, uh_0to2_indices, radmask, uh_0to2_thresh, area_thresh_uh, cont_thresh)

            comp_dz_obj[d,:,:,:] = find_objects_timestep(comp_dz, radmask, newse_dz_thresh, area_thresh_dz)

################################# Make Figure Template: ###################################################

radmask = comp_dz[0,:,:] * 0. #dummy variable (will be used when MRMS obs are included) 
print 'basemap part'

#Load pickled basemap instance for faster plotting: 

map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])
xx, yy = map.makegrid(xlat.shape[1], xlat.shape[0], returnxy=True)[2:4]   #equidistant x/y grid for streamline plots

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

aws_qc = xlat * 0.   ### set aws_qc to 0 since not plotting verification 

### Plot dz objects:
print 't is:  ', t 
painttl_plot(map, fig, ax1, ax2, ax3, x, y, x, y, dz_paint_plot, aws_qc, comp_dz_obj, radmask, t, init_label, valid_label, cb_colors.tl_paintball_colors, fcst_times, domain, outdir, 5, 0, blank='False')

###########################################################################################################
### make new figure for each paint plot ... can't figure out how to remove each members plot from basemap
###########################################################################################################

print 'basemap part, part 2 - the basemappening'

map2, fig2, ax21, ax22, ax23 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

### Plot uh2to5 objects: 

painttl_plot(map2, fig2, ax21, ax22, ax23, x, y, x, y, uh2to5_paint_plot, aws_qc, uh_2to5_obj, radmask, t, init_label, valid_label, cb_colors.tl_paintball_colors, fcst_times, domain, outdir, 5, 0, blank)

print 'basemap part, part 3 - So very tired ...'

map3, fig3, ax31, ax32, ax33 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

### Plot uh0to2 objects: 

painttl_plot(map3, fig3, ax31, ax32, ax33, x, y, x, y, uh0to2_paint_plot, aws_qc, uh_0to2_obj, radmask, t, init_label, valid_label, cb_colors.tl_paintball_colors, fcst_times,  domain, outdir, 5, 0, blank)


