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

#################################### Get input file:  #####################################################

summary_files_temp = os.listdir(summary_dir)
timestep = str(t)
if (len(timestep) == 1):
   timestep = '0' + timestep

for f, file in enumerate(summary_files_temp):
   if ((file[-28:-25] == 'ENV') and (file[-24:-22] == timestep)):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      infile = os.path.join(summary_dir, file)

print 'Matched ENV File: ', timestep
print 'Matched ENV File: ', infile

#################################### Basemap Variables:  #####################################################

resolution      = 'h'
area_thresh     = 1000.

damage_files = '' #['/scratch/skinnerp/2018_newse_post/damage_files/extractDamage_new/extractDamagePolys'] #['/Volumes/fast_scr/pythonletkf/vortex_se/2013-11-17/shapefiles/extractDamage_11-17/extractDamagePaths']

#################################### User-Defined Variables:  #####################################################

domain          = 'full'        #vestigial variable that's still needed in the plotting subroutines ... 
edge            = 7             #number of grid points to remove from near domain boundaries
thin            = 13             #thinning factor for quiver values (e.g. 6 means slice every 6th grid point)

neighborhood    = 15            #grid point radius of prob matched mean neighborhood

plot_alpha      = 0.55          #transparency value for filled contour plots

#################################### Contour Levels:  #####################################################

cape_levels     	= np.arange(250.,4000.,250.)		#(J Kg^-1)
cape_0to3_levels     	= np.arange(25.,400.,25.)		#(J Kg^-1)
cin_levels      	= np.arange(-200.,25.,25.)		#(J Kg^-1)
temp_levels_ugly        = np.arange(50., 110., 5.)
temp_levels             = np.arange(-20., 125., 5.)
temp_levels_850         = np.arange(25., 85., 5.)
temp_levels_700         = np.arange(10., 70., 5.)
the_levels              = np.arange(273., 360., 3.)
td_levels		= np.arange(-16., 88., 4.)		#(deg F)
td_levels_2             = np.arange(12., 80., 4.)             #(deg F) 
td_levels_3             = np.arange(12., 60., 4.)             #(deg F)
td_levels_ugly          = np.arange(32., 80., 4.)             #(deg F) 
#td_levels_ugly          = np.arange(40., 80., 4.)             #(deg F) 
ws_levels_low       	= np.arange(10.,70.,6.)			#(m s^-1)
ws_levels_high       	= np.arange(30.,90.,6.)			#(m s^-1)
ws_levels_low_ms      	= np.arange(5.,35.,3.)			#(m s^-1)
ws_levels_high_ms      	= np.arange(15.,45.,3.)			#(m s^-1)
srh_levels      	= np.arange(40.,640.,40.)		#(m^2 s^-2)
stp_levels      	= np.arange(0.25,7.75,0.5)		#(unitless)
swdown_levels           = np.arange(0.,1300.,100.)              #(W m^-2)
dz_levels_nws           = np.arange(20.0,80.,5.)                #(dBZ)
pbl_levels              = np.arange(0.,2400.,200.)
#mfc_levels              = np.arange(-90.,110.,20.)              #(g/kg*s)
#mfc_levels              = np.arange(-45.,50.,10.)              #(g/kg*s)
mfc_levels              = np.arange(0.,80.,10.)              #(g/kg*s)
corf_levels             = np.arange(10.,60.,5.)         #(m s^-1)
ul_dvg_levels           = np.arange(-6.,7.,1.)          #(10^-5 s^-1)
cp_levels               = np.arange(50,1050,50)
pw_levels               = np.arange(0.0,3.1,0.1)
cwp_levels              = [0.05,0.075,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0,15.0,30.0,60.]
mslp_levels_tc          = np.arange(940.,1016.,4)
mslp_levels             = np.arange(970.,1027.,3)
omega_levels            = np.arange(-9.,11.,2.)
ws_levels_500           = np.arange(10., 90., 8.)

pmm_dz_levels 		= [35., 50.]				#(dBZ) 
pmm_dz_colors_gray	= [cb_colors.gray8, cb_colors.gray8]	#gray contours

pmm_cp_levels           = [40., 1050.]
pmm_cp_colors_gray      = [cb_colors.gray1, cb_colors.gray1]

pmm_cwp_levels          = [0., 60.]
pmm_cwp_colors_gray     = [cb_colors.gray1, cb_colors.gray1]

pmm_pw_levels           = [0., 60.]
pmm_pw_colors_gray      = [cb_colors.gray1, cb_colors.gray1]

#################################### Initialize plot attributes using 'web plot' objects:  #####################################################

ws_plot = web_plot('',                   \
                   '',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   '',            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.orange9,		\
                   'none',		\
                   cb_colors.wind_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)  

cape_plot = web_plot('mlcape',                 \
                   'Ens. Mean 100 hPa MLCAPE (J Kg$^{-1}$)',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   cape_levels,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple3,		\
                   'none',		\
                   cb_colors.cape_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)  

cin_plot = web_plot('mlcin',                  \
                   'Ens. Mean 100 hPa MLCIN (J Kg$^{-1}$)',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   cin_levels,           \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.gray1,		\
                   cb_colors.purple8,		\
                   cb_colors.cin_cmap,             \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)  

pbl_plot = web_plot('pblh',                 \
                   'Ens. Mean PBL Height (m)',                    \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   pbl_levels,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple3,           \
                   'none',              \
                   cb_colors.cape_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)

stp_plot = web_plot('stp',                  \
                   'Ens. Mean Significant Tornado Parameter',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   stp_levels,           \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple3,		\
                   'none',		\
                   cb_colors.cape_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)  

swdown_plot = web_plot('swdown',                  \
                   'Ens. Mean Downward Short Wave Radiation Flux',                   \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   swdown_levels,           \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple3,           \
                   'none',              \
                   cb_colors.wz_cmap_extend,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)

srh_plot = web_plot('',                  \
                   '',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   srh_levels,           \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple3,		\
                   'none',		\
                   cb_colors.cape_cmap,              \
                   'max',                \
                   plot_alpha,               \
                   neighborhood)  

temp_plot = web_plot('temp',                 \
                   'Ens. Mean 2 m Temperature ($^{\circ}$F)',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   temp_levels_ugly,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,		\
                   cb_colors.blue8,		\
                   cb_colors.temp_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)  

temp850_plot = web_plot('temp_850',                 \
                   'Ens. Mean 850 hPa Temperature ($^{\circ}$F)',                       \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   temp_levels_850,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,           \
                   cb_colors.blue8,             \
                   cb_colors.temp_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

temp700_plot = web_plot('temp_700',                 \
                   'Ens. Mean 700 hPa Temperature ($^{\circ}$F)',                       \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   temp_levels_850,          \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,           \
                   cb_colors.blue8,             \
                   cb_colors.temp_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

td_plot = web_plot('td',                   \
                   'Ens. Mean 2 m Dewpoint Temp ($^{\circ}$F)',			\
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   td_levels_ugly,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,		\
                   cb_colors.orange4,		\
                   cb_colors.td_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)  

td850_plot = web_plot('td_850',                   \
                   'Ens. Mean 850 hPa Dewpoint Temp ($^{\circ}$F)',                     \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   td_levels_ugly,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple5,           \
                   cb_colors.orange4,           \
                   cb_colors.td_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

td700_plot = web_plot('td_700',                   \
                   'Ens. Mean 700 hPa Dewpoint Temp ($^{\circ}$F)',                     \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   td_levels_3,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple7,           \
                   cb_colors.orange4,           \
                   cb_colors.td_cmap_ugly,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

mfc_plot = web_plot('pbl_mfc',                  \
                  'Ens. Mean Boundary Layer Moisture Convergence (g kg$^{-1}$ s$^{-1}$)',                       \
                  'Probability Matched Mean - Composite Reflectivity (dBZ)',                    \
                  cb_colors.gray6,      \
                  mfc_levels,           \
                  pmm_dz_levels,                \
                  '',                   \
                  '',                   \
                  pmm_dz_colors_gray,   \
                  cb_colors.red8,            \
                  'none',              \
                  cb_colors.mfc_cmap,   \
                  'max',               \
                  plot_alpha,           \
                  neighborhood)

mslp_plot = web_plot('mslp',                   \
                   'Ens. Mean Mean Sea-Level Pressure (hPa)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   mslp_levels,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.gray1,           \
                   cb_colors.purple8,           \
                   cb_colors.mslp_cmap,              \
                   'both',               \
                   plot_alpha,               \
                   neighborhood)

the_plot = web_plot('the',                   \
                   'Ens. Mean 100 hPa Equivalent Pot. Temp (K)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   the_levels,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.purple2,           \
                   'none',             \
                   cb_colors.temp_cmap,              \
                   'max',               \
                   plot_alpha,               \
                   neighborhood)

dz_plot = web_plot('comp_dz',                   \
                   'Ens. Mean Composite Reflectivity (dBZ)',                  \
                   '',                   \
                   cb_colors.gray6,      \
                   dz_levels_nws,            \
                   [80., 90.],        \
                   '',                   \
                   '',                   \
                   ['none', 'none'],        \
                   cb_colors.c75,               \
                   'none',              \
                   cb_colors.nws_dz_cmap,              \
                   'max',                \
                   0.8,               \
                   neighborhood)

ul_dvg_plot = web_plot('ul_dvg',                \
                  'Ens. Mean 400-250 hPa wind divergence (10$^{-5}$ s$^{-1}$)',         \
                  'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                  cb_colors.gray6,              \
                  ul_dvg_levels,                \
                  pmm_dz_levels,                \
                  '',                           \
                  '',                           \
                  pmm_dz_colors_gray,           \
                  cb_colors.purple9,            \
                  cb_colors.green9,             \
                  cb_colors.ul_dvg_cmap,        \
                  'both',                       \
                  plot_alpha,                   \
                  neighborhood)

corf_plot = web_plot('corf',                  \
                  'Ens. Mean Corfidi Storm Motion (kts)',                       \
                  'Probability Matched Mean - Composite Reflectivity (dBZ)',                    \
                  cb_colors.gray6,      \
                  corf_levels,           \
                  pmm_dz_levels,                \
                  '',                   \
                  '',                   \
                  pmm_dz_colors_gray,   \
                  cb_colors.orange9,            \
                  cb_colors.purple9,             \
                  cb_colors.corf_cmap,   \
                  'both',               \
                  plot_alpha,           \
                  neighborhood)

pw_plot = web_plot('',                  \
                   'Ens. Mean Precipitable Water (in)',         \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,     \
                   pw_levels,           \
                   pmm_dz_levels,       \
                   '',                  \
                   '',                  \
                   pmm_dz_colors_gray,  \
                   cb_colors.red9,              \
                   'none',              \
                   cb_colors.pw_cmap,   \
                   'max',               \
                   plot_alpha,          \
                   neighborhood)

ctp_plot = web_plot('',                 \
                    'Ens. Mean Cloud Top Pressure (hPa)',               \
                    '',                 \
                    cb_colors.gray6,    \
                    cp_levels,          \
                    pmm_cp_levels,      \
                    '',                 \
                    '',                 \
                    pmm_cp_colors_gray, \
                    'none',             \
                    'none',             \
                    cb_colors.cp_cmap,  \
                    'neither',          \
                    plot_alpha,         \
                    neighborhood)

cwp_plot = web_plot('',                 \
                    'Ens. Mean Cloud Water Path (kg m$^{-2}$)',         \
                    '',                 \
                    cb_colors.gray6,    \
                    cwp_levels,         \
                    pmm_cwp_levels,     \
                    '',                 \
                    '',                 \
                    pmm_cwp_colors_gray,        \
                    cb_colors.purple9,          \
                    'none',             \
                    cb_colors.cwp_cmap, \
                    'max',              \
                    plot_alpha,         \
                    neighborhood)

omega_plot = web_plot('w_700p',                   \
                   'Ens. Mean 700 hPa Omega (Pa s$^{-1}$)',                 \
                   'Probability Matched Mean - Composite Reflectivity (dBZ)',                   \
                   cb_colors.gray6,      \
                   omega_levels,            \
                   pmm_dz_levels,        \
                   '',                   \
                   '',                   \
                   pmm_dz_colors_gray,        \
                   cb_colors.green9,           \
                   cb_colors.purple9,             \
                   cb_colors.omega_cmap,              \
                   'both',               \
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

xlat = fin.variables["xlat"][edge:-edge,edge:-edge]                     #latitude (dec deg; Lambert conformal)
xlon = fin.variables["xlon"][edge:-edge,edge:-edge]                     #longitude (dec deg; Lambert conformal)

#xlat = xlat[edge:-edge,edge:-edge]
#xlon = xlon[edge:-edge,edge:-edge]

sw_lat_full = xlat[0,0]
sw_lon_full = xlon[0,0]
ne_lat_full = xlat[-1,-1]
ne_lon_full = xlon[-1,-1]

######################### Read PMM composite reflectivity: ####################################

pmm_dz = fin.variables["comp_dz_pmm"][edge:-edge,edge:-edge]

######################### Initialize variables: ####################################

mslp = fin.variables["mslp"][:,edge:-edge,edge:-edge]
tf_2 = fin.variables["t_2"][:,edge:-edge,edge:-edge]
tdf_2 = fin.variables["td_2"][:,edge:-edge,edge:-edge]
tf_850p = fin.variables["t_850p"][:,edge:-edge,edge:-edge]
tdf_850p = fin.variables["td_850p"][:,edge:-edge,edge:-edge]
tf_700p = fin.variables["t_700p"][:,edge:-edge,edge:-edge]
tdf_700p = fin.variables["td_700p"][:,edge:-edge,edge:-edge]
the_ml = fin.variables["th_e_ml"][:,edge:-edge,edge:-edge]
u_10 = fin.variables["u_10"][:,edge:-edge,edge:-edge]
v_10 = fin.variables["v_10"][:,edge:-edge,edge:-edge]
u_500 = fin.variables["u_500"][:,edge:-edge,edge:-edge]
v_500 = fin.variables["v_500"][:,edge:-edge,edge:-edge]
u_850p = fin.variables["u_850p"][:,edge:-edge,edge:-edge]
v_850p = fin.variables["v_850p"][:,edge:-edge,edge:-edge]
u_700p = fin.variables["u_700p"][:,edge:-edge,edge:-edge]
v_700p = fin.variables["v_700p"][:,edge:-edge,edge:-edge]
u_500p = fin.variables["u_500p"][:,edge:-edge,edge:-edge]
v_500p = fin.variables["v_500p"][:,edge:-edge,edge:-edge]
u_bl = fin.variables["u_bl"][:,edge:-edge,edge:-edge]
v_bl = fin.variables["v_bl"][:,edge:-edge,edge:-edge]
u_dvg = fin.variables["u_dvg"][:,edge:-edge,edge:-edge]
v_dvg = fin.variables["v_dvg"][:,edge:-edge,edge:-edge]
u_850_300p = fin.variables["u_850_300p"][:,edge:-edge,edge:-edge]
v_850_300p = fin.variables["v_850_300p"][:,edge:-edge,edge:-edge]
cape = fin.variables["cape_ml"][:,edge:-edge,edge:-edge]
mucape = fin.variables["cape_mu"][:,edge:-edge,edge:-edge]
sbcape = fin.variables["cape_sfc"][:,edge:-edge,edge:-edge]
#cape_0to3 = fin.variables["cape_0to3_ml"][:,edge:-edge,edge:-edge]
#pbl_hgt = fin.variables["pbl_hgt"][:,edge:-edge,edge:-edge]   #not in WRFWOF yet
pbl_mfc = fin.variables["pbl_mfc"][:,edge:-edge,edge:-edge]
cin = fin.variables["cin_ml"][:,edge:-edge,edge:-edge]
stp = fin.variables["stp_ml"][:,edge:-edge,edge:-edge]
corf_us_u = fin.variables["corf_us_u"][:,edge:-edge,edge:-edge]
corf_us_v = fin.variables["corf_us_v"][:,edge:-edge,edge:-edge]
corf_ds_u = fin.variables["corf_ds_u"][:,edge:-edge,edge:-edge]
corf_ds_v = fin.variables["corf_ds_v"][:,edge:-edge,edge:-edge]
bunk_r_u = fin.variables["bunk_r_u"][:,edge:-edge,edge:-edge]
bunk_r_v = fin.variables["bunk_r_v"][:,edge:-edge,edge:-edge]
srh_0to1 = fin.variables["srh_0to1"][:,edge:-edge,edge:-edge]
srh_0to3 = fin.variables["srh_0to3"][:,edge:-edge,edge:-edge]
shear_u_0to6 = fin.variables["shear_u_0to6"][:,edge:-edge,edge:-edge]
shear_v_0to6 = fin.variables["shear_v_0to6"][:,edge:-edge,edge:-edge]
shear_u_0to1 = fin.variables["shear_u_0to1"][:,edge:-edge,edge:-edge]
shear_v_0to1 = fin.variables["shear_v_0to1"][:,edge:-edge,edge:-edge]
swdown = fin.variables["sw_down"][:,edge:-edge,edge:-edge]
comp_dz = fin.variables["comp_dz"][:,edge:-edge,edge:-edge]
ctp = fin.variables["ctp"][:,edge:-edge,edge:-edge]
pwat = fin.variables["pw"][:,edge:-edge,edge:-edge]
iwp = fin.variables["iwp"][:,edge:-edge,edge:-edge]
lwp = fin.variables["lwp"][:,edge:-edge,edge:-edge]
w_700p = fin.variables["w_700p"][:,edge:-edge,edge:-edge]

### mask areas with no CAPE in cin plot ###: 
masked_cin = np.ma.masked_where((cape <= 0.), (cin))

mean_mslp = np.mean(mslp, axis=0)
mean_tf_2 = np.mean(tf_2, axis=0)
mean_tdf_2 = np.mean(tdf_2, axis=0)
mean_tf_850 = np.mean(tf_850p, axis=0)
mean_tdf_850 = np.mean(tdf_850p, axis=0)
mean_tf_700 = np.mean(tf_700p, axis=0)
mean_tdf_700 = np.mean(tdf_700p, axis=0)
mean_the_ml = np.mean(the_ml, axis=0)
mean_u_10 = np.mean(u_10, axis=0)
mean_v_10 = np.mean(v_10, axis=0)
mean_u_500 = np.mean(u_500, axis=0)
mean_v_500 = np.mean(v_500, axis=0)
mean_u_850p = np.mean(u_850p, axis=0)
mean_v_850p = np.mean(v_850p, axis=0)
mean_u_700p = np.mean(u_700p, axis=0)
mean_v_700p = np.mean(v_700p, axis=0)
mean_u_500p = np.mean(u_500p, axis=0)
mean_v_500p = np.mean(v_500p, axis=0)
mean_u_850_300p = np.mean(u_850_300p, axis=0)
mean_v_850_300p = np.mean(v_850_300p, axis=0)
mean_u_bl = np.mean(u_bl, axis=0)
mean_v_bl = np.mean(v_bl, axis=0)
mean_u_dvg = np.mean(u_dvg, axis=0)
mean_v_dvg = np.mean(v_dvg, axis=0)
mean_cape = np.mean(cape, axis=0)
mean_mucape = np.mean(mucape, axis=0)
mean_sbcape = np.mean(sbcape, axis=0)
#mean_cape_0to3 = np.mean(cape_0to3, axis=0)
#mean_pblh = np.mean(pbl_hgt, axis=0)   #not in WRFWOF yet
mean_pbl_mfc = np.mean(pbl_mfc, axis=0)
mean_cin = np.mean(masked_cin, axis=0)
mean_stp = np.mean(stp, axis=0)
mean_bunk_r_u = np.mean(bunk_r_u, axis=0)
mean_bunk_r_v = np.mean(bunk_r_v, axis=0)
mean_corf_us_u = np.mean(corf_us_u, axis=0)
mean_corf_us_v = np.mean(corf_us_v, axis=0)
mean_corf_ds_u = np.mean(corf_ds_u, axis=0)
mean_corf_ds_v = np.mean(corf_ds_v, axis=0)
mean_srh_0to1 = np.mean(srh_0to1, axis=0)
mean_srh_0to3 = np.mean(srh_0to3, axis=0)
mean_shear_u_0to6 = np.mean(shear_u_0to6, axis=0)
mean_shear_v_0to6 = np.mean(shear_v_0to6, axis=0)
mean_shear_u_0to1 = np.mean(shear_u_0to1, axis=0)
mean_shear_v_0to1 = np.mean(shear_v_0to1, axis=0)
mean_swdown = np.mean(swdown, axis=0)
mean_comp_dz = np.mean(comp_dz, axis=0)
mean_pwat = np.mean(pwat, axis=0)
mean_ctp = np.mean(ctp, axis=0)
mean_cwp = np.mean((iwp+lwp), axis=0) / 1000.   #convert to kg/m^2
mean_omega = np.mean(w_700p, axis=0)

fin.close()
del fin

##################### Calculate ensemble mean values: ###################################

mean_ws_500 = np.sqrt(mean_u_500**2 + mean_v_500**2)
mean_ws_850p = np.sqrt(mean_u_850p**2 + mean_v_850p**2)
mean_ws_700p = np.sqrt(mean_u_700p**2 + mean_v_700p**2)
mean_ws_500p = np.sqrt(mean_u_500p**2 + mean_v_500p**2)
mean_ws_850_300p = np.sqrt(mean_u_850_300p**2 + mean_v_850_300p**2)
mean_ws_bl = np.sqrt(mean_u_bl**2 + mean_v_bl**2)
mean_shear_0to6 = np.sqrt(mean_shear_u_0to6**2 + mean_shear_v_0to6**2)
mean_shear_0to1 = np.sqrt(mean_shear_u_0to1**2 + mean_shear_v_0to1**2)
mean_bunk = np.sqrt(mean_bunk_r_u**2 + mean_bunk_r_v**2)
mean_corf_us = np.sqrt(mean_corf_us_u**2 + mean_corf_us_v**2)
mean_corf_ds = np.sqrt(mean_corf_ds_u**2 + mean_corf_ds_v**2)
mean_ul_dvg = np.sqrt(mean_u_dvg**2 + mean_v_dvg**2)

#################### Thin arrays used for quiver plots (removes extra vectors): #################################

quiv_xlon = xlon[0:-1:thin,0:-1:thin]
quiv_xlat = xlat[0:-1:thin,0:-1:thin]
quiv_u_10 = mean_u_10[0:-1:thin,0:-1:thin]
quiv_v_10 = mean_v_10[0:-1:thin,0:-1:thin]
quiv_u_500 = mean_u_500[0:-1:thin,0:-1:thin]
quiv_v_500 = mean_v_500[0:-1:thin,0:-1:thin]
quiv_u_850p = mean_u_850p[0:-1:thin,0:-1:thin]
quiv_v_850p = mean_v_850p[0:-1:thin,0:-1:thin]
quiv_u_700p = mean_u_700p[0:-1:thin,0:-1:thin]
quiv_v_700p = mean_v_700p[0:-1:thin,0:-1:thin]
quiv_u_500p = mean_u_500p[0:-1:thin,0:-1:thin]
quiv_v_500p = mean_v_500p[0:-1:thin,0:-1:thin]
quiv_u_850_300p = mean_u_850_300p[0:-1:thin,0:-1:thin]
quiv_v_850_300p = mean_v_850_300p[0:-1:thin,0:-1:thin]
quiv_u_bl = mean_u_bl[0:-1:thin,0:-1:thin]
quiv_v_bl = mean_v_bl[0:-1:thin,0:-1:thin]
quiv_u_dvg = mean_u_dvg[0:-1:thin,0:-1:thin]
quiv_v_dvg = mean_v_dvg[0:-1:thin,0:-1:thin]
quiv_shear0to6_u = mean_shear_u_0to6[0:-1:thin,0:-1:thin]
quiv_shear0to6_v = mean_shear_v_0to6[0:-1:thin,0:-1:thin]
quiv_shear0to1_u = mean_shear_u_0to1[0:-1:thin,0:-1:thin]
quiv_shear0to1_v = mean_shear_v_0to1[0:-1:thin,0:-1:thin]
quiv_bunk_u = mean_bunk_r_u[0:-1:thin,0:-1:thin]
quiv_bunk_v = mean_bunk_r_v[0:-1:thin,0:-1:thin]
quiv_corf_us_u = mean_corf_us_u[0:-1:thin,0:-1:thin]
quiv_corf_us_v = mean_corf_us_v[0:-1:thin,0:-1:thin]
quiv_corf_ds_u = mean_corf_ds_u[0:-1:thin,0:-1:thin]
quiv_corf_ds_v = mean_corf_ds_v[0:-1:thin,0:-1:thin]

################################# Make Figure Template: ###################################################

print 'basemap part'

#Load pickled basemap instance for faster plotting: 

fig, ax1, ax2, ax3 = create_fig_nomap()

map_temp = pl.load(open(mapname, 'rb'))

P.sca(ax1)
map = mymap_boundaries(map_temp, damage_files)

map.drawcounties(linewidth=0.7, color=cb_colors.gray4)
map.drawstates(linewidth=1., color=cb_colors.gray6)
map.drawcoastlines(linewidth=1., color=cb_colors.gray5)
map.drawcountries(linewidth=1., color=cb_colors.gray5)

#map, fig, ax1, ax2, ax3 = create_fig(sw_lat_full, sw_lon_full, ne_lat_full, ne_lon_full, true_lat_1, true_lat_2, cen_lat, stand_lon, damage_files, resolution, area_thresh)

x, y = map(xlon[:], xlat[:])

##########################################################################################################
###################################### Make Plots: #######################################################
##########################################################################################################

print 'plot part'

################# Rotate quiver values according to map projection: #####################

#q_u, q_v = map.rotate_vector(quiv_u_10[:,:], quiv_v_10[:,:], quiv_xlon, quiv_xlat, returnxy=False)
#q_u_500, q_v_500 = map.rotate_vector(quiv_u_500[:,:], quiv_v_500[:,:], quiv_xlon, quiv_xlat, returnxy=False)
#shear0to6_q_u, shear0to6_q_v = map.rotate_vector(quiv_shear0to6_u[:,:], quiv_shear0to6_v[:,:], quiv_xlon, quiv_xlat, returnxy=False)
#shear0to1_q_u, shear0to1_q_v = map.rotate_vector(quiv_shear0to1_u[:,:], quiv_shear0to1_v[:,:], quiv_xlon, quiv_xlat, returnxy=False)
#bunk_q_u, bunk_q_v = map.rotate_vector(quiv_bunk_u[:,:], quiv_bunk_v[:,:], quiv_xlon, quiv_xlat, returnxy=False)

q_u = quiv_u_10
q_v = quiv_v_10
q_u_500 = quiv_u_500
q_v_500 = quiv_v_500
q_u_bl = quiv_u_bl
q_v_bl = quiv_v_bl
q_u_850p = quiv_u_850p
q_v_850p = quiv_v_850p
q_u_700p = quiv_u_700p
q_v_700p = quiv_v_700p
q_u_500p = quiv_u_500p
q_v_500p = quiv_v_500p
q_u_850_300p = quiv_u_850_300p
q_v_850_300p = quiv_v_850_300p
q_u_dvg = quiv_u_dvg
q_v_dvg = quiv_v_dvg
shear0to1_q_u = quiv_shear0to1_u
shear0to1_q_v = quiv_shear0to1_v
shear0to6_q_u = quiv_shear0to6_u
shear0to6_q_v = quiv_shear0to6_v
bunk_q_u = quiv_bunk_u
bunk_q_v = quiv_bunk_v
corf_us_q_u = quiv_corf_us_u
corf_us_q_v = quiv_corf_us_v
corf_ds_q_u = quiv_corf_ds_u
corf_ds_q_v = quiv_corf_ds_v

######################## Environment Plots: #####################

env_plot(map, fig, ax1, ax2, ax3, x, y, temp_plot, mean_tf_2[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True') 

env_plot(map, fig, ax1, ax2, ax3, x, y, td_plot, mean_tdf_2[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True') 

env_plot(map, fig, ax1, ax2, ax3, x, y, mslp_plot, mean_mslp[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True')

env_plot(map, fig, ax1, ax2, ax3, x, y, temp850_plot, mean_tf_850[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_850p, q_v_850p, 500, 5, 0, spec='False', quiv='True')

env_plot(map, fig, ax1, ax2, ax3, x, y, td850_plot, mean_tdf_850[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_850p, q_v_850p, 500, 5, 0, spec='False', quiv='True')

env_plot(map, fig, ax1, ax2, ax3, x, y, temp700_plot, mean_tf_700[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_700p, q_v_700p, 500, 5, 0, spec='False', quiv='True')

env_plot(map, fig, ax1, ax2, ax3, x, y, td700_plot, mean_tdf_700[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_700p, q_v_700p, 500, 5, 0, spec='False', quiv='True')

env_plot(map, fig, ax1, ax2, ax3, x, y, the_plot, mean_the_ml[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u, q_v, 500, 5, 0, spec='False', quiv='True') 

env_plot(map, fig, ax1, ax2, ax3, x, y, cape_plot, mean_cape[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

cape_plot.name = 'mucape'
cape_plot.var1_title = 'Ens. Mean MUCAPE (J Kg$^{-1}$)' 
env_plot(map, fig, ax1, ax2, ax3, x, y, cape_plot, mean_mucape[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

cape_plot.name = 'sbcape'
cape_plot.var1_title = 'Ens. Mean SBCAPE (J Kg$^{-1}$)' 
env_plot(map, fig, ax1, ax2, ax3, x, y, cape_plot, mean_sbcape[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

#env_plot(map, fig, ax1, ax2, ax3, x, y, pbl_plot, mean_pblh[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')  #no WRFWOF yet

#cape_plot.name = 'mlcape_0to3'
#cape_plot.var1_title = 'Ens. Mean 100 hPa MLCAPE below 3 km (J Kg$^{-1}$)' 
#cape_plot.var1_levels = cape_0to3_levels
#env_plot(map, fig, ax1, ax2, ax3, x, y, cape_plot, mean_cape_0to3[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

env_plot(map, fig, ax1, ax2, ax3, x, y, cin_plot, mean_cin[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

env_plot(map, fig, ax1, ax2, ax3, x, y, stp_plot, mean_stp[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

env_plot(map, fig, ax1, ax2, ax3, x, y, swdown_plot, mean_swdown[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

env_plot(map, fig, ax1, ax2, ax3, x, y, dz_plot, mean_comp_dz[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

srh_plot.name = 'srh0to1'
srh_plot.var1_title = 'Ens. Mean 0 - 1 km Storm Relative Helicity (m$^{2}$ s$^{-2}$)'
env_plot(map, fig, ax1, ax2, ax3, x, y, srh_plot, mean_srh_0to1[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

srh_plot.name = 'srh0to3'
srh_plot.var1_title = 'Ens. Mean 0 - 3 km Storm Relative Helicity (m$^{2}$ s$^{-2}$)'
env_plot(map, fig, ax1, ax2, ax3, x, y, srh_plot, mean_srh_0to3[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False') 

mfc_plot.name = 'pbl_mfc'
mfc_plot.var1_title = 'Ens. Mean Lowest 75 hPa Max Moisture Convergence (g Kg$^{-1}$ s$^{-1}$)'
env_plot(map, fig, ax1, ax2, ax3, x, y, mfc_plot, mean_pbl_mfc[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_bl, q_v_bl, 500, 5, 0, spec='False', quiv='True')

ws_plot.name = 'ws_500'
ws_plot.var1_title = 'Ens. Mean 500 m Wind Speed (kts)'
ws_plot.var1_levels = ws_levels_low 
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_ws_500[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_500, q_v_500, 1000, 5, 0, spec='False', quiv='True') 

ws_plot.name = 'ws_850hpa'
ws_plot.var1_title = 'Ens. Mean 850 hPa Wind Speed (kts)'
ws_plot.var1_levels = ws_levels_low
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_ws_850p[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_850p, q_v_850p, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'ws_700hpa'
ws_plot.var1_title = 'Ens. Mean 700 hPa Wind Speed (kts)'
ws_plot.var1_levels = ws_levels_low
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_ws_700p[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_700p, q_v_700p, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'ws_500hpa'
ws_plot.var1_title = 'Ens. Mean 500 hPa Wind Speed (kts)'
ws_plot.var1_levels = ws_levels_500
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_ws_500p[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_500p, q_v_500p, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'ws_850_300_hpa'
ws_plot.var1_title = 'Ens. Mean 850-300 hPa Mean Wind Speed (kts)'
ws_plot.var1_levels = ws_levels_low
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_ws_850_300p[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, q_u_850_300p, q_v_850_300p, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'shear0to6'
ws_plot.var1_title = 'Ens. Mean 0 - 6 km Shear (kts)'
ws_plot.var1_levels = ws_levels_high 
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_shear_0to6[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, shear0to6_q_u, shear0to6_q_v, 1000, 5, 0, spec='False', quiv='True') 

ws_plot.name = 'shear0to1'
ws_plot.var1_title = 'Ens. Mean 0 - 1 km Shear (kts)'
ws_plot.var1_levels = ws_levels_low 
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_shear_0to1[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, shear0to1_q_u, shear0to1_q_v, 1000, 5, 0, spec='False', quiv='True') 

ws_plot.name = 'bunk'
ws_plot.var1_title = 'Ens. Mean Bunkers Storm Motion (kts)'
ws_plot.var1_levels = ws_levels_low 
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_bunk[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, bunk_q_u, bunk_q_v, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'corf_us'
ws_plot.var1_title = 'Ens. Mean Corfidi Upshear Motion (kts)'
ws_plot.var1_levels = corf_levels
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_corf_us[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, corf_us_q_u, corf_us_q_v, 1000, 5, 0, spec='False', quiv='True')

ws_plot.name = 'corf_ds'
ws_plot.var1_title = 'Ens. Mean Corfidi Downshear Motion (kts)'
ws_plot.var1_levels = ws_levels_high
env_plot(map, fig, ax1, ax2, ax3, x, y, ws_plot, mean_corf_ds[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, corf_ds_q_u, corf_ds_q_v, 1000, 5, 0, spec='False', quiv='True')

ul_dvg_plot.name = 'ul_dvg'
ul_dvg_plot_var1_title = 'Ens. Mean 400-250 hPa Divergence (s$^{-1}$)'
ul_dvg_plot_var1_levels = ul_dvg_levels
env_plot(map, fig, ax1, ax2, ax3, x, y, ul_dvg_plot, mean_ul_dvg[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

cwp_plot.name = 'cwp'
cwp_plot.var1_title = 'Ens. Mean Cloud Water Path (kg m$^{-2}$)'
sat_plot(map, fig, ax1, ax2, ax3, x, y, cwp_plot, mean_cwp[:,:], mean_cwp[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

pw_plot.name = 'pw'
pw_plot.var1_title = 'Ens. Mean Precipitable Water (in)'
env_plot(map, fig, ax1, ax2, ax3, x, y, pw_plot, mean_pwat[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

ctp_plot.name = 'ctp'
ctp_plot.var1_title = 'Ens. Mean Cloud Top Pressure (hPa)'
sat_plot(map, fig, ax1, ax2, ax3, x, y, ctp_plot, mean_ctp[:,:], mean_ctp[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')

omega_plot.name = 'omega'
omega_plot.var1_title = 'Ens. Mean 700 hPa Omega (Pa s$^{-1}$)'
omega_plot.var1_levels = omega_levels
env_plot(map, fig, ax1, ax2, ax3, x, y, omega_plot, mean_omega[:,:], pmm_dz[:,:], t, init_label, valid_label, domain, outdir, '', '', '', 5, 0, spec='False', quiv='False')
 
