###################################################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
from scipy import spatial
import pylab as P
import numpy as np
import sys
import os
import time
import netCDF4
from optparse import OptionParser
from news_e_post_cbook import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Directory of summary files (both input and output dir)")
parser.add_option("-t", dest="t", type="int", help = "Forecast timestep being processed")
parser.add_option("-n", dest="fcst_nt", type="int", help = "Total number of forecast timesteps")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.t == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   t = options.t
   fcst_nt = options.fcst_nt

print 'SWATH RAINFALL STARTED: ', t

############################ Set Thresholds for convolution and neighborhood probabilities: #################################

radius_max_9km             = 3			#grid point radius for maximum value filter (3x3 square neighborhood)
radius_max_15km            = 5			#grid point radius for maximum value filter (5x5 square neighborhood)
radius_max_27km            = 9			#grid point radius for maximum value filter (9x9 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean

comp_dz_thresh		   = 45.
############################ Set Thresholds for convolution and neighborhood probabilities: #################################
rain_1_thresh        = 0.5		#0.5 inches
rain_1_3hr_thresh	   = 1.			#1 inch, 3 hr
rain_1_6hr_thresh	   = 2			#2 inches, 6 hr

rain_2_thresh        = 1.			#1 inch
rain_2_3hr_thresh	   = 2.			#2 inches, 3 hr
rain_2_6hr_thresh	   = 3.			#3 inches, 6 hr

rain_3_thresh        = 2.			#2 inches
rain_3_3hr_thresh	   = 3.			#3 inches, 3 hr
rain_3_6hr_thresh	   = 5.			#5 inches, 6 hr

rain_4_thresh        = 3.        #3 inches


rain_5_thresh        = 5.        #5 inches
############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ###

ne = 40

ens_files = []
summary_files_temp = os.listdir(summary_dir)

for f, file in enumerate(summary_files_temp):
   if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file)

ens_files.sort()

for tt in range(0, t+1):
   ens_file = ens_files[tt]
   infile = os.path.join(summary_dir, ens_file)

   if (tt == t):
   ### Set output path ###
      outname = ens_file[0:7] + 'PCP' + ens_file[10:]
      output_path = os.path.join(summary_dir, outname)
   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (tt == 0):
      dx = fin.DX                                             #east-west grid spacing (m)
      dy = fin.DY                                             #north-south grid spacing (m)
      cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
      cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
      stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
      true_lat_1 = fin.TRUELAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
      true_lat_2 = fin.TRUELAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)
      projection = fin.PROJECTION                             #domain projection
      init_time_seconds = fin.INIT_TIME_SECONDS               #initialization time in seconds from 0000 UTC of case
      valid_time_seconds = fin.VALID_TIME_SECONDS             #valid time in seconds from 0000 UTC of case
      forecast_timestep = fin.FORECAST_TIME_STEP              #index of forecast time step

      xlat = fin.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["xlon"][:,:]                     #longitude (dec deg; Lambert conformal)
      hgt = fin.variables["hgt"][:,:]                       #terrain height above MSL (m)

      ny = xlat.shape[0]
      nx = xlat.shape[1]

      rain               = fin.variables["rain"][:,:,:]
      comp_dz		 = fin.variables["comp_dz"][:,:,:]
      rain_15m           = fin.variables["rain"][:,:,:]
      rain_hourly        = fin.variables["rain"][:,:,:]
      rain_3hr		 = fin.variables["rain"][:,:,:]
      rain_6hr		 = fin.variables["rain"][:,:,:]
      comp_dz_hourly	 = fin.variables["comp_dz"][:,:,:]

#### Calc probability matched mean for reflectivity only ####

      temp_dz				 = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz			 = np.mean(temp_dz, axis=0)
      pmm_dz				 = temp_mean_dz * 0.

      pmm_dz[:,:]			 = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

   else:

###

      rain_temp                          = fin.variables["rain"][:,:,:]
      rain                               = rain + rain_temp

      if ((tt > 1) and ((tt % 3) == 1)):
         rain_15m                        = fin.variables["rain"][:,:,:]
      else:
         rain_15m                        = rain_15m + rain_temp
      if ((tt > 1) and ((tt % 12) == 1)): 		#Hourly
         rain_hourly                     = fin.variables["rain"][:,:,:]
      else:
         rain_hourly                     = rain_hourly + rain_temp
      if ((tt > 1) and ((tt % 36) == 1)):
	      rain_3hr			 = fin.variables["rain"][:,:,:]
      else:
         rain_3hr			 = rain_3hr + rain_temp
      if ((tt > 1) and ((tt % 72) == 1)):
	      rain_6hr			 = fin.variables["rain"][:,:,:]
      else:
         rain_6hr			 = rain_6hr + rain_temp

      comp_dz_temp 			 = fin.variables["comp_dz"][:,:,:]
#### Calc probability matched mean for reflectivity only ####

      temp_dz                            = np.where(comp_dz > 100000., 0., comp_dz_temp)
      temp_mean_dz                       = np.mean(temp_dz, axis=0)
      pmm_dz                             = temp_mean_dz * 0.

      pmm_dz[:,:]                        = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

      comp_dz                            = np.where((comp_dz_temp > comp_dz), comp_dz_temp, comp_dz)
      if ((tt > 1) and ((tt % 12) == 1)):
         comp_dz_hourly                  = comp_dz_temp
      else:
         comp_dz_hourly                  = np.where((comp_dz_temp > comp_dz_hourly), comp_dz_temp, comp_dz_hourly)


   fin.close()
   del fin
###

rain_half_90, rain_half_max, rain_half_med, rain_half_3km, rain_half_15km, rain_half_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)
rain_half_90_15m, rain_half_max_15m, rain_half_med_15m, rain_half_3km_15m, rain_half_15km_15m, rain_half_27km_15m = calc_ens_products(rain_15m, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)
rain_half_90_hourly, rain_half_max_hourly, rain_half_med_hourly, rain_half_3km_hourly, rain_half_15km_hourly, rain_half_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)
rain_half_90_3hr, rain_half_max_3hr, rain_half_med_3hr, rain_half_3km_3hr, rain_half_15km_3hr, rain_half_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_1_3hr_thresh)
rain_half_90_6hr, rain_half_max_6hr, rain_half_med_6hr, rain_half_3km_6hr, rain_half_15km_6hr, rain_half_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_27km, radius_max_27km, kernel, rain_1_6hr_thresh)

rain_one_90, rain_one_max, rain_one_med, rain_one_3km, rain_one_15km, rain_one_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)
rain_one_90_15m, rain_one_max_15m, rain_one_med_15m, rain_one_3km_15m, rain_one_15km_15m, rain_one_27km_15m = calc_ens_products(rain_15m, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)
rain_one_90_hourly, rain_one_max_hourly, rain_one_med_hourly, rain_one_3km_hourly, rain_one_15km_hourly, rain_one_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)
rain_one_90_3hr, rain_one_max_3hr, rain_one_med_3hr, rain_one_3km_3hr, rain_one_15km_3hr, rain_one_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_2_3hr_thresh)
rain_one_90_6hr, rain_one_max_6hr, rain_one_med_6hr, rain_one_3km_6hr, rain_one_15km_6hr, rain_one_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_15km, radius_max_27km, kernel, rain_2_6hr_thresh)

rain_two_90, rain_two_max, rain_two_med, rain_two_3km, rain_two_15km, rain_two_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)
rain_two_90_15m, rain_two_max_15m, rain_two_med_15m, rain_two_3km_15m, rain_two_15km_15m, rain_two_27km_15m = calc_ens_products(rain_15m, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)
rain_two_90_hourly, rain_two_max_hourly, rain_two_med_hourly, rain_two_3km_hourly, rain_two_15km_hourly, rain_two_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)
rain_two_90_3hr, rain_two_max_3hr, rain_two_med_3hr, rain_two_3km_3hr, rain_two_15km_3hr, rain_two_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_3_3hr_thresh)
rain_two_90_6hr, rain_two_max_6hr, rain_two_med_6hr, rain_two_3km_6hr, rain_two_15km_6hr, rain_two_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_15km, radius_max_27km, kernel, rain_3_6hr_thresh)

rain_three_90, rain_three_max, rain_three_med, rain_three_3km, rain_three_15km, rain_three_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_4_thresh)

rain_five_90, rain_five_max, rain_five_med, rain_five_3km, rain_five_15km, rain_five_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_5_thresh)

comp_dz_90, comp_dz_max, comp_dz_med, comp_dz_3km, comp_dz_15km, comp_dz_27km = calc_ens_products(comp_dz, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)
comp_dz_90_hourly, comp_dz_max_hourly, comp_dz_med_hourly, comp_dz_3km_hourly, comp_dz_15km_hourly, comp_dz_27km_hourly = calc_ens_products(comp_dz_hourly, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   print "Creating %s ...."% output_path
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NE', ne)
fout.createDimension('NX', nx)
fout.createDimension('NY', ny)

### Set Attributes: ###

setattr(fout,'DX',dx)
setattr(fout,'DY',dy)
setattr(fout,'CEN_LAT',cen_lat)
setattr(fout,'CEN_LON',cen_lon)
setattr(fout,'STAND_LON',stand_lon)
setattr(fout,'TRUE_LAT1',true_lat_1)
setattr(fout,'TRUE_LAT2',true_lat_2)
setattr(fout,'PROJECTION','Lambert Conformal')
setattr(fout,'INIT_TIME_SECONDS',init_time_seconds)
setattr(fout,'VALID_TIME_SECONDS',valid_time_seconds)
setattr(fout,'FORECAST_TIME_STEP',t)

### Create variables ###

xlat1 = fout.createVariable('xlat', 'f4', ('NY','NX',))
xlat1.long_name = "Latitude"
xlat1.units = "degrees North"

xlon1 = fout.createVariable('xlon', 'f4', ('NY','NX',))
xlon1.long_name = "Longitude"
xlon1.units = "degrees West"

hgt1 = fout.createVariable('hgt', 'f4', ('NY','NX',))
hgt1.long_name = "Height AGL"
hgt1.units = "m"

### Member Rainfall ###

rain_memp = fout.createVariable('rain_mem', 'f4', ('NE', 'NY', 'NX',))
rain_memp.long_name = "Member accumulated rainfall"
rain_memp.units = "in"

### 90th percentile and max variables ("p" is placeholder to differentiate variable name from one containing data ) ###

rain_90p = fout.createVariable('rain_90', 'f4', ('NY','NX',))
rain_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_90p_15m = fout.createVariable('rain_90_15m', 'f4', ('NY','NX',))
rain_90p_15m.long_name = "15-min Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p_15m.units = "in"

rain_90p_hourly = fout.createVariable('rain_90_hourly', 'f4', ('NY','NX',))
rain_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p_hourly.units = "in"

rain_90p_3hr = fout.createVariable('rain_90_3hr', 'f4', ('NY','NX',))
rain_90p_3hr.long_name = "3-hr Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_90p_6hr = fout.createVariable('rain_90_6hr', 'f4', ('NY','NX',))
rain_90p_6hr.long_name = "6-hr Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_maxp = fout.createVariable('rain_max', 'f4', ('NY','NX',))
rain_maxp.long_name = "Accumulated ensemble max value of accumulated rainfall"
rain_maxp.units = "in"

rain_maxp_15m = fout.createVariable('rain_max_15m', 'f4', ('NY','NX',))
rain_maxp_15m.long_name = "15-min Accumulated ensemble max value of accumulated rainfall"
rain_maxp_15m.units = "in"

rain_maxp_hourly = fout.createVariable('rain_max_hourly', 'f4', ('NY','NX',))
rain_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of accumulated rainfall"
rain_maxp_hourly.units = "in"

rain_maxp_3hr = fout.createVariable('rain_max_3hr', 'f4', ('NY','NX',))
rain_maxp_3hr.long_name = "3-hr Accumulated ensemble max value of accumulated rainfall"
rain_maxp_3hr.units = "in"

rain_maxp_6hr = fout.createVariable('rain_max_6hr', 'f4', ('NY','NX',))
rain_maxp_6hr.long_name = "6-hr Accumulated ensemble max value of accumulated rainfall"
rain_maxp_6hr.units = "in"

rain_medp = fout.createVariable('rain_med', 'f4', ('NY','NX',))
rain_medp.long_name = "Accumulated ensemble median value of accumulated rainfall"
rain_medp.units = "in"

rain_medp_15m = fout.createVariable('rain_med_15m', 'f4', ('NY','NX',))
rain_medp_15m.long_name = "15-min Accumulated ensemble median value of accumulated rainfall"
rain_medp_15m.units = "in"

rain_medp_hourly = fout.createVariable('rain_med_hourly', 'f4', ('NY','NX',))
rain_medp_hourly.long_name = "1-hr Accumulated ensemble median value of accumulated rainfall"
rain_medp_hourly.units = "in"

rain_medp_3hr = fout.createVariable('rain_med_3hr', 'f4', ('NY','NX',))
rain_medp_3hr.long_name = "3-hr Accumulated ensemble median value of accumulated rainfall"
rain_medp_3hr.units = "in"

rain_medp_6hr = fout.createVariable('rain_med_6hr', 'f4', ('NY','NX',))
rain_medp_6hr.long_name = "6-hr Accumulated ensemble median value of accumulated rainfall"
rain_medp_6hr.units = "in"

comp_dz_90p = fout.createVariable('comp_dz_90', 'f4', ('NY','NX',))
comp_dz_90p.long_name = "Accumulated ensemble 90th percentile value of composite reflectivity"
comp_dz_90p.units = "dBZ"

comp_dz_90p_hourly = fout.createVariable('comp_dz_90_hourly', 'f4', ('NY','NX',))
comp_dz_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of composite reflectivity"
comp_dz_90p_hourly.units = "dBZ"

comp_dz_maxp = fout.createVariable('comp_dz_max', 'f4', ('NY','NX',))
comp_dz_maxp.long_name = "Accumulated ensemble max value of composite reflectivity"
comp_dz_maxp.units = "dBZ"

comp_dz_maxp_hourly = fout.createVariable('comp_dz_max_hourly', 'f4', ('NY','NX',))
comp_dz_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of composite reflectivity"
comp_dz_maxp_hourly.units = "dBZ"

comp_dz_medp = fout.createVariable('comp_dz_med', 'f4', ('NY','NX',))
comp_dz_medp.long_name = "Accumulated ensemble max value of composite reflectivity"
comp_dz_medp.units = "dBZ"

comp_dz_medp_hourly = fout.createVariable('comp_dz_med_hourly', 'f4', ('NY','NX',))
comp_dz_medp_hourly.long_name = "1-hr Accumulated ensemble max value of composite reflectivity"
comp_dz_medp_hourly.units = "dBZ"

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

### Probability of Exceedance ###

### Rainfall 0.5 inch threshold ###

rain_half_prob_3 = fout.createVariable('rain_half_prob_3km', 'f4', ('NY','NX',))
rain_half_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_1_thresh
rain_half_prob_3.units = "%"

rain_half_prob_3_15m = fout.createVariable('rain_half_prob_3km_15m', 'f4', ('NY','NX',))
rain_half_prob_3_15m.long_name = "15-min Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_1_thresh
rain_half_prob_3_15m.units = "%"

rain_half_prob_3_hourly = fout.createVariable('rain_half_prob_3km_hourly', 'f4', ('NY','NX',))
rain_half_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_1_thresh
rain_half_prob_3_hourly.units = "%"

rain_half_prob_3_3hr = fout.createVariable('rain_half_prob_3km_3hr', 'f4', ('NY','NX',))
rain_half_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % rain_1_3hr_thresh
rain_half_prob_3_3hr.units = "%"

rain_half_prob_3_6hr = fout.createVariable('rain_half_prob_3km_6hr', 'f4', ('NY','NX',))
rain_half_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_1_6hr_thresh
rain_half_prob_3_6hr.units = "%"

rain_half_prob_15 = fout.createVariable('rain_half_prob_15km', 'f4', ('NY','NX',))
rain_half_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_1_thresh
rain_half_prob_15.units = "%"

rain_half_prob_15_15m = fout.createVariable('rain_half_prob_15km_15m', 'f4', ('NY','NX',))
rain_half_prob_15_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_1_thresh
rain_half_prob_15_15m.units = "%"

rain_half_prob_15_hourly = fout.createVariable('rain_half_prob_15km_hourly', 'f4', ('NY','NX',))
rain_half_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_1_thresh
rain_half_prob_15_hourly.units = "%"

rain_half_prob_15_3hr = fout.createVariable('rain_half_prob_15km_3hr', 'f4', ('NY','NX',))
rain_half_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % rain_1_3hr_thresh
rain_half_prob_15_3hr.units = "%"

rain_half_prob_15_6hr = fout.createVariable('rain_half_prob_15km_6hr', 'f4', ('NY','NX',))
rain_half_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_1_6hr_thresh
rain_half_prob_15_6hr.units = "%"

rain_half_prob_27 = fout.createVariable('rain_half_prob_27km', 'f4', ('NY','NX',))
rain_half_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_1_thresh
rain_half_prob_27.units = "%"

rain_half_prob_27_15m = fout.createVariable('rain_half_prob_27km_15m', 'f4', ('NY','NX',))
rain_half_prob_27_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" %rain_1_thresh
rain_half_prob_27_15m.units = "%"

rain_half_prob_27_hourly = fout.createVariable('rain_half_prob_27km_hourly', 'f4', ('NY','NX',))
rain_half_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" %rain_1_thresh
rain_half_prob_27_hourly.units = "%"

rain_half_prob_27_3hr = fout.createVariable('rain_half_prob_27km_3hr', 'f4', ('NY','NX',))
rain_half_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % rain_1_3hr_thresh
rain_half_prob_27_3hr.units = "%"

rain_half_prob_27_6hr = fout.createVariable('rain_half_prob_27km_6hr', 'f4', ('NY','NX',))
rain_half_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_1_6hr_thresh
rain_half_prob_27_6hr.units = "%"

### Rainfall 1 inch threshold ###

rain_one_prob_3 = fout.createVariable('rain_one_prob_3km', 'f4', ('NY','NX',))
rain_one_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % rain_2_thresh
rain_one_prob_3.units = "%"

rain_one_prob_3_15m = fout.createVariable('rain_one_prob_3km_15m', 'f4', ('NY','NX',))
rain_one_prob_3_15m.long_name = "15-min Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % rain_2_thresh
rain_one_prob_3_15m.units = "%"

rain_one_prob_3_hourly = fout.createVariable('rain_one_prob_3km_hourly', 'f4', ('NY','NX',))
rain_one_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % rain_2_thresh
rain_one_prob_3_hourly.units = "%"

rain_one_prob_3_3hr = fout.createVariable('rain_one_prob_3km_3hr', 'f4', ('NY','NX',))
rain_one_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_2_3hr_thresh
rain_one_prob_3_3hr.units = "%"

rain_one_prob_3_6hr = fout.createVariable('rain_one_prob_3km_6hr', 'f4', ('NY','NX',))
rain_one_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_2_6hr_thresh
rain_one_prob_3_6hr.units = "%"

rain_one_prob_15 = fout.createVariable('rain_one_prob_15km', 'f4', ('NY','NX',))
rain_one_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % rain_2_thresh
rain_one_prob_15.units = "%"

rain_one_prob_15_15m = fout.createVariable('rain_one_prob_15km_15m', 'f4', ('NY','NX',))
rain_one_prob_15_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % rain_2_thresh
rain_one_prob_15_15m.units = "%"

rain_one_prob_15_hourly = fout.createVariable('rain_one_prob_15km_hourly', 'f4', ('NY','NX',))
rain_one_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % rain_2_thresh
rain_one_prob_15_hourly.units = "%"

rain_one_prob_15_3hr = fout.createVariable('rain_one_prob_15km_3hr', 'f4', ('NY','NX',))
rain_one_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_2_3hr_thresh
rain_one_prob_15_3hr.units = "%"

rain_one_prob_15_6hr = fout.createVariable('rain_one_prob_15km_6hr', 'f4', ('NY','NX',))
rain_one_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_2_6hr_thresh
rain_one_prob_15_6hr.units = "%"

rain_one_prob_27 = fout.createVariable('rain_one_prob_27km', 'f4', ('NY','NX',))
rain_one_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % rain_2_thresh
rain_one_prob_27.units = "%"

rain_one_prob_27_15m = fout.createVariable('rain_one_prob_27km_15m', 'f4', ('NY','NX',))
rain_one_prob_27_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % rain_2_thresh
rain_one_prob_27_15m.units = "%"

rain_one_prob_27_hourly = fout.createVariable('rain_one_prob_27km_hourly', 'f4', ('NY','NX',))
rain_one_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % rain_2_thresh
rain_one_prob_27_hourly.units = "%"

rain_one_prob_27_3hr = fout.createVariable('rain_one_prob_27km_3hr', 'f4', ('NY','NX',))
rain_one_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_2_3hr_thresh
rain_one_prob_27_3hr.units = "%"

rain_one_prob_27_6hr = fout.createVariable('rain_one_prob_27km_6hr', 'f4', ('NY','NX',))
rain_one_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_2_6hr_thresh
rain_one_prob_27_6hr.units = "%"

### Rainfall 2 inch threshold ###

rain_two_prob_3 = fout.createVariable('rain_two_prob_3km', 'f4', ('NY','NX',))
rain_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_3_thresh
rain_two_prob_3.units = "%"

rain_two_prob_3_15m = fout.createVariable('rain_two_prob_3km_15m', 'f4', ('NY','NX',))
rain_two_prob_3_15m.long_name = "15-min Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_3_thresh
rain_two_prob_3_15m.units = "%"

rain_two_prob_3_hourly = fout.createVariable('rain_two_prob_3km_hourly', 'f4', ('NY','NX',))
rain_two_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_3_thresh
rain_two_prob_3_hourly.units = "%"

rain_two_prob_3_3hr = fout.createVariable('rain_two_prob_3km_3hr', 'f4', ('NY','NX',))
rain_two_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_3_3hr_thresh
rain_two_prob_3_3hr.units = "%"

rain_two_prob_3_6hr = fout.createVariable('rain_two_prob_3km_6hr', 'f4', ('NY','NX',))
rain_two_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_3_6hr_thresh
rain_two_prob_3_6hr.units = "%"

rain_two_prob_15 = fout.createVariable('rain_two_prob_15km', 'f4', ('NY','NX',))
rain_two_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_3_thresh
rain_two_prob_15.units = "%"

rain_two_prob_15_15m = fout.createVariable('rain_two_prob_15km_15m', 'f4', ('NY','NX',))
rain_two_prob_15_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_3_thresh
rain_two_prob_15_15m.units = "%"

rain_two_prob_15_hourly = fout.createVariable('rain_two_prob_15km_hourly', 'f4', ('NY','NX',))
rain_two_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_3_thresh
rain_two_prob_15_hourly.units = "%"

rain_two_prob_15_3hr = fout.createVariable('rain_two_prob_15km_3hr', 'f4', ('NY','NX',))
rain_two_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_3_3hr_thresh
rain_two_prob_15_3hr.units = "%"

rain_two_prob_15_6hr = fout.createVariable('rain_two_prob_15km_6hr', 'f4', ('NY','NX',))
rain_two_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_3_6hr_thresh
rain_two_prob_15_6hr.units = "%"

rain_two_prob_27 = fout.createVariable('rain_two_prob_27km', 'f4', ('NY','NX',))
rain_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_3_thresh
rain_two_prob_27.units = "%"

rain_two_prob_27_15m = fout.createVariable('rain_two_prob_27km_15m', 'f4', ('NY','NX',))
rain_two_prob_27_15m.long_name = "15-min Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_3_thresh
rain_two_prob_27_15m.units = "%"

rain_two_prob_27_hourly = fout.createVariable('rain_two_prob_27km_hourly', 'f4', ('NY','NX',))
rain_two_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_3_thresh
rain_two_prob_27_hourly.units = "%"

rain_two_prob_27_3hr = fout.createVariable('rain_two_prob_27km_3hr', 'f4', ('NY','NX',))
rain_two_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_3_3hr_thresh
rain_two_prob_27_3hr.units = "%"

rain_two_prob_27_6hr = fout.createVariable('rain_two_prob_27km_6hr', 'f4', ('NY','NX',))
rain_two_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_3_6hr_thresh
rain_two_prob_27_6hr.units = "%"

### Rainfall 3 inch threshold (only 5-minute products) ###

rain_three_prob_3 = fout.createVariable('rain_three_prob_3km', 'f4', ('NY','NX',))
rain_three_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_4_thresh
rain_three_prob_3.units = "%"

rain_three_prob_15 = fout.createVariable('rain_three_prob_15km', 'f4', ('NY','NX',))
rain_three_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_4_thresh
rain_three_prob_15.units = "%"

rain_three_prob_27 = fout.createVariable('rain_three_prob_27km', 'f4', ('NY','NX',))
rain_three_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_4_thresh
rain_three_prob_27.units = "%"

### Rainfall 5 inch threshold (only 5-minute products) ###

rain_five_prob_3 = fout.createVariable('rain_five_prob_3km', 'f4', ('NY','NX',))
rain_five_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % rain_5_thresh
rain_five_prob_3.units = "%"

rain_five_prob_15 = fout.createVariable('rain_five_prob_15km', 'f4', ('NY','NX',))
rain_five_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % rain_5_thresh
rain_five_prob_15.units = "%"

rain_five_prob_27 = fout.createVariable('rain_five_prob_27km', 'f4', ('NY','NX',))
rain_five_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % rain_5_thresh
rain_five_prob_27.units = "%"

### Comp DZ ###

comp_dz_prob_3 = fout.createVariable('comp_dz_prob_3km', 'f4', ('NY','NX',))
comp_dz_prob_3.long_name = "Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3.units = "%"

comp_dz_prob_3_hourly = fout.createVariable('comp_dz_prob_3km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3_hourly.units = "%"

comp_dz_prob_15 = fout.createVariable('comp_dz_prob_15km', 'f4', ('NY','NX',))
comp_dz_prob_15.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (15 km neighborhood)"
comp_dz_prob_15.units = "%"

comp_dz_prob_15_hourly = fout.createVariable('comp_dz_prob_15km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (15 km neighborhood)"
comp_dz_prob_15_hourly.units = "%"

comp_dz_prob_27 = fout.createVariable('comp_dz_prob_27km', 'f4', ('NY','NX',))
comp_dz_prob_27.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27.units = "%"

comp_dz_prob_27_hourly = fout.createVariable('comp_dz_prob_27km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27_hourly.units = "%"


### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['rain_mem'][:] = rain

fout.variables['rain_90'][:] = rain_half_90
fout.variables['rain_90_15m'][:] = rain_half_90_15m
fout.variables['rain_90_hourly'][:] = rain_half_90_hourly
fout.variables['rain_90_3hr'][:] = rain_half_90_3hr
fout.variables['rain_90_6hr'][:] = rain_half_90_6hr

fout.variables['rain_max'][:] = rain_half_max
fout.variables['rain_max_15m'][:] = rain_half_max_15m
fout.variables['rain_max_hourly'][:] = rain_half_max_hourly
fout.variables['rain_max_3hr'][:] = rain_half_max_3hr
fout.variables['rain_max_6hr'][:] = rain_half_max_6hr

fout.variables['rain_med'][:] = rain_half_med
fout.variables['rain_med_15m'][:] = rain_half_med_15m
fout.variables['rain_med_hourly'][:] = rain_half_med_hourly
fout.variables['rain_med_3hr'][:] = rain_half_med_3hr
fout.variables['rain_med_6hr'][:] = rain_half_med_6hr

fout.variables['comp_dz_90'][:] = comp_dz_90
fout.variables['comp_dz_90_hourly'][:] = comp_dz_90_hourly
fout.variables['comp_dz_max'][:] = comp_dz_max
fout.variables['comp_dz_max_hourly'][:] = comp_dz_max_hourly
fout.variables['comp_dz_med'][:] = comp_dz_med
fout.variables['comp_dz_med_hourly'][:] = comp_dz_med_hourly

fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['rain_half_prob_3km'][:] = rain_half_3km
fout.variables['rain_half_prob_3km_15m'][:] = rain_half_3km_15m
fout.variables['rain_half_prob_3km_hourly'][:] = rain_half_3km_hourly
fout.variables['rain_half_prob_3km_3hr'][:] = rain_half_3km_3hr
fout.variables['rain_half_prob_3km_6hr'][:] = rain_half_3km_6hr

fout.variables['rain_half_prob_15km'][:] = rain_half_15km
fout.variables['rain_half_prob_15km_15m'][:] = rain_half_15km_15m
fout.variables['rain_half_prob_15km_hourly'][:] = rain_half_15km_hourly
fout.variables['rain_half_prob_15km_3hr'][:] = rain_half_15km_3hr
fout.variables['rain_half_prob_15km_6hr'][:] = rain_half_15km_6hr

fout.variables['rain_half_prob_27km'][:] = rain_half_27km
fout.variables['rain_half_prob_27km_15m'][:] = rain_half_27km_15m
fout.variables['rain_half_prob_27km_hourly'][:] = rain_half_27km_hourly
fout.variables['rain_half_prob_27km_3hr'][:] = rain_half_27km_3hr
fout.variables['rain_half_prob_27km_6hr'][:] = rain_half_27km_6hr

fout.variables['rain_one_prob_3km'][:] = rain_one_3km
fout.variables['rain_one_prob_3km_15m'][:] = rain_one_3km_15m
fout.variables['rain_one_prob_3km_hourly'][:] = rain_one_3km_hourly
fout.variables['rain_one_prob_3km_3hr'][:] = rain_one_3km_3hr
fout.variables['rain_one_prob_3km_6hr'][:] = rain_one_3km_6hr

fout.variables['rain_one_prob_15km'][:] = rain_one_15km
fout.variables['rain_one_prob_15km_15m'][:] = rain_one_15km_15m
fout.variables['rain_one_prob_15km_hourly'][:] = rain_one_15km_hourly
fout.variables['rain_one_prob_15km_3hr'][:] = rain_one_15km_3hr
fout.variables['rain_one_prob_15km_6hr'][:] = rain_one_15km_6hr

fout.variables['rain_one_prob_27km'][:] = rain_one_27km
fout.variables['rain_one_prob_27km_15m'][:] = rain_one_27km_15m
fout.variables['rain_one_prob_27km_hourly'][:] = rain_one_27km_hourly
fout.variables['rain_one_prob_27km_3hr'][:] = rain_one_27km_3hr
fout.variables['rain_one_prob_27km_6hr'][:] = rain_one_27km_6hr

fout.variables['rain_two_prob_3km'][:] = rain_two_3km
fout.variables['rain_two_prob_3km_15m'][:] = rain_two_3km_15m
fout.variables['rain_two_prob_3km_hourly'][:] = rain_two_3km_hourly
fout.variables['rain_two_prob_3km_3hr'][:] = rain_two_3km_3hr
fout.variables['rain_two_prob_3km_6hr'][:] = rain_two_3km_6hr

fout.variables['rain_two_prob_15km'][:] = rain_two_15km
fout.variables['rain_two_prob_15km_15m'][:] = rain_two_15km_15m
fout.variables['rain_two_prob_15km_hourly'][:] = rain_two_15km_hourly
fout.variables['rain_two_prob_15km_3hr'][:] = rain_two_15km_3hr
fout.variables['rain_two_prob_15km_6hr'][:] = rain_two_15km_6hr

fout.variables['rain_two_prob_27km'][:] = rain_two_27km
fout.variables['rain_two_prob_27km_15m'][:] = rain_two_27km_15m
fout.variables['rain_two_prob_27km_hourly'][:] = rain_two_27km_hourly
fout.variables['rain_two_prob_27km_3hr'][:] = rain_two_27km_3hr
fout.variables['rain_two_prob_27km_6hr'][:] = rain_two_27km_6hr

fout.variables['rain_three_prob_3km'][:] = rain_three_3km
fout.variables['rain_three_prob_15km'][:] = rain_three_15km
fout.variables['rain_three_prob_27km'][:] = rain_three_27km

fout.variables['rain_five_prob_3km'][:] = rain_five_3km
fout.variables['rain_five_prob_15km'][:] = rain_five_15km
fout.variables['rain_five_prob_27km'][:] = rain_five_27km

fout.variables['comp_dz_prob_3km'][:] = comp_dz_3km
fout.variables['comp_dz_prob_3km_hourly'][:] = comp_dz_3km_hourly
fout.variables['comp_dz_prob_15km'][:] = comp_dz_15km
fout.variables['comp_dz_prob_15km_hourly'][:] = comp_dz_15km_hourly
fout.variables['comp_dz_prob_27km'][:] = comp_dz_27km
fout.variables['comp_dz_prob_27km_hourly'][:] = comp_dz_27km_hourly

### Close output file ###

fout.close()
del fout
