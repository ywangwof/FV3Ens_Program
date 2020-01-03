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

comp_dz_thresh		   = 40.
############################ Set Thresholds for convolution and neighborhood probabilities: #################################
rain_1_thresh              = 0.5		#0.5 inches
rain_1_3hr_thresh	   = 1.			#1 inch, 3 hr
rain_1_6hr_thresh	   = 2			#2 inches, 6 hr

rain_2_thresh              = 1.			#1 inch
rain_2_3hr_thresh	   = 2.			#2 inches, 3 hr
rain_2_6hr_thresh	   = 3.			#3 inches, 6 hr

rain_3_thresh              = 2.			#2 inches
rain_3_3hr_thresh	   = 3.			#3 inches, 3 hr
rain_3_6hr_thresh	   = 5.			#5 inches, 6 hr

############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18
start_t = 0

ens_files = []
swt_files = []
summary_files_temp = os.listdir(summary_dir)

print(summary_files_temp)
for swt_t in range(1, fcst_nt):
   swt_found = 0
   str_t = str(swt_t)
   if (len(str_t) == 1):
      str_t = '0' + str_t

   for f, file in enumerate(summary_files_temp): 
      if ((file[-28:-25] == 'PCP') and (file[-24:-22] == str_t)):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
         swt_found = 1

   if (swt_found == 0):
      break
   else:
      if (((swt_t % 12) == 0) and (swt_t < fcst_nt)):
         start_t = swt_t
         time.sleep(2)

################################# Find and process ENS and PCP files: ################################

for f, file in enumerate(summary_files_temp):
   if (file[-28:-25] == 'ENS'):
      ens_files.append(file)
   if (file[-28:-25] == 'PCP'):
      swt_files.append(file)
      
ens_files.sort() 
swt_files.sort()

print start_t, 'START'

for tt in range(start_t, t+1):

   if (start_t > 0):
      swt_file = swt_files[start_t]
      swt_infile = os.path.join(summary_dir, swt_file)

   print(ens_files)
   ens_file = ens_files[tt]
   infile = os.path.join(summary_dir, ens_file)
   
   if (tt == t): 
   ### Set output path ###
      outname = ens_file[0:7] + 'PCP' + ens_file[10:]
      summry_dir = '/scratch/brian.matilla/'
      output_path = summry_dir + outname

   print(start_t)
   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (tt == start_t):

      if (start_t > 0):
         try:
            swt_fin = netCDF4.Dataset(swt_infile, "r")
            print "Opening %s \n" % swt_infile
         except:
            print "%s does not exist! \n" %swt_infile
            sys.exit(1)

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

      xlat = fin.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["xlon"][:,:]                     #longitude (dec deg; Lambert conformal)
      hgt = fin.variables["hgt"][:,:]                       #terrain height above MSL (m)

      ny = xlat.shape[0]
      nx = xlat.shape[1]

      if (start_t == 0): 
         rain                            = fin.variables["rain"][:,:,:]
         comp_dz                         = fin.variables["comp_dz"][:,:,:]
      else:
         rain                            = swt_fin.variables["rain_out"][:,:,:]
         comp_dz                         = swt_fin.variables["comp_dz_out"][:,:,:]

      rain_hourly                        = fin.variables["rain"][:,:,:]
      rain_3hr				 = fin.variables["rain"][:,:,:]
      rain_6hr				 = fin.variables["rain"][:,:,:]
      comp_dz_hourly			 = fin.variables["comp_dz"][:,:,:]

#### Calc probability matched mean for reflectivity only ####
      temp_dz = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

      if (start_t > 0):
         swt_fin.close()
         del swt_fin

   else:
      comp_dz_temp			 = fin.variables["comp_dz"][:,:,:]

### Calc probability matched mean for reflectivity only ###

      temp_dz = np.where(comp_dz_temp > 100000., 0., comp_dz_temp)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

      comp_dz                            = np.where((comp_dz_temp > comp_dz), comp_dz_temp, comp_dz)
      if ((tt > 1) and ((tt % 12) == 1)):
         comp_dz_hourly                  = comp_dz_temp
      else:
         comp_dz_hourly                  = np.where((comp_dz_temp > comp_dz_hourly), comp_dz_temp, comp_dz_hourly)

###
      rain_temp                          = fin.variables["rain"][:,:,:]
      rain                               = rain + rain_temp
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
   

   fin.close()
   del fin   
###

### Output swaths at hourly intervals to speed up processing: 

if ((t > 1) and ((t % 12) == 0)):
   rain_out = rain
   comp_dz_out = comp_dz

#### Calculate ensemble output variables: ####

rain_half_90, rain_half_max, rain_half_med, rain_half_3km, rain_half_15km, rain_half_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)
rain_half_90_hourly, rain_half_max_hourly, rain_half_med_hourly, rain_half_3km_hourly, rain_half_15km_hourly, rain_half_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)
rain_half_90_3hr, rain_half_max_3hr, rain_half_med_3hr, rain_half_3km_3hr, rain_half_15km_3hr, rain_half_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_1_3hr_thresh)
rain_half_90_6hr, rain_half_max_6hr, rain_half_med_6hr, rain_half_3km_6hr, rain_half_15km_6hr, rain_half_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_27km, radius_max_27km, kernel, rain_1_6hr_thresh)

rain_one_90, rain_one_max, rain_one_med, rain_one_3km, rain_one_15km, rain_one_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)
rain_one_90_hourly, rain_one_max_hourly, rain_one_med_hourly, rain_one_3km_hourly, rain_one_15km_hourly, rain_one_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)
rain_one_90_3hr, rain_one_max_3hr, rain_one_med_3hr, rain_one_3km_3hr, rain_one_15km_3hr, rain_one_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_2_3hr_thresh)
rain_one_90_6hr, rain_one_max_6hr, rain_one_med_6hr, rain_one_3km_6hr, rain_one_15km_6hr, rain_one_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_15km, radius_max_27km, kernel, rain_2_6hr_thresh)

rain_two_90, rain_two_max, rain_two_med, rain_two_3km, rain_two_15km, rain_two_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)
rain_two_90_hourly, rain_two_max_hourly, rain_two_med_hourly, rain_two_3km_hourly, rain_two_15km_hourly, rain_two_27km_hourly = calc_ens_products(rain_hourly, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)
rain_two_90_3hr, rain_two_max_3hr, rain_two_med_3hr, rain_two_3km_3hr, rain_two_15km_3hr, rain_two_27km_3hr = calc_ens_products(rain_3hr, 1, radius_max_15km, radius_max_27km, kernel, rain_3_3hr_thresh)
rain_two_90_6hr, rain_two_max_6hr, rain_two_med_6hr, rain_two_3km_6hr, rain_two_15km_6hr, rain_two_27km_6hr = calc_ens_products(rain_6hr, 1, radius_max_15km, radius_max_27km, kernel, rain_3_6hr_thresh)

comp_dz_90, comp_dz_max, comp_dz_3km, comp_dz_15km, comp_dz_27km = calc_ens_products(comp_dz, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)
comp_dz_90_hourly, comp_dz_max_hourly, comp_dz_3km_hourly, comp_dz_15km_hourly, comp_dz_27km_hourly = calc_ens_products(comp_dz_hourly, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)

#### handle hourly swath files ####

if ((t > 1) and ((t % 12) == 0)):
   fout.createDimension('NE', ne)

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

### handle hourly swath files:

if ((t > 1) and ((t % 12) == 0)):
   comp_dz_hourly_out = fout.createVariable('comp_dz_out', 'f4', ('NE', 'NY', 'NX',))
   rain_hourly_out = fout.createVariable('rain_hourly_out', 'f4', ('NE', 'NY', 'NX',))
   rain_3hr_out = fout.createVariable('rain_3hr_out', 'f4', ('NE', 'NY', 'NX',))
   rain_6hr_out = fout.createVariable('rain_6hr_out', 'f4', ('NE', 'NY', 'NX',))

### 90th percentile and max variables ("p" is placeholder to differentiate variable name from one containing data ) ###

rain_90p = fout.createVariable('rain_90', 'f4', ('NY','NX',))
rain_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

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

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

### Probability of Exceedance ###

### Rainfall 0.5 inch threshold ###

rain_half_prob_3 = fout.createVariable('rain_half_prob_3km', 'f4', ('NY','NX',))
rain_half_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_1_thresh)
rain_half_prob_3.units = "%"

rain_half_prob_3_hourly = fout.createVariable('rain_half_prob_3km_hourly', 'f4', ('NY','NX',))
rain_half_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_1_thresh)
rain_half_prob_3_hourly.units = "%"

rain_half_prob_3_3hr = fout.createVariable('rain_half_prob_3km_3hr', 'f4', ('NY','NX',))
rain_half_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % str(rain_1_3hr_thresh)
rain_half_prob_3_3hr.units = "%"

rain_half_prob_3_6hr = fout.createVariable('rain_half_prob_3km_6hr', 'f4', ('NY','NX',))
rain_half_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_1_6hr_thresh)
rain_half_prob_3_6hr.units = "%"

rain_half_prob_15 = fout.createVariable('rain_half_prob_15km', 'f4', ('NY','NX',))
rain_half_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_1_thresh) 
rain_half_prob_15.units = "%"

rain_half_prob_15_hourly = fout.createVariable('rain_half_prob_15km_hourly', 'f4', ('NY','NX',))
rain_half_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_1_thresh)
rain_half_prob_15_hourly.units = "%"

rain_half_prob_15_3hr = fout.createVariable('rain_half_prob_15km_3hr', 'f4', ('NY','NX',))
rain_half_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % str(rain_1_3hr_thresh)
rain_half_prob_15_3hr.units = "%"

rain_half_prob_15_6hr = fout.createVariable('rain_half_prob_15km_6hr', 'f4', ('NY','NX',))
rain_half_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_1_6hr_thresh)
rain_half_prob_15_6hr.units = "%"

rain_half_prob_27 = fout.createVariable('rain_half_prob_27km', 'f4', ('NY','NX',))
rain_half_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_1_thresh)
rain_half_prob_27.units = "%"

rain_half_prob_27_hourly = fout.createVariable('rain_half_prob_27km_hourly', 'f4', ('NY','NX',))
rain_half_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_1_thresh)
rain_half_prob_27_hourly.units = "%"

rain_half_prob_27_3hr = fout.createVariable('rain_half_prob_27km_3hr', 'f4', ('NY','NX',))
rain_half_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % str(rain_1_3hr_thresh)
rain_half_prob_27_3hr.units = "%"

rain_half_prob_27_6hr = fout.createVariable('rain_half_prob_27km_6hr', 'f4', ('NY','NX',))
rain_half_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_1_6hr_thresh)
rain_half_prob_27_6hr.units = "%"

### Rainfall 1 inch threshold ###

rain_one_prob_3 = fout.createVariable('rain_one_prob_3km', 'f4', ('NY','NX',))
rain_one_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % str(rain_2_thresh)
rain_one_prob_3.units = "%"

rain_one_prob_3_hourly = fout.createVariable('rain_one_prob_3km_hourly', 'f4', ('NY','NX',))
rain_one_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inch" % str(rain_2_thresh)
rain_one_prob_3_hourly.units = "%"

rain_one_prob_3_3hr = fout.createVariable('rain_one_prob_3km_3hr', 'f4', ('NY','NX',))
rain_one_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_2_3hr_thresh)
rain_one_prob_3_3hr.units = "%"

rain_one_prob_3_6hr = fout.createVariable('rain_one_prob_3km_6hr', 'f4', ('NY','NX',))
rain_one_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_2_6hr_thresh)
rain_one_prob_3_6hr.units = "%"

rain_one_prob_15 = fout.createVariable('rain_one_prob_15km', 'f4', ('NY','NX',))
rain_one_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % str(rain_2_thresh)
rain_one_prob_15.units = "%"

rain_one_prob_15_hourly = fout.createVariable('rain_one_prob_15km_hourly', 'f4', ('NY','NX',))
rain_one_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inch (15 km neighborhood)" % str(rain_2_thresh)
rain_one_prob_15_hourly.units = "%"

rain_one_prob_15_3hr = fout.createVariable('rain_one_prob_15km_3hr', 'f4', ('NY','NX',))
rain_one_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_2_3hr_thresh)
rain_one_prob_15_3hr.units = "%"

rain_one_prob_15_6hr = fout.createVariable('rain_one_prob_15km_6hr', 'f4', ('NY','NX',))
rain_one_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_2_6hr_thresh)
rain_one_prob_15_6hr.units = "%"

rain_one_prob_27 = fout.createVariable('rain_one_prob_27km', 'f4', ('NY','NX',))
rain_one_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % str(rain_2_thresh)
rain_one_prob_27.units = "%"

rain_one_prob_27_hourly = fout.createVariable('rain_one_prob_27km_hourly', 'f4', ('NY','NX',))
rain_one_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inch (27 km neighborhood)" % str(rain_2_thresh)
rain_one_prob_27_hourly.units = "%"

rain_one_prob_27_3hr = fout.createVariable('rain_one_prob_27km_3hr', 'f4', ('NY','NX',))
rain_one_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_2_3hr_thresh)
rain_one_prob_27_3hr.units = "%"

rain_one_prob_27_6hr = fout.createVariable('rain_one_prob_27km_6hr', 'f4', ('NY','NX',))
rain_one_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_2_6hr_thresh)
rain_one_prob_27_6hr.units = "%"

### Rainfall 2 inch threshold ###

rain_two_prob_3 = fout.createVariable('rain_two_prob_3km', 'f4', ('NY','NX',))
rain_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_3_thresh)
rain_two_prob_3.units = "%"

rain_two_prob_3_hourly = fout.createVariable('rain_two_prob_3km_hourly', 'f4', ('NY','NX',))
rain_two_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_3_thresh)
rain_two_prob_3_hourly.units = "%"

rain_two_prob_3_3hr = fout.createVariable('rain_two_prob_3km_3hr', 'f4', ('NY','NX',))
rain_two_prob_3_3hr.long_name = "3-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_3_3hr_thresh)
rain_two_prob_3_3hr.units = "%"

rain_two_prob_3_6hr = fout.createVariable('rain_two_prob_3km_6hr', 'f4', ('NY','NX',))
rain_two_prob_3_6hr.long_name = "6-hr Accumulated ensemble gridpoint probability of rainfall greater than %s inches" % str(rain_3_6hr_thresh)
rain_two_prob_3_6hr.units = "%"

rain_two_prob_15 = fout.createVariable('rain_two_prob_15km', 'f4', ('NY','NX',))
rain_two_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_3_thresh)
rain_two_prob_15.units = "%"

rain_two_prob_15_hourly = fout.createVariable('rain_two_prob_15km_hourly', 'f4', ('NY','NX',))
rain_two_prob_15_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_3_thresh)
rain_two_prob_15_hourly.units = "%"

rain_two_prob_15_3hr = fout.createVariable('rain_two_prob_15km_3hr', 'f4', ('NY','NX',))
rain_two_prob_15_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_3_3hr_thresh)
rain_two_prob_15_3hr.units = "%"

rain_two_prob_15_6hr = fout.createVariable('rain_two_prob_15km_6hr', 'f4', ('NY','NX',))
rain_two_prob_15_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (15 km neighborhood)" % str(rain_3_6hr_thresh)
rain_two_prob_15_6hr.units = "%"

rain_two_prob_27 = fout.createVariable('rain_two_prob_27km', 'f4', ('NY','NX',))
rain_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_3_thresh)
rain_two_prob_27.units = "%"

rain_two_prob_27_hourly = fout.createVariable('rain_two_prob_27km_hourly', 'f4', ('NY','NX',))
rain_two_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_3_thresh)
rain_two_prob_27_hourly.units = "%"

rain_two_prob_27_3hr = fout.createVariable('rain_two_prob_27km_3hr', 'f4', ('NY','NX',))
rain_two_prob_27_3hr.long_name = "3-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_3_3hr_thresh) 
rain_two_prob_27_3hr.units = "%"

rain_two_prob_27_6hr = fout.createVariable('rain_two_prob_27km_6hr', 'f4', ('NY','NX',))
rain_two_prob_27_6hr.long_name = "6-hr Accumulated ensemble probability of rainfall greater than %s inches (27 km neighborhood)" % str(rain_3_6hr_thresh)
rain_two_prob_27_6hr.units = "%"

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

### handle hourly swath files:
if ((t > 1) and ((t % 12) == 0)):
   fout.variables['rain_3hr_out'][:] = rain_3hr_out
   fout.variables['rain_6hr_out'][:] = rain_6hr_out
   fout.variables['rain_out'] = rain_out
   fout.variables['comp_dz_out'] = comp_dz_out

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['rain_90'][:] = rain_half_90
fout.variables['rain_90_hourly'][:] = rain_half_90_hourly
fout.variables['rain_90_3hr'][:] = rain_half_90_3hr
fout.variables['rain_90_6hr'][:] = rain_half_90_6hr

fout.variables['rain_max'][:] = rain_half_max
fout.variables['rain_max_hourly'][:] = rain_half_max_hourly
fout.variables['rain_max_3hr'][:] = rain_half_max_3hr
fout.variables['rain_max_6hr'][:] = rain_half_max_6hr

fout.variables['rain_med'][:] = rain_half_med
fout.variables['rain_med_hourly'][:] = rain_half_med_hourly
fout.variables['rain_med_3hr'][:] = rain_half_med_3hr
fout.variables['rain_med_6hr'][:] = rain_half_med_6hr

fout.variables['comp_dz_90'][:] = comp_dz_90
fout.variables['comp_dz_90_hourly'][:] = comp_dz_90_hourly
fout.variables['comp_dz_max'][:] = comp_dz_max
fout.variables['comp_dz_max_hourly'][:] = comp_dz_max_hourly
fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['rain_half_prob_3km'][:] = rain_half_3km
fout.variables['rain_half_prob_3km_hourly'][:] = rain_half_3km_hourly
fout.variables['rain_half_prob_3km_3hr'][:] = rain_half_3km_3hr
fout.variables['rain_half_prob_3km_6hr'][:] = rain_half_3km_6hr

fout.variables['rain_half_prob_15km'][:] = rain_half_15km
fout.variables['rain_half_prob_15km_hourly'][:] = rain_half_15km_hourly
fout.variables['rain_half_prob_15km_3hr'][:] = rain_half_15km_3hr
fout.variables['rain_half_prob_15km_6hr'][:] = rain_half_15km_6hr

fout.variables['rain_half_prob_27km'][:] = rain_half_27km
fout.variables['rain_half_prob_27km_hourly'][:] = rain_half_27km_hourly
fout.variables['rain_half_prob_27km_3hr'][:] = rain_half_27km_3hr
fout.variables['rain_half_prob_27km_6hr'][:] = rain_half_27km_6hr

fout.variables['rain_one_prob_3km'][:] = rain_one_3km
fout.variables['rain_one_prob_3km_hourly'][:] = rain_one_3km_hourly
fout.variables['rain_one_prob_3km_3hr'][:] = rain_one_3km_3hr
fout.variables['rain_one_prob_3km_6hr'][:] = rain_one_3km_6hr

fout.variables['rain_one_prob_15km'][:] = rain_one_15km
fout.variables['rain_one_prob_15km_hourly'][:] = rain_one_15km_hourly
fout.variables['rain_one_prob_15km_3hr'][:] = rain_one_15km_3hr
fout.variables['rain_one_prob_15km_6hr'][:] = rain_one_15km_6hr

fout.variables['rain_one_prob_27km'][:] = rain_one_27km
fout.variables['rain_one_prob_27km_hourly'][:] = rain_one_27km_hourly
fout.variables['rain_one_prob_27km_3hr'][:] = rain_one_27km_3hr
fout.variables['rain_one_prob_27km_6hr'][:] = rain_one_27km_6hr

fout.variables['rain_two_prob_3km'][:] = rain_two_3km
fout.variables['rain_two_prob_3km_hourly'][:] = rain_two_3km_hourly
fout.variables['rain_two_prob_3km_3hr'][:] = rain_two_3km_3hr
fout.variables['rain_two_prob_3km_6hr'][:] = rain_two_3km_6hr

fout.variables['rain_two_prob_15km'][:] = rain_two_15km
fout.variables['rain_two_prob_15km_hourly'][:] = rain_two_15km_hourly
fout.variables['rain_two_prob_15km_3hr'][:] = rain_two_15km_3hr
fout.variables['rain_two_prob_15km_6hr'][:] = rain_two_15km_6hr

fout.variables['rain_two_prob_27km'][:] = rain_two_27km
fout.variables['rain_two_prob_27km_hourly'][:] = rain_two_27km_hourly
fout.variables['rain_two_prob_27km_3hr'][:] = rain_two_27km_3hr
fout.variables['rain_two_prob_27km_6hr'][:] = rain_two_27km_6hr

fout.variables['comp_dz_prob_3km'][:] = comp_dz_3km
fout.variables['comp_dz_prob_3km_hourly'][:] = comp_dz_3km_hourly
fout.variables['comp_dz_prob_15km'][:] = comp_dz_15km
fout.variables['comp_dz_prob_15km_hourly'][:] = comp_dz_15km_hourly
fout.variables['comp_dz_prob_27km'][:] = comp_dz_27km
fout.variables['comp_dz_prob_27km_hourly'][:] = comp_dz_27km_hourly

### Close output file ### 

fout.close()
del fout
