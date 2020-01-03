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

print 'SWATH STARTED: ', t

############################ Set Thresholds for convolution and neighborhood probabilities: #################################

radius_max_9km             = 3			#grid point radius for maximum value filter (3x3 square neighborhood)
radius_max_15km            = 5			#grid point radius for maximum value filter (9x9 square neighborhood)
radius_max_27km            = 9			#grid point radius for maximum value filter (14x14 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

comp_dz_thresh             = 45.		#40 dBZ
rain_1_thresh              = 0.5		#0.5 inches
rain_2_thresh              = 1.			#1 inch
rain_3_thresh              = 2.			#2 inches
soil_moisture_thresh       = .95		#95% soil moisture
hail_thresh                = 1. 		#1 inch
ws_80_thresh               = 50.		#58 mph
w_up_thresh                = 10.                #10 m/s
wz_0to2_thresh             = 0.003		#0.003 s^-1
uh_0to2_thresh             = 30.		#30 m^2/s^2
uh_2to5_thresh             = 60.		#60 m^2/s^2


############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18
#start_t = 0 

ens_files = []
#swt_files = []
summary_files_temp = os.listdir(summary_dir)

############################ Find latest available hourly SWT file to speed up realtime processing: #################################

#for swt_t in range(1, fcst_nt): #skip first time
#   swt_found = 0
#   str_t = str(swt_t)
#   if (len(str_t) == 1):
#      str_t = '0' + str_t
# 
#   for f, file in enumerate(summary_files_temp): 
#      if ((file[-28:-25] == 'SWT') and (file[-24:-22] == str_t)):       #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
#         swt_found = 1
#
#   if ((swt_found == 1) and ((swt_t % 12) == 0) and (swt_t < fcst_nt)): 
#      start_t = swt_t
#      time.sleep(2)

############################ Find and process ENS and SWT files: #################################
   
for f, file in enumerate(summary_files_temp): 
   if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
      ens_files.append(file) 
#   if (file[-28:-25] == 'SWT'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
#      swt_files.append(file) 
      
ens_files.sort() 
#swt_files.sort() 
#print start_t, 'START'

for tt in range(23, t+1):

#   if (start_t > 0): 
#      swt_file = swt_files[start_t]
#      swt_infile = os.path.join(summary_dir, swt_file)

   ens_file = ens_files[tt]
   infile = os.path.join(summary_dir, ens_file)

   if (tt == t): 
   ### Set output path ###
      outname = ens_file[0:7] + '4HR' + ens_file[10:]
      output_path = summary_dir + outname

   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (tt == 23):

#      if (start_t > 0): 
#         try:                                                 #open WRFOUT file
#            swt_fin = netCDF4.Dataset(swt_infile, "r")
#            print "Opening %s \n" % swt_infile
#         except:
#            print "%s does not exist! \n" %swt_infile
#            sys.exit(1)

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

#      if (start_t == 0): 
      uh_2to5                     = fin.variables["uh_2to5"][:,:,:]
      uh_0to2                     = fin.variables["uh_0to2"][:,:,:]
      wz_0to2                     = fin.variables["wz_0to2"][:,:,:]
      rain                        = fin.variables["rain"][:,:,:]
      comp_dz                     = fin.variables["comp_dz"][:,:,:]
      w_up                        = fin.variables["w_up"][:,:,:]
      ws_80                       = fin.variables["ws_80"][:,:,:]
      hail                        = fin.variables["hail"][:,:,:]
      hailcast                    = fin.variables["hailcast"][:,:,:]
#      else: 
#         uh_2to5                            = swt_fin.variables["uh_2to5_out"][:,:,:]
#         uh_0to2                            = swt_fin.variables["uh_0to2_out"][:,:,:]
#         wz_0to2                            = swt_fin.variables["wz_0to2_out"][:,:,:]
#         rain                               = swt_fin.variables["rain_out"][:,:,:]
#         comp_dz                            = swt_fin.variables["comp_dz_out"][:,:,:]
#         w_up                               = swt_fin.variables["w_up_out"][:,:,:]
#         ws_80                              = swt_fin.variables["ws_80_out"][:,:,:]
#         hail                               = swt_fin.variables["hail_out"][:,:,:]
#         hailcast                           = swt_fin.variables["hailcast_out"][:,:,:]

### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

#      if (start_t > 0):
#         swt_fin.close()
#         del swt_fin

   else: 

      uh_2to5                            = np.where((fin.variables["uh_2to5"][:,:,:] > uh_2to5), fin.variables["uh_2to5"][:,:,:], uh_2to5)

###

      uh_0to2                            = np.where((fin.variables["uh_0to2"][:,:,:] > uh_0to2), fin.variables["uh_0to2"][:,:,:], uh_0to2)

###

      wz_0to2                            = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2), fin.variables["wz_0to2"][:,:,:], wz_0to2)

###

      comp_dz_temp                  = fin.variables["comp_dz"][:,:,:]

### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz_temp > 100000., 0., comp_dz_temp)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

      comp_dz                            = np.where((comp_dz_temp > comp_dz), comp_dz_temp, comp_dz)

###

      rain_temp                          = fin.variables["rain"][:,:,:]
      rain                               = rain + rain_temp

###

      w_up                               = np.where((fin.variables["w_up"][:,:,:] > w_up), fin.variables["w_up"][:,:,:], w_up)

###

      ws_80                              = np.where((fin.variables["ws_80"][:,:,:] > ws_80), fin.variables["ws_80"][:,:,:], ws_80)

###

      hail                               = np.where((fin.variables["hail"][:,:,:] > hail), fin.variables["hail"][:,:,:], hail)

###

      hailcast                           = np.where((fin.variables["hailcast"][:,:,:] > hailcast), fin.variables["hailcast"][:,:,:], hailcast)

   fin.close()
   del fin

### Output swaths at hourly intervals to speed up processing: 

#if ((t > 1) and ((t % 12) == 0)):
#   uh_2to5_out = uh_2to5         
#   uh_0to2_out = uh_0to2         
#   wz_0to2_out = wz_0to2         
#   rain_out = rain         
#   comp_dz_out = comp_dz         
#   w_up_out = w_up         
#   ws_80_out = ws_80         
#   hail_out = hail         
#   hailcast_out = hailcast         

##################### Calculate ensemble output variables: ########################

uh_2to5_90, uh_2to5_max, uh_2to5_9km, uh_2to5_15km, uh_2to5_27km = calc_ens_products(uh_2to5, radius_max_9km, radius_max_15km, radius_max_27km, kernel, uh_2to5_thresh)

uh_0to2_90, uh_0to2_max, uh_0to2_9km, uh_0to2_15km, uh_0to2_27km = calc_ens_products(uh_0to2, radius_max_9km, radius_max_15km, radius_max_27km, kernel, uh_0to2_thresh)

wz_0to2_90, wz_0to2_max, wz_0to2_9km, wz_0to2_15km, wz_0to2_27km = calc_ens_products(wz_0to2, radius_max_9km, radius_max_15km, radius_max_27km, kernel, wz_0to2_thresh)

comp_dz_90, comp_dz_max, comp_dz_3km, comp_dz_15km, comp_dz_27km = calc_ens_products(comp_dz, 1, radius_max_15km, radius_max_27km, kernel, comp_dz_thresh)

rain_half_90, rain_half_max, rain_half_3km, rain_half_15km, rain_half_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_1_thresh)

rain_one_90, rain_one_max, rain_one_3km, rain_one_15km, rain_one_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_2_thresh)

rain_two_90, rain_two_max, rain_two_3km, rain_two_15km, rain_two_27km = calc_ens_products(rain, 1, radius_max_15km, radius_max_27km, kernel, rain_3_thresh)

w_up_90, w_up_max, w_up_9km, w_up_15km, w_up_27km = calc_ens_products(w_up, radius_max_9km, radius_max_15km, radius_max_27km, kernel, w_up_thresh)

ws_80_90, ws_80_max, ws_80_9km, ws_80_15km, ws_80_27km = calc_ens_products(ws_80, radius_max_9km, radius_max_15km, radius_max_27km, kernel, ws_80_thresh)

hail_90, hail_max, hail_9km, hail_15km, hail_27km = calc_ens_products(hail, radius_max_9km, radius_max_15km, radius_max_27km, kernel, hail_thresh)

hailcast_90, hailcast_max, hailcast_9km, hailcast_15km, hailcast_27km = calc_ens_products(hailcast, radius_max_9km, radius_max_15km, radius_max_27km, kernel, hail_thresh)


##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)

#handle hourly swath files: 

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

#handle hourly swath files: 

#if ((t > 1) and ((t % 12) == 0)):
#   uh_2to5_hourly_out = fout.createVariable('uh_2to5_out', 'f4', ('NE','NY','NX',))
#   uh_0to2_hourly_out = fout.createVariable('uh_0to2_out', 'f4', ('NE','NY','NX',))
#   wz_0to2_hourly_out = fout.createVariable('wz_0to2_out', 'f4', ('NE','NY','NX',))
#   rain_hourly_out = fout.createVariable('rain_out', 'f4', ('NE','NY','NX',))
#   comp_dz_hourly_out = fout.createVariable('comp_dz_out', 'f4', ('NE','NY','NX',))
#   w_up_hourly_out = fout.createVariable('w_up_out', 'f4', ('NE','NY','NX',))
#   ws_80_hourly_out = fout.createVariable('ws_80_out', 'f4', ('NE','NY','NX',))
#   hail_hourly_out = fout.createVariable('hail_out', 'f4', ('NE','NY','NX',))
#   hailcast_hourly_out = fout.createVariable('hailcast_out', 'f4', ('NE','NY','NX',))

### 90th percentile and max variables ("p" is placeholder to differentiate variable name from one containing data ) ###

uh_2to5_90p = fout.createVariable('uh_2to5_90', 'f4', ('NY','NX',))
uh_2to5_90p.long_name = "Accumulated ensemble 90th percentile value of 2-5 km updraft helicity"
uh_2to5_90p.units = "m^2/s^2"

uh_2to5_maxp = fout.createVariable('uh_2to5_max', 'f4', ('NY','NX',))
uh_2to5_maxp.long_name = "Accumulated ensemble max value of 2-5 km updraft helicity"
uh_2to5_maxp.units = "m^2/s^2"

uh_0to2_90p = fout.createVariable('uh_0to2_90', 'f4', ('NY','NX',))
uh_0to2_90p.long_name = "Accumulated ensemble 90th percentile value of 0-2 km updraft helicity"
uh_0to2_90p.units = "m^2/s^2"

uh_0to2_maxp = fout.createVariable('uh_0to2_max', 'f4', ('NY','NX',))
uh_0to2_maxp.long_name = "Accumulated ensemble max value of 0-2 km updraft helicity"
uh_0to2_maxp.units = "m^2/s^2"

wz_0to2_90p = fout.createVariable('wz_0to2_90', 'f4', ('NY','NX',))
wz_0to2_90p.long_name = "Accumulated ensemble 90th percentile value of average 0-2 km vertical vorticity"
wz_0to2_90p.units = "s^-1"

wz_0to2_maxp = fout.createVariable('wz_0to2_max', 'f4', ('NY','NX',))
wz_0to2_maxp.long_name = "Accumulated ensemble max value of average 0-2 km vertical vorticity"
wz_0to2_maxp.units = "s^-1"

rain_90p = fout.createVariable('rain_90', 'f4', ('NY','NX',))
rain_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_maxp = fout.createVariable('rain_max', 'f4', ('NY','NX',))
rain_maxp.long_name = "Accumulated ensemble max value of accumulated rainfall"
rain_maxp.units = "in"

comp_dz_90p = fout.createVariable('comp_dz_90', 'f4', ('NY','NX',))
comp_dz_90p.long_name = "Accumulated ensemble 90th percentile value of composite reflectivity"
comp_dz_90p.units = "dBZ"

comp_dz_maxp = fout.createVariable('comp_dz_max', 'f4', ('NY','NX',))
comp_dz_maxp.long_name = "Accumulated ensemble max value of composite reflectivity"
comp_dz_maxp.units = "dBZ"

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

w_up_90p = fout.createVariable('w_up_90', 'f4', ('NY','NX',))
w_up_90p.long_name = "Accumulated ensemble 90th percentile value of updraft velocity"
w_up_90p.units = "m/s"

w_up_maxp = fout.createVariable('w_up_max', 'f4', ('NY','NX',))
w_up_maxp.long_name = "Accumulated ensemble max value of updraft velocity"
w_up_maxp.units = "m/s"

ws_80_90p = fout.createVariable('ws_80_90', 'f4', ('NY','NX',))
ws_80_90p.long_name = "Accumulated ensemble 90th percentile value of 80-m wind speed"
ws_80_90p.units = "kts"

ws_80_maxp = fout.createVariable('ws_80_max', 'f4', ('NY','NX',))
ws_80_maxp.long_name = "Accumulated ensemble max value of 80-m wind speed"
ws_80_maxp.units = "kts"

hail_90p = fout.createVariable('hail_90', 'f4', ('NY','NX',))
hail_90p.long_name = "Accumulated ensemble 90th percentile value of max hail size at the surface (NSSL 2-moment)"
hail_90p.units = "in"

hail_maxp = fout.createVariable('hail_max', 'f4', ('NY','NX',))
hail_maxp.long_name = "Accumulated ensemble max value of max hail size at the surface (NSSL 2-moment)"
hail_maxp.units = "in"

hailcast_90p = fout.createVariable('hailcast_90', 'f4', ('NY','NX',))
hailcast_90p.long_name = "Accumulated ensemble 90th percentile value of max hail size at the surface (Hailcast)"
hailcast_90p.units = "in"

hailcast_maxp = fout.createVariable('hailcast_max', 'f4', ('NY','NX',))
hailcast_maxp.long_name = "Accumulated ensemble max value of max hail size at the surface (Hailcast)"
hailcast_maxp.units = "in"

### Probability of exceedence ###

uh_2to5_prob_9 = fout.createVariable('uh_2to5_prob_9km', 'f4', ('NY','NX',))
uh_2to5_prob_9.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (9 km neighborhod)"
uh_2to5_prob_9.units = "%"

uh_2to5_prob_15 = fout.createVariable('uh_2to5_prob_15km', 'f4', ('NY','NX',))
uh_2to5_prob_15.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (15 km neighborhod)"
uh_2to5_prob_15.units = "%"

uh_2to5_prob_27 = fout.createVariable('uh_2to5_prob_27km', 'f4', ('NY','NX',))
uh_2to5_prob_27.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (27 km neighborhod)"
uh_2to5_prob_27.units = "%"

uh_0to2_prob_9 = fout.createVariable('uh_0to2_prob_9km', 'f4', ('NY','NX',))
uh_0to2_prob_9.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (9 km neighborhod)"
uh_0to2_prob_9.units = "%"

uh_0to2_prob_15 = fout.createVariable('uh_0to2_prob_15km', 'f4', ('NY','NX',))
uh_0to2_prob_15.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (15 km neighborhod)"
uh_0to2_prob_15.units = "%"

uh_0to2_prob_27 = fout.createVariable('uh_0to2_prob_27km', 'f4', ('NY','NX',))
uh_0to2_prob_27.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (27 km neighborhod)"
uh_0to2_prob_27.units = "%"

wz_0to2_prob_9 = fout.createVariable('wz_0to2_prob_9km', 'f4', ('NY','NX',))
wz_0to2_prob_9.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (9 km neighborhod)"
wz_0to2_prob_9.units = "%"

wz_0to2_prob_15 = fout.createVariable('wz_0to2_prob_15km', 'f4', ('NY','NX',))
wz_0to2_prob_15.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (15 km neighborhod)"
wz_0to2_prob_15.units = "%"

wz_0to2_prob_27 = fout.createVariable('wz_0to2_prob_27km', 'f4', ('NY','NX',))
wz_0to2_prob_27.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (27 km neighborhod)"
wz_0to2_prob_27.units = "%"

rain_half_prob_3 = fout.createVariable('rain_half_prob_3km', 'f4', ('NY','NX',))
rain_half_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 0.5 inches"
rain_half_prob_3.units = "%"

rain_half_prob_15 = fout.createVariable('rain_half_prob_15km', 'f4', ('NY','NX',))
rain_half_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches (15 km neighborhood)"
rain_half_prob_15.units = "%"

rain_half_prob_27 = fout.createVariable('rain_half_prob_27km', 'f4', ('NY','NX',))
rain_half_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches (27 km neighborhood)"
rain_half_prob_27.units = "%"

rain_one_prob_3 = fout.createVariable('rain_one_prob_3km', 'f4', ('NY','NX',))
rain_one_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 1 inch"
rain_one_prob_3.units = "%"

rain_one_prob_15 = fout.createVariable('rain_one_prob_15km', 'f4', ('NY','NX',))
rain_one_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch (15 km neighborhood)"
rain_one_prob_15.units = "%"

rain_one_prob_27 = fout.createVariable('rain_one_prob_27km', 'f4', ('NY','NX',))
rain_one_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch (27 km neighborhood)"
rain_one_prob_27.units = "%"

rain_two_prob_3 = fout.createVariable('rain_two_prob_3km', 'f4', ('NY','NX',))
rain_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 2 inches"
rain_two_prob_3.units = "%"

rain_two_prob_15 = fout.createVariable('rain_two_prob_15km', 'f4', ('NY','NX',))
rain_two_prob_15.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (15 km neighborhood)"
rain_two_prob_15.units = "%"

rain_two_prob_27 = fout.createVariable('rain_two_prob_27km', 'f4', ('NY','NX',))
rain_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (27 km neighborhood)"
rain_two_prob_27.units = "%"

comp_dz_prob_3 = fout.createVariable('comp_dz_prob_3km', 'f4', ('NY','NX',))
comp_dz_prob_3.long_name = "Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3.units = "%"

comp_dz_prob_15 = fout.createVariable('comp_dz_prob_15km', 'f4', ('NY','NX',))
comp_dz_prob_15.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (15 km neighborhood)"
comp_dz_prob_15.units = "%"

comp_dz_prob_27 = fout.createVariable('comp_dz_prob_27km', 'f4', ('NY','NX',))
comp_dz_prob_27.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27.units = "%"

ws_80_prob_9 = fout.createVariable('ws_80_prob_9km', 'f4', ('NY','NX',))
ws_80_prob_9.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 50 kts (9 km neighborhod)"
ws_80_prob_9.units = "%"

ws_80_prob_15 = fout.createVariable('ws_80_prob_15km', 'f4', ('NY','NX',))
ws_80_prob_15.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 50 kts (15 km neighborhod)"
ws_80_prob_15.units = "%"

ws_80_prob_27 = fout.createVariable('ws_80_prob_27km', 'f4', ('NY','NX',))
ws_80_prob_27.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 50 kts (27 km neighborhod)"
ws_80_prob_27.units = "%"

w_up_prob_9 = fout.createVariable('w_up_prob_9km', 'f4', ('NY','NX',))
w_up_prob_9.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (9 km neighborhod)"
w_up_prob_9.units = "%"

w_up_prob_15 = fout.createVariable('w_up_prob_15km', 'f4', ('NY','NX',))
w_up_prob_15.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (15 km neighborhod)"
w_up_prob_15.units = "%"

w_up_prob_27 = fout.createVariable('w_up_prob_27km', 'f4', ('NY','NX',))
w_up_prob_27.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (27 km neighborhod)"
w_up_prob_27.units = "%"

hail_prob_9 = fout.createVariable('hail_prob_9km', 'f4', ('NY','NX',))
hail_prob_9.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (9 km neighborhod)"
hail_prob_9.units = "%"

hail_prob_15 = fout.createVariable('hail_prob_15km', 'f4', ('NY','NX',))
hail_prob_15.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (15 km neighborhod)"
hail_prob_15.units = "%"

hail_prob_27 = fout.createVariable('hail_prob_27km', 'f4', ('NY','NX',))
hail_prob_27.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (27 km neighborhod)"
hail_prob_27.units = "%"

hailcast_prob_9 = fout.createVariable('hailcast_prob_9km', 'f4', ('NY','NX',))
hailcast_prob_9.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (9 km neighborhod)"
hailcast_prob_9.units = "%"

hailcast_prob_15 = fout.createVariable('hailcast_prob_15km', 'f4', ('NY','NX',))
hailcast_prob_15.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (15 km neighborhod)"
hailcast_prob_15.units = "%"

hailcast_prob_27 = fout.createVariable('hailcast_prob_27km', 'f4', ('NY','NX',))
hailcast_prob_27.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (27 km neighborhod)"
hailcast_prob_27.units = "%"

### Write variables ###

#handle hourly swath files: 

#if ((t > 1) and ((t % 12) == 0)):
#   fout.variables['uh_2to5_out'][:] = uh_2to5_out
#   fout.variables['uh_0to2_out'][:] = uh_0to2_out
#   fout.variables['wz_0to2_out'][:] = wz_0to2_out
#   fout.variables['rain_out'][:] = rain_out
#   fout.variables['comp_dz_out'][:] = comp_dz_out
#   fout.variables['w_up_out'][:] = w_up_out
#   fout.variables['ws_80_out'][:] = ws_80_out
#   fout.variables['hail_out'][:] = hail_out

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['uh_2to5_90'][:] = uh_2to5_90
fout.variables['uh_2to5_max'][:] = uh_2to5_max

fout.variables['uh_0to2_90'][:] = uh_0to2_90
fout.variables['uh_0to2_max'][:] = uh_0to2_max

fout.variables['wz_0to2_90'][:] = wz_0to2_90
fout.variables['wz_0to2_max'][:] = wz_0to2_max

fout.variables['comp_dz_90'][:] = comp_dz_90
fout.variables['comp_dz_max'][:] = comp_dz_max
fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['rain_90'][:] = rain_half_90
fout.variables['rain_max'][:] = rain_half_max

fout.variables['w_up_90'][:] = w_up_90
fout.variables['w_up_max'][:] = w_up_max

fout.variables['ws_80_90'][:] = ws_80_90
fout.variables['ws_80_max'][:] = ws_80_max

fout.variables['hail_90'][:] = hail_90
fout.variables['hail_max'][:] = hail_max

fout.variables['hailcast_90'][:] = hailcast_90
fout.variables['hailcast_max'][:] = hailcast_max

### 

fout.variables['uh_2to5_prob_9km'][:] = uh_2to5_9km
fout.variables['uh_2to5_prob_15km'][:] = uh_2to5_15km
fout.variables['uh_2to5_prob_27km'][:] = uh_2to5_27km
fout.variables['uh_0to2_prob_9km'][:] = uh_0to2_9km
fout.variables['uh_0to2_prob_15km'][:] = uh_0to2_15km
fout.variables['uh_0to2_prob_27km'][:] = uh_0to2_27km
fout.variables['wz_0to2_prob_9km'][:] = wz_0to2_9km
fout.variables['wz_0to2_prob_15km'][:] = wz_0to2_15km
fout.variables['wz_0to2_prob_27km'][:] = wz_0to2_27km
fout.variables['comp_dz_prob_3km'][:] = comp_dz_3km
fout.variables['comp_dz_prob_15km'][:] = comp_dz_15km
fout.variables['comp_dz_prob_27km'][:] = comp_dz_27km
fout.variables['rain_half_prob_3km'][:] = rain_half_3km
fout.variables['rain_half_prob_15km'][:] = rain_half_15km
fout.variables['rain_half_prob_27km'][:] = rain_half_27km
fout.variables['rain_one_prob_3km'][:] = rain_one_3km
fout.variables['rain_one_prob_15km'][:] = rain_one_15km
fout.variables['rain_one_prob_27km'][:] = rain_one_27km
fout.variables['rain_two_prob_3km'][:] = rain_two_3km
fout.variables['rain_two_prob_15km'][:] = rain_two_15km
fout.variables['rain_two_prob_27km'][:] = rain_two_27km
fout.variables['w_up_prob_9km'][:] = w_up_9km
fout.variables['w_up_prob_15km'][:] = w_up_15km
fout.variables['w_up_prob_27km'][:] = w_up_27km
fout.variables['ws_80_prob_9km'][:] = ws_80_9km
fout.variables['ws_80_prob_15km'][:] = ws_80_15km
fout.variables['ws_80_prob_27km'][:] = ws_80_27km
fout.variables['hail_prob_9km'][:] = hail_9km
fout.variables['hail_prob_15km'][:] = hail_15km
fout.variables['hail_prob_27km'][:] = hail_27km
fout.variables['hailcast_prob_9km'][:] = hailcast_9km
fout.variables['hailcast_prob_15km'][:] = hailcast_15km
fout.variables['hailcast_prob_27km'][:] = hailcast_27km

### Close output file ### 

fout.close()
del fout




