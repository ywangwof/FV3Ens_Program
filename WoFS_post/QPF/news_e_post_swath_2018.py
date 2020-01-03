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
radius_max_27km            = 9			#grid point radius for maximum value filter (9x9 square neighborhood)
radius_max_42km            = 14			#grid point radius for maximum value filter (14x14 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

comp_dz_thresh             = 40.		#40 dBZ
rain_1_thresh              = 0.5		#0.5 inches
rain_2_thresh              = 1.			#1 inch
rain_3_thresh              = 2.			#2 inches
soil_moisture_thresh       = .95		#95% soil moisture
hail_thresh                = 1. 		#1 inch
ws_80_thresh               = 50.		#58 kts
w_up_thresh                = 10.                #10 m/s
wz_0to2_thresh             = 0.003		#0.003 s^-1
uh_0to2_thresh             = 30.		#30 m^2/s^2
uh_2to5_thresh             = 60.		#60 m^2/s^2


############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18

############hack to try and fix realtime for the last timestep: 
#if (t == fcst_nt): 
#   print 'ENTERING WHILE LOOP (NEWS-E-POST-SWATH): ', t, fcst_nt
#   lastfile = 0
#   while (lastfile == 0):  
#      ens_files = []
#      summary_files_temp = os.listdir(summary_dir)
#      for f, file in enumerate(summary_files_temp):
#         if (file[-28:-25] == 'ENS'):                               #assumes filename format of: news-e_ENS_20170516_2200_0000.nc
#            ens_files.append(file)
#      ens_files.sort()
#      tempfile = ens_files[-1]
#      if (tempfile[-24:-22] == '36'): 
#         print 'swathswathswathswathswath', tempfile
#         lastfile == 1
#         time.sleep(10) #give it 10 s to finish writing ... just in case
#      else: 
#         print 'Not Yet, SWATH'
#         time.sleep(10) #try waiting 30 s for final ens. file
#############

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
      outname = ens_file[0:7] + 'SWT' + ens_file[10:]
      output_path = summary_dir + outname

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

      uh_2to5                            = fin.variables["uh_2to5"][:,:,:]
      uh_0to2                            = fin.variables["uh_0to2"][:,:,:]
      wz_0to2                            = fin.variables["wz_0to2"][:,:,:]
      rain                               = fin.variables["rain"][:,:,:]
#      soil_moisture                      = fin.variables["soil_moisture"][:,:,:]
      comp_dz                            = fin.variables["comp_dz"][:,:,:]
      w_up                               = fin.variables["w_up"][:,:,:]
      ws_80                              = fin.variables["ws_80"][:,:,:]
      hail                               = fin.variables["hail"][:,:,:]
      hailcast                           = fin.variables["hailcast"][:,:,:]
#      rain_sat                           = np.where(soil_moisture > 0.95, rain, 0.)

      uh_2to5_hourly                     = fin.variables["uh_2to5"][:,:,:]
      uh_0to2_hourly                     = fin.variables["uh_0to2"][:,:,:]
      wz_0to2_hourly                     = fin.variables["wz_0to2"][:,:,:]
      rain_hourly                        = fin.variables["rain"][:,:,:]
#      soil_moisture_hourly               = fin.variables["soil_moisture"][:,:,:]
      comp_dz_hourly                     = fin.variables["comp_dz"][:,:,:]
      w_up_hourly                        = fin.variables["w_up"][:,:,:]
      ws_80_hourly                       = fin.variables["ws_80"][:,:,:]
      hail_hourly                        = fin.variables["hail"][:,:,:]
      hailcast_hourly                    = fin.variables["hailcast"][:,:,:]
#      rain_sat_hourly                    = np.where(soil_moisture_hourly > 0.95, rain_hourly, 0.)

      ### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

   else: 

      uh_2to5                            = np.where((fin.variables["uh_2to5"][:,:,:] > uh_2to5), fin.variables["uh_2to5"][:,:,:], uh_2to5)
      if ((tt > 1) and ((tt % 12) == 1)):
         uh_2to5_hourly                  = fin.variables["uh_2to5"][:,:,:]
      else:  
         uh_2to5_hourly                  = np.where((fin.variables["uh_2to5"][:,:,:] > uh_2to5_hourly), fin.variables["uh_2to5"][:,:,:], uh_2to5_hourly)

###

      uh_0to2                            = np.where((fin.variables["uh_0to2"][:,:,:] > uh_0to2), fin.variables["uh_0to2"][:,:,:], uh_0to2)
      if ((tt > 1) and ((tt % 12) == 1)):
         uh_0to2_hourly                  = fin.variables["uh_0to2"][:,:,:]
      else:
         uh_0to2_hourly                  = np.where((fin.variables["uh_0to2"][:,:,:] > uh_0to2_hourly), fin.variables["uh_0to2"][:,:,:], uh_0to2_hourly)

###

      wz_0to2                            = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2), fin.variables["wz_0to2"][:,:,:], wz_0to2)
      if ((tt > 1) and ((tt % 12) == 1)):
         wz_0to2_hourly                  = fin.variables["wz_0to2"][:,:,:]
      else:
         wz_0to2_hourly                  = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2_hourly), fin.variables["wz_0to2"][:,:,:], wz_0to2_hourly)

###

      comp_dz_temp                  = fin.variables["comp_dz"][:,:,:]

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
      if ((tt > 1) and ((tt % 12) == 1)):
         rain_hourly                     = fin.variables["rain"][:,:,:]
      else:
         rain_hourly                     = rain_hourly + rain_temp

###

#      soil_moisture_temp                 = fin.variables["soil_moisture"][:,:,:]
#      soil_moisture                      = np.where((fin.variables["soil_moisture"][:,:,:] > soil_moisture), fin.variables["soil_moisture"][:,:,:], soil_moisture)
#      if ((tt > 1) and ((tt % 12) == 1)):
#         soil_moisture_hourly            = fin.variables["soil_moisture"][:,:,:]
#      else:
#         soil_moisture_hourly            = np.where((fin.variables["soil_moisture"][:,:,:] > soil_moisture_hourly), fin.variables["soil_moisture"][:,:,:], soil_moisture_hourly)

###

#      rain_sat_temp                          = np.where(soil_moisture_temp > 0.95, rain_temp, 0.)
#      rain_sat                               = rain_sat + rain_sat_temp
#      if ((tt > 1) and ((tt % 12) == 1)):
#         rain_sat_hourly                     = np.where(soil_moisture_temp > 0.95, rain_temp, 0.)
#      else:
#         rain_sat_hourly                     = rain_sat_hourly + rain_sat_temp

###

      w_up                               = np.where((fin.variables["w_up"][:,:,:] > w_up), fin.variables["w_up"][:,:,:], w_up)
      if ((tt > 1) and ((tt % 12) == 1)):
         w_up_hourly                     = fin.variables["w_up"][:,:,:]
      else:
         w_up_hourly                     = np.where((fin.variables["w_up"][:,:,:] > w_up_hourly), fin.variables["w_up"][:,:,:], w_up_hourly)

###

      ws_80                              = np.where((fin.variables["ws_80"][:,:,:] > ws_80), fin.variables["ws_80"][:,:,:], ws_80)
      if ((tt > 1) and ((tt % 12) == 1)):
         ws_80_hourly                    = fin.variables["ws_80"][:,:,:]
      else:
         ws_80_hourly                    = np.where((fin.variables["ws_80"][:,:,:] > ws_80_hourly), fin.variables["ws_80"][:,:,:], ws_80_hourly)

###

      hail                               = np.where((fin.variables["hail"][:,:,:] > hail), fin.variables["hail"][:,:,:], hail)
      if ((tt > 1) and ((tt % 12) == 1)):
         hail_hourly                     = fin.variables["hail"][:,:,:]
      else:
         hail_hourly                     = np.where((fin.variables["hail"][:,:,:] > hail_hourly), fin.variables["hail"][:,:,:], hail_hourly)

###

      hailcast                           = np.where((fin.variables["hailcast"][:,:,:] > hailcast), fin.variables["hailcast"][:,:,:], hailcast)
      if ((tt > 1) and ((tt % 12) == 1)):
         hailcast_hourly                 = fin.variables["hailcast"][:,:,:]
      else:
         hailcast_hourly                 = np.where((fin.variables["hailcast"][:,:,:] > hailcast_hourly), fin.variables["hailcast"][:,:,:], hailcast_hourly)

   fin.close()
   del fin

##################### Calculate ensemble output variables: ########################

uh_2to5_90, uh_2to5_max, uh_2to5_9km, uh_2to5_27km, uh_2to5_42km = calc_ens_products(uh_2to5, radius_max_9km, radius_max_27km, radius_max_42km, kernel, uh_2to5_thresh)
uh_2to5_90_hourly, uh_2to5_max_hourly, uh_2to5_9km_hourly, uh_2to5_27km_hourly, uh_2to5_42km_hourly = calc_ens_products(uh_2to5_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, uh_2to5_thresh)

uh_0to2_90, uh_0to2_max, uh_0to2_9km, uh_0to2_27km, uh_0to2_42km = calc_ens_products(uh_0to2, radius_max_9km, radius_max_27km, radius_max_42km, kernel, uh_0to2_thresh)
uh_0to2_90_hourly, uh_0to2_max_hourly, uh_0to2_9km_hourly, uh_0to2_27km_hourly, uh_0to2_42km_hourly = calc_ens_products(uh_0to2_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, uh_0to2_thresh)

wz_0to2_90, wz_0to2_max, wz_0to2_9km, wz_0to2_27km, wz_0to2_42km = calc_ens_products(wz_0to2, radius_max_9km, radius_max_27km, radius_max_42km, kernel, wz_0to2_thresh)
wz_0to2_90_hourly, wz_0to2_max_hourly, wz_0to2_9km_hourly, wz_0to2_27km_hourly, wz_0to2_42km_hourly = calc_ens_products(wz_0to2_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, wz_0to2_thresh)

comp_dz_90, comp_dz_max, comp_dz_3km, comp_dz_27km, comp_dz_42km = calc_ens_products(comp_dz, 1, radius_max_27km, radius_max_42km, kernel, comp_dz_thresh)
comp_dz_90_hourly, comp_dz_max_hourly, comp_dz_3km_hourly, comp_dz_27km_hourly, comp_dz_42km_hourly = calc_ens_products(comp_dz_hourly, 1, radius_max_27km, radius_max_42km, kernel, comp_dz_thresh)

#soil_moisture_90, soil_moisture_max, soil_moisture_3km, soil_moisture_27km, soil_moisture_42km = calc_ens_products(soil_moisture, 1, radius_max_27km, radius_max_42km, kernel, soil_moisture_thresh)
#soil_moisture_90_hourly, soil_moisture_max_hourly, soil_moisture_3km_hourly, soil_moisture_27km_hourly, soil_moisture_42km_hourly = calc_ens_products(soil_moisture_hourly, 1, radius_max_27km, radius_max_42km, kernel, soil_moisture_thresh)

rain_half_90, rain_half_max, rain_half_3km, rain_half_27km, rain_half_42km = calc_ens_products(rain, 1, radius_max_27km, radius_max_42km, kernel, rain_1_thresh)
rain_half_90_hourly, rain_half_max_hourly, rain_half_3km_hourly, rain_half_27km_hourly, rain_half_42km_hourly = calc_ens_products(rain_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_1_thresh)

rain_one_90, rain_one_max, rain_one_3km, rain_one_27km, rain_one_42km = calc_ens_products(rain, 1, radius_max_27km, radius_max_42km, kernel, rain_2_thresh)
rain_one_90_hourly, rain_one_max_hourly, rain_one_3km_hourly, rain_one_27km_hourly, rain_one_42km_hourly = calc_ens_products(rain_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_2_thresh)

rain_two_90, rain_two_max, rain_two_3km, rain_two_27km, rain_two_42km = calc_ens_products(rain, 1, radius_max_27km, radius_max_42km, kernel, rain_3_thresh)
rain_two_90_hourly, rain_two_max_hourly, rain_two_3km_hourly, rain_two_27km_hourly, rain_two_42km_hourly = calc_ens_products(rain_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_3_thresh)

#rain_sat_half_90, rain_sat_half_max, rain_sat_half_3km, rain_sat_half_27km, rain_sat_half_42km = calc_ens_products(rain_sat, 1, radius_max_27km, radius_max_42km, kernel, rain_1_thresh)
#rain_sat_half_90_hourly, rain_sat_half_max_hourly, rain_sat_half_3km_hourly, rain_sat_half_27km_hourly, rain_sat_half_42km_hourly = calc_ens_products(rain_sat_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_1_thresh)

#rain_sat_one_90, rain_sat_one_max, rain_sat_one_3km, rain_sat_one_27km, rain_sat_one_42km = calc_ens_products(rain_sat, 1, radius_max_27km, radius_max_42km, kernel, rain_2_thresh)
#rain_sat_one_90_hourly, rain_sat_one_max_hourly, rain_sat_one_3km_hourly, rain_sat_one_27km_hourly, rain_sat_one_42km_hourly = calc_ens_products(rain_sat_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_2_thresh)

#rain_sat_two_90, rain_sat_two_max, rain_sat_two_3km, rain_sat_two_27km, rain_sat_two_42km = calc_ens_products(rain_sat, 1, radius_max_27km, radius_max_42km, kernel, rain_3_thresh)
#rain_sat_two_90_hourly, rain_sat_two_max_hourly, rain_sat_two_3km_hourly, rain_sat_two_27km_hourly, rain_sat_two_42km_hourly = calc_ens_products(rain_sat_hourly, 1, radius_max_27km, radius_max_42km, kernel, rain_3_thresh)

w_up_90, w_up_max, w_up_9km, w_up_27km, w_up_42km = calc_ens_products(w_up, radius_max_9km, radius_max_27km, radius_max_42km, kernel, w_up_thresh)
w_up_90_hourly, w_up_max_hourly, w_up_9km_hourly, w_up_27km_hourly, w_up_42km_hourly = calc_ens_products(w_up_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, w_up_thresh)

ws_80_90, ws_80_max, ws_80_9km, ws_80_27km, ws_80_42km = calc_ens_products(ws_80, radius_max_9km, radius_max_27km, radius_max_42km, kernel, ws_80_thresh)
ws_80_90_hourly, ws_80_max_hourly, ws_80_9km_hourly, ws_80_27km_hourly, ws_80_42km_hourly = calc_ens_products(ws_80_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, ws_80_thresh)

hail_90, hail_max, hail_9km, hail_27km, hail_42km = calc_ens_products(hail, radius_max_9km, radius_max_27km, radius_max_42km, kernel, hail_thresh)
hail_90_hourly, hail_max_hourly, hail_9km_hourly, hail_27km_hourly, hail_42km_hourly = calc_ens_products(hail_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, hail_thresh)

hailcast_90, hailcast_max, hailcast_9km, hailcast_27km, hailcast_42km = calc_ens_products(hailcast, radius_max_9km, radius_max_27km, radius_max_42km, kernel, hail_thresh)
hailcast_90_hourly, hailcast_max_hourly, hailcast_9km_hourly, hailcast_27km_hourly, hailcast_42km_hourly = calc_ens_products(hailcast_hourly, radius_max_9km, radius_max_27km, radius_max_42km, kernel, hail_thresh)


##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

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

### 90th percentile and max variables ("p" is placeholder to differentiate variable name from one containing data ) ###

uh_2to5_90p = fout.createVariable('uh_2to5_90', 'f4', ('NY','NX',))
uh_2to5_90p.long_name = "Accumulated ensemble 90th percentile value of 2-5 km updraft helicity"
uh_2to5_90p.units = "m^2/s^2"

uh_2to5_90p_hourly = fout.createVariable('uh_2to5_90_hourly', 'f4', ('NY','NX',))
uh_2to5_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of 2-5 km updraft helicity"
uh_2to5_90p_hourly.units = "m^2/s^2"

uh_2to5_maxp = fout.createVariable('uh_2to5_max', 'f4', ('NY','NX',))
uh_2to5_maxp.long_name = "Accumulated ensemble max value of 2-5 km updraft helicity"
uh_2to5_maxp.units = "m^2/s^2"

uh_2to5_maxp_hourly = fout.createVariable('uh_2to5_max_hourly', 'f4', ('NY','NX',))
uh_2to5_maxp_hourly.long_name = "1-hr Accumulated max percentile value of 2-5 km updraft helicity"
uh_2to5_maxp_hourly.units = "m^2/s^2"

uh_0to2_90p = fout.createVariable('uh_0to2_90', 'f4', ('NY','NX',))
uh_0to2_90p.long_name = "Accumulated ensemble 90th percentile value of 0-2 km updraft helicity"
uh_0to2_90p.units = "m^2/s^2"

uh_0to2_90p_hourly = fout.createVariable('uh_0to2_90_hourly', 'f4', ('NY','NX',))
uh_0to2_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of 0-2 km updraft helicity"
uh_0to2_90p_hourly.units = "m^2/s^2"

uh_0to2_maxp = fout.createVariable('uh_0to2_max', 'f4', ('NY','NX',))
uh_0to2_maxp.long_name = "Accumulated ensemble max value of 0-2 km updraft helicity"
uh_0to2_maxp.units = "m^2/s^2"

uh_0to2_maxp_hourly = fout.createVariable('uh_0to2_max_hourly', 'f4', ('NY','NX',))
uh_0to2_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of 0-2 km updraft helicity"
uh_0to2_maxp_hourly.units = "m^2/s^2"

wz_0to2_90p = fout.createVariable('wz_0to2_90', 'f4', ('NY','NX',))
wz_0to2_90p.long_name = "Accumulated ensemble 90th percentile value of average 0-2 km vertical vorticity"
wz_0to2_90p.units = "s^-1"

wz_0to2_90p_hourly = fout.createVariable('wz_0to2_90_hourly', 'f4', ('NY','NX',))
wz_0to2_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of average 0-2 km vertical vorticity"
wz_0to2_90p_hourly.units = "s^-1"

wz_0to2_maxp = fout.createVariable('wz_0to2_max', 'f4', ('NY','NX',))
wz_0to2_maxp.long_name = "Accumulated ensemble max value of average 0-2 km vertical vorticity"
wz_0to2_maxp.units = "s^-1"

wz_0to2_maxp_hourly = fout.createVariable('wz_0to2_max_hourly', 'f4', ('NY','NX',))
wz_0to2_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of average 0-2 km vertical vorticity"
wz_0to2_maxp_hourly.units = "s^-1"

rain_90p = fout.createVariable('rain_90', 'f4', ('NY','NX',))
rain_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p.units = "in"

rain_90p_hourly = fout.createVariable('rain_90_hourly', 'f4', ('NY','NX',))
rain_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of accumulated rainfall"
rain_90p_hourly.units = "in"

rain_maxp = fout.createVariable('rain_max', 'f4', ('NY','NX',))
rain_maxp.long_name = "Accumulated ensemble max value of accumulated rainfall"
rain_maxp.units = "in"

rain_maxp_hourly = fout.createVariable('rain_max_hourly', 'f4', ('NY','NX',))
rain_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of accumulated rainfall"
rain_maxp_hourly.units = "in"

#rain_sat_90p = fout.createVariable('rain_sat_90', 'f4', ('NY','NX',))
#rain_sat_90p.long_name = "Accumulated ensemble 90th percentile value of accumulated rainfall on saturated soil"
#rain_sat_90p.units = "in"

#rain_sat_90p_hourly = fout.createVariable('rain_sat_90_hourly', 'f4', ('NY','NX',))
#rain_sat_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of accumulated rainfall on saturated soil"
#rain_sat_90p_hourly.units = "in"

#rain_sat_maxp = fout.createVariable('rain_sat_max', 'f4', ('NY','NX',))
#rain_sat_maxp.long_name = "Accumulated ensemble max value of accumulated rainfall on saturated soil"
#rain_sat_maxp.units = "in"

#rain_sat_maxp_hourly = fout.createVariable('rain_sat_max_hourly', 'f4', ('NY','NX',))
#rain_sat_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of accumulated rainfall on saturated soil"
#rain_sat_maxp_hourly.units = "in"

#soil_moisture_90p = fout.createVariable('soil_moisture_90', 'f4', ('NY','NX',))
#soil_moisture_90p.long_name = "Accumulated ensemble 90th percentile value of soil_moisture"
#soil_moisture_90p.units = "%"

#soil_moisture_90p_hourly = fout.createVariable('soil_moisture_90_hourly', 'f4', ('NY','NX',))
#soil_moisture_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of soil_moisture"
#soil_moisture_90p_hourly.units = "%"

#soil_moisture_maxp = fout.createVariable('soil_moisture_max', 'f4', ('NY','NX',))
#soil_moisture_maxp.long_name = "Accumulated ensemble max value of soil_moisture"
#soil_moisture_maxp.units = "%"

#soil_moisture_maxp_hourly = fout.createVariable('soil_moisture_max_hourly', 'f4', ('NY','NX',))
#soil_moisture_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of soil_moisture"
#soil_moisture_maxp_hourly.units = "%"

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

w_up_90p = fout.createVariable('w_up_90', 'f4', ('NY','NX',))
w_up_90p.long_name = "Accumulated ensemble 90th percentile value of updraft velocity"
w_up_90p.units = "m/s"

w_up_90p_hourly = fout.createVariable('w_up_90_hourly', 'f4', ('NY','NX',))
w_up_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of updraft velocity"
w_up_90p_hourly.units = "m/s"

w_up_maxp = fout.createVariable('w_up_max', 'f4', ('NY','NX',))
w_up_maxp.long_name = "Accumulated ensemble max value of updraft velocity"
w_up_maxp.units = "m/s"

w_up_maxp_hourly = fout.createVariable('w_up_max_hourly', 'f4', ('NY','NX',))
w_up_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of updraft velocity"
w_up_maxp_hourly.units = "m/s"

ws_80_90p = fout.createVariable('ws_80_90', 'f4', ('NY','NX',))
ws_80_90p.long_name = "Accumulated ensemble 90th percentile value of 80-m wind speed"
ws_80_90p.units = "kts"

ws_80_90p_hourly = fout.createVariable('ws_80_90_hourly', 'f4', ('NY','NX',))
ws_80_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of 80-m wind speed"
ws_80_90p_hourly.units = "kts"

ws_80_maxp = fout.createVariable('ws_80_max', 'f4', ('NY','NX',))
ws_80_maxp.long_name = "Accumulated ensemble max value of 80-m wind speed"
ws_80_maxp.units = "kts"

ws_80_maxp_hourly = fout.createVariable('ws_80_max_hourly', 'f4', ('NY','NX',))
ws_80_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of 80-m wind speed"
ws_80_maxp_hourly.units = "kts"

hail_90p = fout.createVariable('hail_90', 'f4', ('NY','NX',))
hail_90p.long_name = "Accumulated ensemble 90th percentile value of max hail size at the surface (NSSL 2-moment)"
hail_90p.units = "in"

hail_90p_hourly = fout.createVariable('hail_90_hourly', 'f4', ('NY','NX',))
hail_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of max hail size at the surface (NSSL 2-moment)"
hail_90p_hourly.units = "in"

hail_maxp = fout.createVariable('hail_max', 'f4', ('NY','NX',))
hail_maxp.long_name = "Accumulated ensemble max value of max hail size at the surface (NSSL 2-moment)"
hail_maxp.units = "in"

hail_maxp_hourly = fout.createVariable('hail_max_hourly', 'f4', ('NY','NX',))
hail_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of max hail size at the surface (NSSL 2-moment)"
hail_maxp_hourly.units = "in"

hailcast_90p = fout.createVariable('hailcast_90', 'f4', ('NY','NX',))
hailcast_90p.long_name = "Accumulated ensemble 90th percentile value of max hail size at the surface (Hailcast)"
hailcast_90p.units = "in"

hailcast_90p_hourly = fout.createVariable('hailcast_90_hourly', 'f4', ('NY','NX',))
hailcast_90p_hourly.long_name = "1-hr Accumulated ensemble 90th percentile value of max hail size at the surface (Hailcast)"
hailcast_90p_hourly.units = "in"

hailcast_maxp = fout.createVariable('hailcast_max', 'f4', ('NY','NX',))
hailcast_maxp.long_name = "Accumulated ensemble max value of max hail size at the surface (Hailcast)"
hailcast_maxp.units = "in"

hailcast_maxp_hourly = fout.createVariable('hailcast_max_hourly', 'f4', ('NY','NX',))
hailcast_maxp_hourly.long_name = "1-hr Accumulated ensemble max value of max hail size at the surface (Hailcast)"
hailcast_maxp_hourly.units = "in"

### Probability of exceedence ###

uh_2to5_prob_9 = fout.createVariable('uh_2to5_prob_9km', 'f4', ('NY','NX',))
uh_2to5_prob_9.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (9 km neighborhod)"
uh_2to5_prob_9.units = "%"

uh_2to5_prob_9_hourly = fout.createVariable('uh_2to5_prob_9km_hourly', 'f4', ('NY','NX',))
uh_2to5_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (9 km neighborhod)"
uh_2to5_prob_9_hourly.units = "%"

uh_2to5_prob_27 = fout.createVariable('uh_2to5_prob_27km', 'f4', ('NY','NX',))
uh_2to5_prob_27.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (27 km neighborhod)"
uh_2to5_prob_27.units = "%"

uh_2to5_prob_27_hourly = fout.createVariable('uh_2to5_prob_27km_hourly', 'f4', ('NY','NX',))
uh_2to5_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (27 km neighborhod)"
uh_2to5_prob_27_hourly.units = "%"

uh_2to5_prob_42 = fout.createVariable('uh_2to5_prob_42km', 'f4', ('NY','NX',))
uh_2to5_prob_42.long_name = "Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (42 km neighborhod)"
uh_2to5_prob_42.units = "%"

uh_2to5_prob_42_hourly = fout.createVariable('uh_2to5_prob_42km_hourly', 'f4', ('NY','NX',))
uh_2to5_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of 2-5 km updraft helicity greater than 60 m^2/S^2 (42 km neighborhod)"
uh_2to5_prob_42_hourly.units = "%"

uh_0to2_prob_9 = fout.createVariable('uh_0to2_prob_9km', 'f4', ('NY','NX',))
uh_0to2_prob_9.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (9 km neighborhod)"
uh_0to2_prob_9.units = "%"

uh_0to2_prob_9_hourly = fout.createVariable('uh_0to2_prob_9km_hourly', 'f4', ('NY','NX',))
uh_0to2_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (9 km neighborhod)"
uh_0to2_prob_9_hourly.units = "%"

uh_0to2_prob_27 = fout.createVariable('uh_0to2_prob_27km', 'f4', ('NY','NX',))
uh_0to2_prob_27.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (27 km neighborhod)"
uh_0to2_prob_27.units = "%"

uh_0to2_prob_27_hourly = fout.createVariable('uh_0to2_prob_27km_hourly', 'f4', ('NY','NX',))
uh_0to2_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (27 km neighborhod)"
uh_0to2_prob_27_hourly.units = "%"

uh_0to2_prob_42 = fout.createVariable('uh_0to2_prob_42km', 'f4', ('NY','NX',))
uh_0to2_prob_42.long_name = "Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (42 km neighborhod)"
uh_0to2_prob_42.units = "%"

uh_0to2_prob_42_hourly = fout.createVariable('uh_0to2_prob_42km_hourly', 'f4', ('NY','NX',))
uh_0to2_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km updraft helicity greater than 30 m^2/S^2 (42 km neighborhod)"
uh_0to2_prob_42_hourly.units = "%"

wz_0to2_prob_9 = fout.createVariable('wz_0to2_prob_9km', 'f4', ('NY','NX',))
wz_0to2_prob_9.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (9 km neighborhod)"
wz_0to2_prob_9.units = "%"

wz_0to2_prob_9_hourly = fout.createVariable('wz_0to2_prob_9km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (9 km neighborhod)"
wz_0to2_prob_9_hourly.units = "%"

wz_0to2_prob_27 = fout.createVariable('wz_0to2_prob_27km', 'f4', ('NY','NX',))
wz_0to2_prob_27.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (27 km neighborhod)"
wz_0to2_prob_27.units = "%"

wz_0to2_prob_27_hourly = fout.createVariable('wz_0to2_prob_27km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (27 km neighborhod)"
wz_0to2_prob_27_hourly.units = "%"

wz_0to2_prob_42 = fout.createVariable('wz_0to2_prob_42km', 'f4', ('NY','NX',))
wz_0to2_prob_42.long_name = "Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (42 km neighborhod)"
wz_0to2_prob_42.units = "%"

wz_0to2_prob_42_hourly = fout.createVariable('wz_0to2_prob_42km_hourly', 'f4', ('NY','NX',))
wz_0to2_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of 0-2 km average vertical vorticity greater than 0.003 s^-1 (42 km neighborhod)"
wz_0to2_prob_42_hourly.units = "%"

rain_half_prob_3 = fout.createVariable('rain_half_prob_3km', 'f4', ('NY','NX',))
rain_half_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 0.5 inches"
rain_half_prob_3.units = "%"

rain_half_prob_3_hourly = fout.createVariable('rain_half_prob_3km_hourly', 'f4', ('NY','NX',))
rain_half_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 0.5 inches"
rain_half_prob_3_hourly.units = "%"

rain_half_prob_27 = fout.createVariable('rain_half_prob_27km', 'f4', ('NY','NX',))
rain_half_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches (27 km neighborhood)"
rain_half_prob_27.units = "%"

rain_half_prob_27_hourly = fout.createVariable('rain_half_prob_27km_hourly', 'f4', ('NY','NX',))
rain_half_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 0.5 inches (27 km neighborhood)"
rain_half_prob_27_hourly.units = "%"

rain_half_prob_42 = fout.createVariable('rain_half_prob_42km', 'f4', ('NY','NX',))
rain_half_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches (42 km neighborhood)"
rain_half_prob_42.units = "%"

rain_half_prob_42_hourly = fout.createVariable('rain_half_prob_42km_hourly', 'f4', ('NY','NX',))
rain_half_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 0.5 inches (42 km neighborhood)"
rain_half_prob_42_hourly.units = "%"

rain_one_prob_3 = fout.createVariable('rain_one_prob_3km', 'f4', ('NY','NX',))
rain_one_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 1 inch"
rain_one_prob_3.units = "%"

rain_one_prob_3_hourly = fout.createVariable('rain_one_prob_3km_hourly', 'f4', ('NY','NX',))
rain_one_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 1 inch"
rain_one_prob_3_hourly.units = "%"

rain_one_prob_27 = fout.createVariable('rain_one_prob_27km', 'f4', ('NY','NX',))
rain_one_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch (27 km neighborhood)"
rain_one_prob_27.units = "%"

rain_one_prob_27_hourly = fout.createVariable('rain_one_prob_27km_hourly', 'f4', ('NY','NX',))
rain_one_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 1 inch (27 km neighborhood)"
rain_one_prob_27_hourly.units = "%"

rain_one_prob_42 = fout.createVariable('rain_one_prob_42km', 'f4', ('NY','NX',))
rain_one_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch (42 km neighborhood)"
rain_one_prob_42.units = "%"

rain_one_prob_42_hourly = fout.createVariable('rain_one_prob_42km_hourly', 'f4', ('NY','NX',))
rain_one_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 1 inch (42 km neighborhood)"
rain_one_prob_42_hourly.units = "%"

rain_two_prob_3 = fout.createVariable('rain_two_prob_3km', 'f4', ('NY','NX',))
rain_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 2 inches"
rain_two_prob_3.units = "%"

rain_two_prob_3_hourly = fout.createVariable('rain_two_prob_3km_hourly', 'f4', ('NY','NX',))
rain_two_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 2 inches"
rain_two_prob_3_hourly.units = "%"

rain_two_prob_27 = fout.createVariable('rain_two_prob_27km', 'f4', ('NY','NX',))
rain_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (27 km neighborhood)"
rain_two_prob_27.units = "%"

rain_two_prob_27_hourly = fout.createVariable('rain_two_prob_27km_hourly', 'f4', ('NY','NX',))
rain_two_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches (27 km neighborhood)"
rain_two_prob_27_hourly.units = "%"

rain_two_prob_42 = fout.createVariable('rain_two_prob_42km', 'f4', ('NY','NX',))
rain_two_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches (42 km neighborhood)"
rain_two_prob_42.units = "%"

rain_two_prob_42_hourly = fout.createVariable('rain_two_prob_42km_hourly', 'f4', ('NY','NX',))
rain_two_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches (42 km neighborhood)"
rain_two_prob_42_hourly.units = "%"

#rain_sat_half_prob_3 = fout.createVariable('rain_sat_half_prob_3km', 'f4', ('NY','NX',))
#rain_sat_half_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 0.5 inches on saturated soil"
#rain_sat_half_prob_3.units = "%"

#rain_sat_half_prob_3_hourly = fout.createVariable('rain_sat_half_prob_3km_hourly', 'f4', ('NY','NX',))
#rain_sat_half_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 0.5 inches on saturated soil"
#rain_sat_half_prob_3_hourly.units = "%"

#rain_sat_half_prob_27 = fout.createVariable('rain_sat_half_prob_27km', 'f4', ('NY','NX',))
#rain_sat_half_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches on saturated soil (27 km neighborhood)"
#rain_sat_half_prob_27.units = "%"

#rain_sat_half_prob_27_hourly = fout.createVariable('rain_sat_half_prob_27km_hourly', 'f4', ('NY','NX',))
#rain_sat_half_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 0.5 inches on saturated soil (27 km neighborhood)"
#rain_sat_half_prob_27_hourly.units = "%"

#rain_sat_half_prob_42 = fout.createVariable('rain_sat_half_prob_42km', 'f4', ('NY','NX',))
#rain_sat_half_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 0.5 inches on saturated soil (42 km neighborhood)"
#rain_sat_half_prob_42.units = "%"

#rain_sat_half_prob_42_hourly = fout.createVariable('rain_sat_half_prob_42km_hourly', 'f4', ('NY','NX',))
#rain_sat_half_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 0.5 inches on saturated soil (42 km neighborhood)"
#rain_sat_half_prob_42_hourly.units = "%"

#rain_sat_one_prob_3 = fout.createVariable('rain_sat_one_prob_3km', 'f4', ('NY','NX',))
#rain_sat_one_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 1 inch on saturated soil"
#rain_sat_one_prob_3.units = "%"

#rain_sat_one_prob_3_hourly = fout.createVariable('rain_sat_one_prob_3km_hourly', 'f4', ('NY','NX',))
#rain_sat_one_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 1 inch on saturated soil"
#rain_sat_one_prob_3_hourly.units = "%"

#rain_sat_one_prob_27 = fout.createVariable('rain_sat_one_prob_27km', 'f4', ('NY','NX',))
#rain_sat_one_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch on saturated soil (27 km neighborhood)"
#rain_sat_one_prob_27.units = "%"

#rain_sat_one_prob_27_hourly = fout.createVariable('rain_sat_one_prob_27km_hourly', 'f4', ('NY','NX',))
#rain_sat_one_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 1 inch on saturated soil (27 km neighborhood)"
#rain_sat_one_prob_27_hourly.units = "%"

#rain_sat_one_prob_42 = fout.createVariable('rain_sat_one_prob_42km', 'f4', ('NY','NX',))
#rain_sat_one_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 1 inch on saturated soil (42 km neighborhood)"
#rain_sat_one_prob_42.units = "%"

#rain_sat_one_prob_42_hourly = fout.createVariable('rain_sat_one_prob_42km_hourly', 'f4', ('NY','NX',))
#rain_sat_one_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 1 inch on saturated soil (42 km neighborhood)"
#rain_sat_one_prob_42_hourly.units = "%"

#rain_sat_two_prob_3 = fout.createVariable('rain_sat_two_prob_3km', 'f4', ('NY','NX',))
#rain_sat_two_prob_3.long_name = "Accumulated ensemble gridpoint probability of rainfall greater than 2 inches on saturated soil"
#rain_sat_two_prob_3.units = "%"

#rain_sat_two_prob_3_hourly = fout.createVariable('rain_sat_two_prob_3km_hourly', 'f4', ('NY','NX',))
#rain_sat_two_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of rainfall greater than 2 inches on saturated soil"
#rain_sat_two_prob_3_hourly.units = "%"

#rain_sat_two_prob_27 = fout.createVariable('rain_sat_two_prob_27km', 'f4', ('NY','NX',))
#rain_sat_two_prob_27.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches on saturated soil (27 km neighborhood)"
#rain_sat_two_prob_27.units = "%"

#rain_sat_two_prob_27_hourly = fout.createVariable('rain_sat_two_prob_27km_hourly', 'f4', ('NY','NX',))
#rain_sat_two_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches on saturated soil (27 km neighborhood)"
#rain_sat_two_prob_27_hourly.units = "%"

#rain_sat_two_prob_42 = fout.createVariable('rain_sat_two_prob_42km', 'f4', ('NY','NX',))
#rain_sat_two_prob_42.long_name = "Accumulated ensemble probability of rainfall greater than 2 inches on saturated soil (42 km neighborhood)"
#rain_sat_two_prob_42.units = "%"

#rain_sat_two_prob_42_hourly = fout.createVariable('rain_sat_two_prob_42km_hourly', 'f4', ('NY','NX',))
#rain_sat_two_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of rainfall greater than 2 inches on saturated soil (42 km neighborhood)"
#rain_sat_two_prob_42_hourly.units = "%"

#soil_moisture_prob_3 = fout.createVariable('soil_moisture_prob_3km', 'f4', ('NY','NX',))
#soil_moisture_prob_3.long_name = "Accumulated ensemble gridpoint probability of top layer soil moisture greater than 95%"
#soil_moisture_prob_3.units = "%"

#soil_moisture_prob_3_hourly = fout.createVariable('soil_moisture_prob_3km_hourly', 'f4', ('NY','NX',))
#soil_moisture_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of top layer soil moisture greater than 95%"
#soil_moisture_prob_3_hourly.units = "%"

#soil_moisture_prob_27 = fout.createVariable('soil_moisture_prob_27km', 'f4', ('NY','NX',))
#soil_moisture_prob_27.long_name = "Accumulated ensemble probability of top layer soil moisture greater than 95% (27 km neighborhood)"
#soil_moisture_prob_27.units = "%"

#soil_moisture_prob_27_hourly = fout.createVariable('soil_moisture_prob_27km_hourly', 'f4', ('NY','NX',))
#soil_moisture_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of top layer soil moisture greater than 95% (27 km neighborhood)"
#soil_moisture_prob_27_hourly.units = "%"

#soil_moisture_prob_42 = fout.createVariable('soil_moisture_prob_42km', 'f4', ('NY','NX',))
#soil_moisture_prob_42.long_name = "Accumulated ensemble probability of top layer soil moisture greater than 95% (42 km neighborhood)"
#soil_moisture_prob_42.units = "%"

#soil_moisture_prob_42_hourly = fout.createVariable('soil_moisture_prob_42km_hourly', 'f4', ('NY','NX',))
#soil_moisture_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of top layer soil moisture greater than 95% (42 km neighborhood)"
#soil_moisture_prob_42_hourly.units = "%"

comp_dz_prob_3 = fout.createVariable('comp_dz_prob_3km', 'f4', ('NY','NX',))
comp_dz_prob_3.long_name = "Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3.units = "%"

comp_dz_prob_3_hourly = fout.createVariable('comp_dz_prob_3km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_3_hourly.long_name = "1-hr Accumulated ensemble gridpoint probability of composite reflectivity greater than 40 dBZ"
comp_dz_prob_3_hourly.units = "%"

comp_dz_prob_27 = fout.createVariable('comp_dz_prob_27km', 'f4', ('NY','NX',))
comp_dz_prob_27.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27.units = "%"

comp_dz_prob_27_hourly = fout.createVariable('comp_dz_prob_27km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (27 km neighborhood)"
comp_dz_prob_27_hourly.units = "%"

comp_dz_prob_42 = fout.createVariable('comp_dz_prob_42km', 'f4', ('NY','NX',))
comp_dz_prob_42.long_name = "Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (42 km neighborhood)"
comp_dz_prob_42.units = "%"

comp_dz_prob_42_hourly = fout.createVariable('comp_dz_prob_42km_hourly', 'f4', ('NY','NX',))
comp_dz_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of composite reflectivity greater than 40 dBZ (42 km neighborhood)"
comp_dz_prob_42_hourly.units = "%"

ws_80_prob_9 = fout.createVariable('ws_80_prob_9km', 'f4', ('NY','NX',))
ws_80_prob_9.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 58 kts (9 km neighborhod)"
ws_80_prob_9.units = "%"

ws_80_prob_9_hourly = fout.createVariable('ws_80_prob_9km_hourly', 'f4', ('NY','NX',))
ws_80_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 58 kts (9 km neighborhod)"
ws_80_prob_9_hourly.units = "%"

ws_80_prob_27 = fout.createVariable('ws_80_prob_27km', 'f4', ('NY','NX',))
ws_80_prob_27.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 58 kts (27 km neighborhod)"
ws_80_prob_27.units = "%"

ws_80_prob_27_hourly = fout.createVariable('ws_80_prob_27km_hourly', 'f4', ('NY','NX',))
ws_80_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 58 kts (27 km neighborhod)"
ws_80_prob_27_hourly.units = "%"

ws_80_prob_42 = fout.createVariable('ws_80_prob_42km', 'f4', ('NY','NX',))
ws_80_prob_42.long_name = "Accumulated ensemble probability of 80-m wind speed greater than 58 kts (42 km neighborhod)"
ws_80_prob_42.units = "%"

ws_80_prob_42_hourly = fout.createVariable('ws_80_prob_42km_hourly', 'f4', ('NY','NX',))
ws_80_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of 80-m wind speed greater than 58 kts (42 km neighborhod)"
ws_80_prob_42_hourly.units = "%"

w_up_prob_9 = fout.createVariable('w_up_prob_9km', 'f4', ('NY','NX',))
w_up_prob_9.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (9 km neighborhod)"
w_up_prob_9.units = "%"

w_up_prob_9_hourly = fout.createVariable('w_up_prob_9km_hourly', 'f4', ('NY','NX',))
w_up_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of updraft velocity greater than 10 m/s (9 km neighborhod)"
w_up_prob_9_hourly.units = "%"

w_up_prob_27 = fout.createVariable('w_up_prob_27km', 'f4', ('NY','NX',))
w_up_prob_27.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (27 km neighborhod)"
w_up_prob_27.units = "%"

w_up_prob_27_hourly = fout.createVariable('w_up_prob_27km_hourly', 'f4', ('NY','NX',))
w_up_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of updraft velocity greater than 10 m/s (27 km neighborhod)"
w_up_prob_27_hourly.units = "%"

w_up_prob_42 = fout.createVariable('w_up_prob_42km', 'f4', ('NY','NX',))
w_up_prob_42.long_name = "Accumulated ensemble probability of updraft velocity greater than 10 m/s (42 km neighborhod)"
w_up_prob_42.units = "%"

w_up_prob_42_hourly = fout.createVariable('w_up_prob_42km_hourly', 'f4', ('NY','NX',))
w_up_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of updraft velocity greater than 10 m/s (42 km neighborhod)"
w_up_prob_42_hourly.units = "%"

hail_prob_9 = fout.createVariable('hail_prob_9km', 'f4', ('NY','NX',))
hail_prob_9.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (9 km neighborhod)"
hail_prob_9.units = "%"

hail_prob_9_hourly = fout.createVariable('hail_prob_9km_hourly', 'f4', ('NY','NX',))
hail_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (9 km neighborhod)"
hail_prob_9_hourly.units = "%"

hail_prob_27 = fout.createVariable('hail_prob_27km', 'f4', ('NY','NX',))
hail_prob_27.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (27 km neighborhod)"
hail_prob_27.units = "%"

hail_prob_27_hourly = fout.createVariable('hail_prob_27km_hourly', 'f4', ('NY','NX',))
hail_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (27 km neighborhod)"
hail_prob_27_hourly.units = "%"

hail_prob_42 = fout.createVariable('hail_prob_42km', 'f4', ('NY','NX',))
hail_prob_42.long_name = "Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (42 km neighborhod)"
hail_prob_42.units = "%"

hail_prob_42_hourly = fout.createVariable('hail_prob_42km_hourly', 'f4', ('NY','NX',))
hail_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (NSSL 2-Moment) greater than 1 in (42 km neighborhod)"
hail_prob_42_hourly.units = "%"

hailcast_prob_9 = fout.createVariable('hailcast_prob_9km', 'f4', ('NY','NX',))
hailcast_prob_9.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (9 km neighborhod)"
hailcast_prob_9.units = "%"

hailcast_prob_9_hourly = fout.createVariable('hailcast_prob_9km_hourly', 'f4', ('NY','NX',))
hailcast_prob_9_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (9 km neighborhod)"
hailcast_prob_9_hourly.units = "%"

hailcast_prob_27 = fout.createVariable('hailcast_prob_27km', 'f4', ('NY','NX',))
hailcast_prob_27.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (27 km neighborhod)"
hailcast_prob_27.units = "%"

hailcast_prob_27_hourly = fout.createVariable('hailcast_prob_27km_hourly', 'f4', ('NY','NX',))
hailcast_prob_27_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (27 km neighborhod)"
hailcast_prob_27_hourly.units = "%"

hailcast_prob_42 = fout.createVariable('hailcast_prob_42km', 'f4', ('NY','NX',))
hailcast_prob_42.long_name = "Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (42 km neighborhod)"
hailcast_prob_42.units = "%"

hailcast_prob_42_hourly = fout.createVariable('hailcast_prob_42km_hourly', 'f4', ('NY','NX',))
hailcast_prob_42_hourly.long_name = "1-hr Accumulated ensemble probability of maximum hail size at the surface (Hailcast) greater than 1 in (42 km neighborhod)"
hailcast_prob_42_hourly.units = "%"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['uh_2to5_90'][:] = uh_2to5_90
fout.variables['uh_2to5_90_hourly'][:] = uh_2to5_90_hourly
fout.variables['uh_2to5_max'][:] = uh_2to5_max
fout.variables['uh_2to5_max_hourly'][:] = uh_2to5_max_hourly

fout.variables['uh_0to2_90'][:] = uh_0to2_90
fout.variables['uh_0to2_90_hourly'][:] = uh_0to2_90_hourly
fout.variables['uh_0to2_max'][:] = uh_0to2_max
fout.variables['uh_0to2_max_hourly'][:] = uh_0to2_max_hourly

fout.variables['wz_0to2_90'][:] = wz_0to2_90
fout.variables['wz_0to2_90_hourly'][:] = wz_0to2_90_hourly
fout.variables['wz_0to2_max'][:] = wz_0to2_max
fout.variables['wz_0to2_max_hourly'][:] = wz_0to2_max_hourly

fout.variables['comp_dz_90'][:] = comp_dz_90
fout.variables['comp_dz_90_hourly'][:] = comp_dz_90_hourly
fout.variables['comp_dz_max'][:] = comp_dz_max
fout.variables['comp_dz_max_hourly'][:] = comp_dz_max_hourly
fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['rain_90'][:] = rain_half_90
fout.variables['rain_90_hourly'][:] = rain_half_90_hourly
fout.variables['rain_max'][:] = rain_half_max
fout.variables['rain_max_hourly'][:] = rain_half_max_hourly

#fout.variables['rain_sat_90'][:] = rain_sat_half_90
#fout.variables['rain_sat_90_hourly'][:] = rain_sat_half_90_hourly
#fout.variables['rain_sat_max'][:] = rain_sat_half_max
#fout.variables['rain_sat_max_hourly'][:] = rain_sat_half_max_hourly

#fout.variables['soil_moisture_90'][:] = soil_moisture_90
#fout.variables['soil_moisture_90_hourly'][:] = soil_moisture_90_hourly
#fout.variables['soil_moisture_max'][:] = soil_moisture_max
#fout.variables['soil_moisture_max_hourly'][:] = soil_moisture_max_hourly

fout.variables['w_up_90'][:] = w_up_90
fout.variables['w_up_90_hourly'][:] = w_up_90_hourly
fout.variables['w_up_max'][:] = w_up_max
fout.variables['w_up_max_hourly'][:] = w_up_max_hourly

fout.variables['ws_80_90'][:] = ws_80_90
fout.variables['ws_80_90_hourly'][:] = ws_80_90_hourly
fout.variables['ws_80_max'][:] = ws_80_max
fout.variables['ws_80_max_hourly'][:] = ws_80_max_hourly

fout.variables['hail_90'][:] = hail_90
fout.variables['hail_90_hourly'][:] = hail_90_hourly
fout.variables['hail_max'][:] = hail_max
fout.variables['hail_max_hourly'][:] = hail_max_hourly

fout.variables['hailcast_90'][:] = hailcast_90
fout.variables['hailcast_90_hourly'][:] = hailcast_90_hourly
fout.variables['hailcast_max'][:] = hailcast_max
fout.variables['hailcast_max_hourly'][:] = hailcast_max_hourly

### 

fout.variables['uh_2to5_prob_9km'][:] = uh_2to5_9km
fout.variables['uh_2to5_prob_9km_hourly'][:] = uh_2to5_9km_hourly
fout.variables['uh_2to5_prob_27km'][:] = uh_2to5_27km
fout.variables['uh_2to5_prob_27km_hourly'][:] = uh_2to5_27km_hourly
fout.variables['uh_2to5_prob_42km'][:] = uh_2to5_42km
fout.variables['uh_2to5_prob_42km_hourly'][:] = uh_2to5_42km_hourly
fout.variables['uh_0to2_prob_9km'][:] = uh_0to2_9km
fout.variables['uh_0to2_prob_9km_hourly'][:] = uh_0to2_9km_hourly
fout.variables['uh_0to2_prob_27km'][:] = uh_0to2_27km
fout.variables['uh_0to2_prob_27km_hourly'][:] = uh_0to2_27km_hourly
fout.variables['uh_0to2_prob_42km'][:] = uh_0to2_42km
fout.variables['uh_0to2_prob_42km_hourly'][:] = uh_0to2_42km_hourly
fout.variables['wz_0to2_prob_9km'][:] = wz_0to2_9km
fout.variables['wz_0to2_prob_9km_hourly'][:] = wz_0to2_9km_hourly
fout.variables['wz_0to2_prob_27km'][:] = wz_0to2_27km
fout.variables['wz_0to2_prob_27km_hourly'][:] = wz_0to2_27km_hourly
fout.variables['wz_0to2_prob_42km'][:] = wz_0to2_42km
fout.variables['wz_0to2_prob_42km_hourly'][:] = wz_0to2_42km_hourly
fout.variables['comp_dz_prob_3km'][:] = comp_dz_3km
fout.variables['comp_dz_prob_3km_hourly'][:] = comp_dz_3km_hourly
fout.variables['comp_dz_prob_27km'][:] = comp_dz_27km
fout.variables['comp_dz_prob_27km_hourly'][:] = comp_dz_27km_hourly
fout.variables['comp_dz_prob_42km'][:] = comp_dz_42km
fout.variables['comp_dz_prob_42km_hourly'][:] = comp_dz_42km_hourly
fout.variables['rain_half_prob_3km'][:] = rain_half_3km
fout.variables['rain_half_prob_3km_hourly'][:] = rain_half_3km_hourly
fout.variables['rain_half_prob_27km'][:] = rain_half_27km
fout.variables['rain_half_prob_27km_hourly'][:] = rain_half_27km_hourly
fout.variables['rain_half_prob_42km'][:] = rain_half_42km
fout.variables['rain_half_prob_42km_hourly'][:] = rain_half_42km_hourly
fout.variables['rain_one_prob_3km'][:] = rain_one_3km
fout.variables['rain_one_prob_3km_hourly'][:] = rain_one_3km_hourly
fout.variables['rain_one_prob_27km'][:] = rain_one_27km
fout.variables['rain_one_prob_27km_hourly'][:] = rain_one_27km_hourly
fout.variables['rain_one_prob_42km'][:] = rain_one_42km
fout.variables['rain_one_prob_42km_hourly'][:] = rain_one_42km_hourly
fout.variables['rain_two_prob_3km'][:] = rain_two_3km
fout.variables['rain_two_prob_3km_hourly'][:] = rain_two_3km_hourly
fout.variables['rain_two_prob_27km'][:] = rain_two_27km
fout.variables['rain_two_prob_27km_hourly'][:] = rain_two_27km_hourly
fout.variables['rain_two_prob_42km'][:] = rain_two_42km
fout.variables['rain_two_prob_42km_hourly'][:] = rain_two_42km_hourly
fout.variables['rain_sat_half_prob_3km'][:] = rain_sat_half_3km
#fout.variables['rain_sat_half_prob_3km_hourly'][:] = rain_sat_half_3km_hourly
#fout.variables['rain_sat_half_prob_27km'][:] = rain_sat_half_27km
#fout.variables['rain_sat_half_prob_27km_hourly'][:] = rain_sat_half_27km_hourly
#fout.variables['rain_sat_half_prob_42km'][:] = rain_sat_half_42km
#fout.variables['rain_sat_half_prob_42km_hourly'][:] = rain_sat_half_42km_hourly
#fout.variables['rain_sat_one_prob_3km'][:] = rain_sat_one_3km
#fout.variables['rain_sat_one_prob_3km_hourly'][:] = rain_sat_one_3km_hourly
#fout.variables['rain_sat_one_prob_27km'][:] = rain_sat_one_27km
#fout.variables['rain_sat_one_prob_27km_hourly'][:] = rain_sat_one_27km_hourly
#fout.variables['rain_sat_one_prob_42km'][:] = rain_sat_one_42km
#fout.variables['rain_sat_one_prob_42km_hourly'][:] = rain_sat_one_42km_hourly
#fout.variables['rain_sat_two_prob_3km'][:] = rain_sat_two_3km
#fout.variables['rain_sat_two_prob_3km_hourly'][:] = rain_sat_two_3km_hourly
#fout.variables['rain_sat_two_prob_27km'][:] = rain_sat_two_27km
#fout.variables['rain_sat_two_prob_27km_hourly'][:] = rain_sat_two_27km_hourly
#fout.variables['rain_sat_two_prob_42km'][:] = rain_sat_two_42km
#fout.variables['rain_sat_two_prob_42km_hourly'][:] = rain_sat_two_42km_hourly
#fout.variables['soil_moisture_prob_3km'][:] = soil_moisture_3km
#fout.variables['soil_moisture_prob_3km_hourly'][:] = soil_moisture_3km_hourly
#fout.variables['soil_moisture_prob_27km'][:] = soil_moisture_27km
#fout.variables['soil_moisture_prob_27km_hourly'][:] = soil_moisture_27km_hourly
#fout.variables['soil_moisture_prob_42km'][:] = soil_moisture_42km
#fout.variables['soil_moisture_prob_42km_hourly'][:] = soil_moisture_42km_hourly
fout.variables['w_up_prob_9km'][:] = w_up_9km
fout.variables['w_up_prob_9km_hourly'][:] = w_up_9km_hourly
fout.variables['w_up_prob_27km'][:] = w_up_27km
fout.variables['w_up_prob_27km_hourly'][:] = w_up_27km_hourly
fout.variables['w_up_prob_42km'][:] = w_up_42km
fout.variables['w_up_prob_42km_hourly'][:] = w_up_42km_hourly
fout.variables['ws_80_prob_9km'][:] = ws_80_9km
fout.variables['ws_80_prob_9km_hourly'][:] = ws_80_9km_hourly
fout.variables['ws_80_prob_27km'][:] = ws_80_27km
fout.variables['ws_80_prob_27km_hourly'][:] = ws_80_27km_hourly
fout.variables['ws_80_prob_42km'][:] = ws_80_42km
fout.variables['ws_80_prob_42km_hourly'][:] = ws_80_42km_hourly
fout.variables['hail_prob_9km'][:] = hail_9km
fout.variables['hail_prob_9km_hourly'][:] = hail_9km_hourly
fout.variables['hail_prob_27km'][:] = hail_27km
fout.variables['hail_prob_27km_hourly'][:] = hail_27km_hourly
fout.variables['hail_prob_42km'][:] = hail_42km
fout.variables['hail_prob_42km_hourly'][:] = hail_42km_hourly
fout.variables['hailcast_prob_9km'][:] = hailcast_9km
fout.variables['hailcast_prob_9km_hourly'][:] = hailcast_9km_hourly
fout.variables['hailcast_prob_27km'][:] = hailcast_27km
fout.variables['hailcast_prob_27km_hourly'][:] = hailcast_27km_hourly
fout.variables['hailcast_prob_42km'][:] = hailcast_42km
fout.variables['hailcast_prob_42km_hourly'][:] = hailcast_42km_hourly

### Close output file ### 

fout.close()
del fout




