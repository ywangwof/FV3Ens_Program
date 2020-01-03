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

############################ Set Thresholds for convolution and neighborhood probabilities: #################################

radius_max_9km             = 3			#grid point radius for maximum value filter (3x3 square neighborhood)
radius_max_27km            = 9			#grid point radius for maximum value filter (9x9 square neighborhood)
radius_max_42km            = 14			#grid point radius for maximum value filter (14x14 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

comp_dz_thresh             = [35., 40., 45., 50., 55., 60.]		#40 dBZ
rain_thresh                = [0.01, 0.25, 0.5, 1., 2., 3.]		#0.5 inches
soil_moisture_thresh       = [.7, .75, .8, .85, .9, .95]		#95% soil moisture
hail_thresh                = [0.5, 1., 1.5, 2., 2.5, 3.] 		#1 inch
ws_80_thresh               = [38., 44., 50., 56., 62., 68.]		#58 kts
w_up_thresh                = [4., 10., 16., 22., 28., 34.]                #10 m/s
wz_0to2_thresh             = [0.002, 0.003, 0.004, 0.005, 0.006, 0.007]		#0.003 s^-1
uh_0to2_thresh             = [10., 25., 40., 65., 70., 85.]		#30 m^2/s^2
uh_2to5_thresh             = [20., 60., 100., 140., 180., 220.]		#60 m^2/s^2

perc                       = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]             #Ens. percentiles to calc/plot

############################ Find WRFOUT files to process: #################################

### Find ENS Summary files ### 

ne = 18

############hack to try and fix realtime for the last timestep: 
#if (t == fcst_nt):
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
#         print 'GOGOGOGOGOGOGOGO', tempfile
#         lastfile == 1
#         time.sleep(10) #give it 10 s to finish writing ... just in case
#      else:
#         print 'Not Yet, SUMMARY'
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
      outname = ens_file[0:7] + 'SUM' + ens_file[10:]
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

      ### Calc probability matched mean for reflectivity only ###
      temp_dz = np.where(comp_dz > 100000., 0., comp_dz)
      temp_mean_dz = np.mean(temp_dz, axis=0)
      pmm_dz = temp_mean_dz * 0.

      pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

   else: 

      uh_2to5                            = np.where((fin.variables["uh_2to5"][:,:,:] > uh_2to5), fin.variables["uh_2to5"][:,:,:], uh_2to5)
      uh_0to2                            = np.where((fin.variables["uh_0to2"][:,:,:] > uh_0to2), fin.variables["uh_0to2"][:,:,:], uh_0to2)
      wz_0to2                            = np.where((fin.variables["wz_0to2"][:,:,:] > wz_0to2), fin.variables["wz_0to2"][:,:,:], wz_0to2)
      comp_dz                            = np.where((fin.variables["comp_dz"][:,:,:] > comp_dz), fin.variables["comp_dz"][:,:,:], comp_dz)

      rain_temp                          = fin.variables["rain"][:,:,:]
      rain                               = rain + rain_temp

#      soil_moisture_temp                 = fin.variables["soil_moisture"][:,:,:]
#      soil_moisture                      = np.where((fin.variables["soil_moisture"][:,:,:] > soil_moisture), fin.variables["soil_moisture"][:,:,:], soil_moisture)

#      rain_sat_temp                      = np.where(soil_moisture_temp > 0.95, rain_temp, 0.)
#      rain_sat                           = rain_sat + rain_sat_temp

      w_up                               = np.where((fin.variables["w_up"][:,:,:] > w_up), fin.variables["w_up"][:,:,:], w_up)
      ws_80                              = np.where((fin.variables["ws_80"][:,:,:] > ws_80), fin.variables["ws_80"][:,:,:], ws_80)
      hail                               = np.where((fin.variables["hail"][:,:,:] > hail), fin.variables["hail"][:,:,:], hail)
      hailcast                           = np.where((fin.variables["hailcast"][:,:,:] > hailcast), fin.variables["hailcast"][:,:,:], hailcast)

   fin.close()
   del fin

##################### Calculate ensemble output variables: ########################

uh_2to5_perc = calc_perc_products(uh_2to5, perc, radius_max_9km, kernel)
uh_0to2_perc = calc_perc_products(uh_0to2, perc, radius_max_9km, kernel)
wz_0to2_perc = calc_perc_products(wz_0to2, perc, radius_max_9km, kernel)
comp_dz_perc = calc_perc_products(comp_dz, perc, 1, kernel)
rain_perc = calc_perc_products(rain, perc, 1, kernel)
#soil_moisture_perc = calc_perc_products(soil_moisture, perc, 1, kernel)
#rain_sat_perc = calc_perc_products(rain_sat, perc, 1, kernel)
ws_80_perc = calc_perc_products(ws_80, perc, radius_max_9km, kernel)
w_up_perc = calc_perc_products(w_up, perc, radius_max_9km, kernel)
hail_perc = calc_perc_products(hail, perc, radius_max_9km, kernel)
hailcast_perc = calc_perc_products(hailcast, perc, radius_max_9km, kernel)

uh_2to5_prob = calc_prob_products(uh_2to5, uh_2to5_thresh, radius_max_9km, kernel)
uh_0to2_prob = calc_prob_products(uh_0to2, uh_0to2_thresh, radius_max_9km, kernel)
wz_0to2_prob = calc_prob_products(wz_0to2, wz_0to2_thresh, radius_max_9km, kernel)
comp_dz_prob = calc_prob_products(comp_dz, comp_dz_thresh, 1, kernel)
rain_prob = calc_prob_products(rain, rain_thresh, 1, kernel)
#soil_moisture_prob = calc_prob_products(soil_moisture, soil_moisture_thresh, 1, kernel)
#rain_sat_prob = calc_prob_products(rain_sat, rain_thresh, 1, kernel)
ws_80_prob = calc_prob_products(ws_80, ws_80_thresh, radius_max_9km, kernel)
w_up_prob = calc_prob_products(w_up, w_up_thresh, radius_max_9km, kernel)
hail_prob = calc_prob_products(hail, hail_thresh, radius_max_9km, kernel)
hailcast_prob = calc_prob_products(hailcast, hail_thresh, radius_max_9km, kernel)

mean_uh_2to5 = np.mean(uh_2to5, axis=0)
mean_uh_0to2 = np.mean(uh_0to2, axis=0)
mean_wz_0to2 = np.mean(wz_0to2, axis=0)
mean_comp_dz = np.mean(comp_dz, axis=0)
mean_rain = np.mean(rain, axis=0)
#mean_soil_moisture = np.mean(soil_moisture, axis=0)
#mean_rain_sat = np.mean(rain_sat, axis=0)
mean_ws_80 = np.mean(ws_80, axis=0)
mean_w_up = np.mean(w_up, axis=0)
mean_hail = np.mean(hail, axis=0)
mean_hailcast = np.mean(hailcast, axis=0)

uh_2to5_pmm = prob_match_mean(uh_2to5, mean_uh_2to5, neighborhood)
uh_0to2_pmm = prob_match_mean(uh_0to2, mean_uh_0to2, neighborhood)
wz_0to2_pmm = prob_match_mean(wz_0to2, mean_wz_0to2, neighborhood)
comp_dz_pmm = prob_match_mean(comp_dz, mean_comp_dz, neighborhood)
rain_pmm = prob_match_mean(rain, mean_rain, neighborhood)
#soil_moisture_pmm = prob_match_mean(soil_moisture, mean_soil_moisture, neighborhood)
#rain_sat_pmm = prob_match_mean(rain_sat, mean_rain_sat, neighborhood)
ws_80_pmm = prob_match_mean(ws_80, mean_ws_80, neighborhood)
w_up_pmm = prob_match_mean(w_up, mean_w_up, neighborhood)
hail_pmm = prob_match_mean(hail, mean_hail, neighborhood)
hailcast_pmm = prob_match_mean(hailcast, mean_hailcast, neighborhood)

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NP', len(perc))
fout.createDimension('NR', len(uh_0to2_thresh))
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

### percentile variables ("p" is placeholder to differentiate variable name from one containing data ) ###

uh_2to5_pperc = fout.createVariable('uh_2to5_perc', 'f4', ('NP','NY','NX',))
uh_2to5_pperc.long_name = "Accumulated ensemble percentile values of 2-5 km updraft helicity"
uh_2to5_pperc.units = "m^2/s^2"

uh_0to2_pperc = fout.createVariable('uh_0to2_perc', 'f4', ('NP','NY','NX',))
uh_0to2_pperc.long_name = "Accumulated ensemble percentile values of 0-2 km updraft helicity"
uh_0to2_pperc.units = "m^2/s^2"

wz_0to2_pperc = fout.createVariable('wz_0to2_perc', 'f4', ('NP','NY','NX',))
wz_0to2_pperc.long_name = "Accumulated ensemble percentile values of average 0-2 km vertical vorticity"
wz_0to2_pperc.units = "s^-1"

rain_pperc = fout.createVariable('rain_perc', 'f4', ('NP','NY','NX',))
rain_pperc.long_name = "Accumulated ensemble percentile values of accumulated rainfall"
rain_pperc.units = "in"

#rain_sat_pperc = fout.createVariable('rain_sat_perc', 'f4', ('NP','NY','NX',))
#rain_sat_pperc.long_name = "Accumulated ensemble percentile values of accumulated rainfall on saturated soil"
#rain_sat_pperc.units = "in"

#soil_moisture_pperc = fout.createVariable('soil_moisture_perc', 'f4', ('NP','NY','NX',))
#soil_moisture_pperc.long_name = "Accumulated ensemble percentile values of soil_moisture"
#soil_moisture_pperc.units = "%"

comp_dz_pperc = fout.createVariable('comp_dz_perc', 'f4', ('NP','NY','NX',))
comp_dz_pperc.long_name = "Accumulated ensemble percentile values of composite reflectivity"
comp_dz_pperc.units = "dBZ"

w_up_pperc = fout.createVariable('w_up_perc', 'f4', ('NP','NY','NX',))
w_up_pperc.long_name = "Accumulated ensemble percentile values of updraft velocity"
w_up_pperc.units = "m/s"

ws_80_pperc = fout.createVariable('ws_80_perc', 'f4', ('NP','NY','NX',))
ws_80_pperc.long_name = "Accumulated ensemble percentile values of 80-m wind speed"
ws_80_pperc.units = "kts"

hail_pperc = fout.createVariable('hail_perc', 'f4', ('NP','NY','NX',))
hail_pperc.long_name = "Accumulated ensemble percentile values of max hail size at the surface (NSSL 2-moment)"
hail_pperc.units = "in"

hailcast_pperc = fout.createVariable('hailcast_perc', 'f4', ('NP','NY','NX',))
hailcast_pperc.long_name = "Accumulated ensemble percentile values of max hail size at the surface (Hailcast)"
hailcast_pperc.units = "in"

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

### probability variables ("p" is placeholder to differentiate variable name from one containing data ) ###

uh_2to5_pprob = fout.createVariable('uh_2to5_prob', 'f4', ('NR','NY','NX',))
uh_2to5_pprob.long_name = "Accumulated ensemble probabilities of 2-5 km updraft helicity"
uh_2to5_pprob.units = "%"

uh_0to2_pprob = fout.createVariable('uh_0to2_prob', 'f4', ('NR','NY','NX',))
uh_0to2_pprob.long_name = "Accumulated ensemble probabilities of 0-2 km updraft helicity"
uh_0to2_pprob.units = "%"

wz_0to2_pprob = fout.createVariable('wz_0to2_prob', 'f4', ('NR','NY','NX',))
wz_0to2_pprob.long_name = "Accumulated ensemble probabilities of average 0-2 km vertical vorticity"
wz_0to2_pprob.units = "%"

rain_pprob = fout.createVariable('rain_prob', 'f4', ('NR','NY','NX',))
rain_pprob.long_name = "Accumulated ensemble probabilities of accumulated rainfall"
rain_pprob.units = "%"

#rain_sat_pprob = fout.createVariable('rain_sat_prob', 'f4', ('NR','NY','NX',))
#rain_sat_pprob.long_name = "Accumulated ensemble probabilities of accumulated rainfall on saturated soil"
#rain_sat_pprob.units = "%"

#soil_moisture_pprob = fout.createVariable('soil_moisture_prob', 'f4', ('NR','NY','NX',))
#soil_moisture_pprob.long_name = "Accumulated ensemble probabilities of soil_moisture"
#soil_moisture_pprob.units = "%"

comp_dz_pprob = fout.createVariable('comp_dz_prob', 'f4', ('NR','NY','NX',))
comp_dz_pprob.long_name = "Accumulated ensemble probabilites of composite reflectivity"
comp_dz_pprob.units = "%"

w_up_pprob = fout.createVariable('w_up_prob', 'f4', ('NR','NY','NX',))
w_up_pprob.long_name = "Accumulated ensemble probabilities of updraft velocity"
w_up_pprob.units = "%"

ws_80_pprob = fout.createVariable('ws_80_prob', 'f4', ('NR','NY','NX',))
ws_80_pprob.long_name = "Accumulated ensemble probabilites of 80-m wind speed"
ws_80_pprob.units = "%"

hail_pprob = fout.createVariable('hail_prob', 'f4', ('NR','NY','NX',))
hail_pprob.long_name = "Accumulated ensemble probabilites of max hail size at the surface (NSSL 2-moment)"
hail_pprob.units = "%"

hailcast_pprob = fout.createVariable('hailcast_prob', 'f4', ('NR','NY','NX',))
hailcast_pprob.long_name = "Accumulated ensemble probabilities of max hail size at the surface (Hailcast)"
hailcast_pprob.units = "%"

### probability matched mean variables ("p" is placeholder to differentiate variable name from one containing data ) ###

uh_2to5_ppmm = fout.createVariable('uh_2to5_pmm', 'f4', ('NY','NX',))
uh_2to5_ppmm.long_name = "2-5 km updraft helicity probability matched mean"
uh_2to5_ppmm.units = "m^2/s^2"

uh_0to2_ppmm = fout.createVariable('uh_0to2_pmm', 'f4', ('NY','NX',))
uh_0to2_ppmm.long_name = "0-2 km updraft helicity probability matched mean"
uh_0to2_ppmm.units = "m^2/s^2"

wz_0to2_ppmm = fout.createVariable('wz_0to2_pmm', 'f4', ('NY','NX',))
wz_0to2_ppmm.long_name = "0-2 km vertical vorticity probability matched mean"
wz_0to2_ppmm.units = "s^-1"

rain_ppmm = fout.createVariable('rain_pmm', 'f4', ('NY','NX',))
rain_ppmm.long_name = "Accumulated rainfall probability matched mean"
rain_ppmm.units = "in"

#rain_sat_ppmm = fout.createVariable('rain_sat_pmm', 'f4', ('NY','NX',))
#rain_sat_ppmm.long_name = "Accumulated rainfall on saturated soil probability matched mean"
#rain_sat_ppmm.units = "in"

#soil_moisture_ppmm = fout.createVariable('soil_moisture_pmm', 'f4', ('NY','NX',))
#soil_moisture_ppmm.long_name = "Soil_moisture probability matched mean"
#soil_moisture_ppmm.units = "%"

comp_dz_ppmm_full = fout.createVariable('comp_dz_pmm_full', 'f4', ('NY','NX',))
comp_dz_ppmm_full.long_name = "Composite reflectivity probability matched mean"
comp_dz_ppmm_full.units = "dBZ"

w_up_ppmm = fout.createVariable('w_up_pmm', 'f4', ('NY','NX',))
w_up_ppmm.long_name = "Updraft velocity probability matched mean"
w_up_ppmm.units = "m/s"

ws_80_ppmm = fout.createVariable('ws_80_pmm', 'f4', ('NY','NX',))
ws_80_ppmm.long_name = "80-m wind speed probability matched mean"
ws_80_ppmm.units = "m/s"

hail_ppmm = fout.createVariable('hail_pmm', 'f4', ('NY','NX',))
hail_ppmm.long_name = "Max hail size at the surface (NSSL 2-moment) probability matched mean"
hail_ppmm.units = "in"

hailcast_ppmm = fout.createVariable('hailcast_pmm', 'f4', ('NY','NX',))
hailcast_ppmm.long_name = "Max hail size at the surface (Hailcast) probability matched mean"
hailcast_ppmm.units = "in"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['uh_2to5_perc'][:] = uh_2to5_perc
fout.variables['uh_0to2_perc'][:] = uh_0to2_perc
fout.variables['wz_0to2_perc'][:] = wz_0to2_perc
fout.variables['comp_dz_perc'][:] = comp_dz_perc
fout.variables['rain_perc'][:] = rain_perc
#fout.variables['rain_sat_perc'][:] = rain_sat_perc
#fout.variables['soil_moisture_perc'][:] = soil_moisture_perc
fout.variables['w_up_perc'][:] = w_up_perc
fout.variables['ws_80_perc'][:] = ws_80_perc
fout.variables['hail_perc'][:] = hail_perc
fout.variables['hailcast_perc'][:] = hailcast_perc

fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['uh_2to5_prob'][:] = uh_2to5_prob
fout.variables['uh_0to2_prob'][:] = uh_0to2_prob
fout.variables['wz_0to2_prob'][:] = wz_0to2_prob
fout.variables['comp_dz_prob'][:] = comp_dz_prob
fout.variables['rain_prob'][:] = rain_prob
#fout.variables['rain_sat_prob'][:] = rain_sat_prob
#fout.variables['soil_moisture_prob'][:] = soil_moisture_prob
fout.variables['w_up_prob'][:] = w_up_prob
fout.variables['ws_80_prob'][:] = ws_80_prob
fout.variables['hail_prob'][:] = hail_prob
fout.variables['hailcast_prob'][:] = hailcast_prob

fout.variables['uh_2to5_pmm'][:] = uh_2to5_pmm
fout.variables['uh_0to2_pmm'][:] = uh_0to2_pmm
fout.variables['wz_0to2_pmm'][:] = wz_0to2_pmm
fout.variables['comp_dz_pmm_full'][:] = comp_dz_pmm
fout.variables['rain_pmm'][:] = rain_pmm
#fout.variables['rain_sat_pmm'][:] = rain_sat_pmm
#fout.variables['soil_moisture_pmm'][:] = soil_moisture_pmm
fout.variables['w_up_pmm'][:] = w_up_pmm
fout.variables['ws_80_pmm'][:] = ws_80_pmm
fout.variables['hail_pmm'][:] = hail_pmm
fout.variables['hailcast_pmm'][:] = hailcast_pmm

### Close output file ### 

fout.close()
del fout




