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
import netCDF4
from optparse import OptionParser
from news_e_post_cbook import *

####################################### File Variables: ######################################################

parser = OptionParser()
parser.add_option("-d", dest="indir", type="string", default= None, help="Input directory of member directories containing WRFOUT files to process")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for summary file)")
parser.add_option("-t", dest="t", type="int", help = "Forecast timestep being processed")

(options, args) = parser.parse_args()

if ((options.indir == None) or (options.outdir == None) or (options.t == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   indir = options.indir
   outdir = options.outdir
   t = options.t

############################ Set Thresholds for convolution and neighborhood probabilities: #################################

radius_max_9km             = 3			#grid point radius for maximum value filter (3x3 square neighborhood)
radius_max_15km            = 5			#grid point radius for maximum value filter (5x5 square neighborhood)
radius_max_27km            = 9			#grid point radius for maximum value filter (9x9 square neighborhood)

radius_gauss               = 2			#grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

comp_dz_thresh             = 45.		#45 dBZ
rain_1_thresh              = 0.5		#0.5 inches
rain_1_3hr_thresh	   = 1.
rain_1_6hr_thresh	   = 1.5
rain_2_thresh              = 1.			#1 inch
rain_2_3hr_thresh	   = 2.
rain_2_6hr_thresh	   = 3.
rain_3_thresh              = 2.			#2 inches
rain_3_3hr_thresh	   = 3.
rain_3_6hr_thresh	   = 5.

#soil_moisture_thresh       = .95		#95% soil moisture

############################ Find WRFOUT files to process: #################################

### Find member dirs ### 

ne = 18
member_dirs = []

member_dirs_temp = os.listdir(indir)

for d, dir in enumerate(member_dirs_temp):
   if (dir[0:3] == 'ENS'):
      member_dirs.append(dir)

member_dirs.sort() #sorts as members [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]

files = []

for n in range(0, len(member_dirs)): 
   temp_dir = os.path.join(indir, member_dirs[n])

   member_files = []
   temp_files = os.listdir(temp_dir) 

   for f, file in enumerate(temp_files): 
      if (file[0:6] == 'wrfout'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
#      if (file[0:6] == 'wrfwof'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
         member_files.append(file)
   member_files.sort()

   files.append(os.path.join(temp_dir, member_files[t]))

files.sort()  #should have sorted directory paths to each ensemble file to be processed


##################### Process WRFOUT Files ########################

for f, infile in enumerate(files): 

   try:							#open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (f == 0):

##################### Get dimensions and attributes using first WRFOUT file ########################

      ### Get init and valid times ###

      start_date = fin.START_DATE				#Read initilization time string
      init_year = start_date[0:4]
      init_mon = start_date[5:7]
      init_day = start_date[8:10]
      init_hr = start_date[11:13]
      init_min = start_date[14:16]

      init_date = init_year + init_mon + init_day		#YYYYMMDD string for output file
      init_time = init_hr + init_min				#HHMM string for initialization time

      valid_hr = infile[-8:-6]				#Parse valid time from WRFOUT filename
      valid_min = infile[-5:-3]

      valid_time = valid_hr + valid_min			#HHMM string for valid time

      ### Set output path ###
      timestep = str(t) 
      if (len(timestep) == 1): 
         timestep = '0' + timestep
      outname = "news-e_ESR_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + ".nc"         #output file
      output_path = outdir + outname

      ### Get grid/projection info ### 

      dx = fin.DX                                             #east-west grid spacing (m)
      dy = fin.DY                                             #north-south grid spacing (m)
      cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
      cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
      stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
      true_lat_1 = fin.TRUELAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
      true_lat_2 = fin.TRUELAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)

      xlat = fin.variables["XLAT"][0,:,:]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["XLONG"][0,:,:]                    #longitude (dec deg; Lambert conformal)
      hgt = fin.variables["HGT"][0,:,:]                       #terrain height above MSL (m)

      ny = xlat.shape[0]
      nx = xlat.shape[1]

      ### Calculate initial and valid time in seconds ### 

      init_time_seconds = int(init_hr) * 3600. + int(init_min) * 60.
      valid_time_seconds = int(valid_hr) * 3600. + int(valid_min) * 60.

      if (valid_time_seconds < 32400.): 	#Convert values past 0000 UTC, assumes forecast not run past ~9 UTC 
         valid_time_seconds = valid_time_seconds + 86400. 

      ### Initialize ESR variables: ###

      comp_dz = np.zeros((ne,ny,nx))
      rain = np.zeros((ne,ny,nx))
#      soil_moisture = np.zeros((ne,ny,nx))
############################### Read WRFOUT variables: ##################################

   rain_temp = fin.variables["PREC_ACC_NC"][0,:,:]			#5-min accumulated precip
   rain[f,:,:] = rain_temp / 25.4					#convert to inches

#   soil_moisture[f,:,:] = fin.variables["SMOIS"][0,0,:,:]

   dbz = fin.variables["REFL_10CM"][0,:,:,:]

### Close WRFOUT file ###

   fin.close()
   del fin

########################## Calculate derived variables (using news_e_post_cbook.py): ##############################

   ######### Storm-scale values #########

   comp_dz[f,:,:] = np.max(dbz, axis=0)

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
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

### 90th percentile variables ###

rain_var = fout.createVariable('rain', 'f4', ('NE','NY','NX',))
rain_var.long_name = "5-minute Accumulated rainfall"
rain_var.units = "in"

#soil_moisture_var = fout.createVariable('soil_moisture', 'f4', ('NE','NY','NX',))
#soil_moisture_var.long_name = "Top layer soil_moisture"
#soil_moisture_var.units = "%"

comp_dz_var = fout.createVariable('comp_dz', 'f4', ('NE','NY','NX',))
comp_dz_var.long_name = "Composite reflectivity"
comp_dz_var.units = "dBZ"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['comp_dz'][:] = comp_dz
fout.variables['rain'][:] = rain
#fout.variables['soil_moisture'][:] = soil_moisture

### Close output file ### 

fout.close()
del fout




