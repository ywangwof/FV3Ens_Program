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
ws_80_thresh1              = 50.		#58 mph
w_up_thresh                = 10.                #10 m/s
wz_0to2_thresh             = 0.003		#0.003 s^-1
uh_0to2_thresh             = 30.		#30 m^2/s^2
uh_2to5_thresh             = 60.		#60 m^2/s^2


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
      if (file[0:9] == 'wrfout_d0'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
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
      outname = "news-e_ENS_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + ".nc"         #output file
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

      if (valid_time_seconds < 25000.): 	#Convert values past 0000 UTC, assumes forecast not run past ~7 UTC 
         valid_time_seconds = valid_time_seconds + 86400. 

      ### Initialize ENS variables: ###

      uh_2to5 = np.zeros((ne,ny,nx))
      uh_0to2 = np.zeros((ne,ny,nx))
      wz_0to2 = np.zeros((ne,ny,nx))
      w_up = np.zeros((ne,ny,nx))
      ws_80 = np.zeros((ne,ny,nx))
      hail = np.zeros((ne,ny,nx))
      hailcast = np.zeros((ne,ny,nx))
      comp_dz = np.zeros((ne,ny,nx))
      rain = np.zeros((ne,ny,nx))
      soil_moisture = np.zeros((ne,ny,nx))


############################### Read WRFOUT variables: ##################################

   ws_80_temp = fin.variables["WSPD80"][0,:,:]
   ws_80[f,:,:] = ws_80_temp * 1.943844                               #convert to kts

   w_up[f,:,:] = fin.variables["W_UP_MAX"][0,:,:]

   rain_temp = fin.variables["PREC_ACC_NC"][0,:,:]			#5-min accumulated precip
   rain[f,:,:] = rain_temp / 25.4					#convert to inches

   soil_moisture[f,:,:] = fin.variables["SMOIS"][0,0,:,:]

   hail_temp = fin.variables["HAIL_MAXK1"][0,:,:]
   hail[f,:,:] = hail_temp * 100. / 2.54                               #convert to inches

   hailcast_temp = fin.variables["HAILCAST_DIAM_MAX"][0,:,:]
   hailcast[f,:,:] = hailcast_temp / 10. / 2.54                        #convert to inches

  #u, v not needed if "REL_VORT" output available
#   u = fin.variables["U"][0,:,:,:]		                #expects var dimensions of (nt, nz, ny, nx) with nt = 1
#   v = fin.variables["V"][0,:,:,:]
   w = fin.variables["W"][0,:,:,:]


   ph = fin.variables["PH"][0,:,:,:]
   phb = fin.variables["PHB"][0,:,:,:]
#   p = fin.variables["P"][0,:,:,:]
#   pb = fin.variables["PB"][0,:,:,:]

#   uc = (u[:,:,:-1]+u[:,:,1:])/2.                          #convert staggered grids to centered
#   vc = (v[:,:-1,:]+v[:,1:,:])/2.
   wc = (w[:-1,:,:]+w[1:,:,:])/2.

   nz = wc.shape[0]                         	        #get number of vertical grid levels

   dbz = fin.variables["REFL_10CM"][0,:,:,:]
   wz = fin.variables["REL_VORT"][0,:,:,:]

###ffair
#   wz = vc * 0.
#   wz[:,:-1,:-1] = (vc[:,:-1,1:]-vc[:,:-1,:-1])/dx - (uc[:,1:,:-1]-uc[:,:-1,:-1])/dy

### Close WRFOUT file ###

   fin.close()
   del fin

########################## Calculate derived variables (using news_e_post_cbook.py): ##############################

   ######### Calculate vertical grid values #########

   z, dz = calc_height(ph, phb)                            #height and layer thickness (m)
   z_agl = z - hgt

   ######### Storm-scale values #########

   comp_dz[f,:,:] = np.max(dbz, axis=0)

   wz_0to2[f,:,:] = calc_avg_vort(wz, z_agl, dz, 0., 2000.)
   uh_0to2[f,:,:] = calc_uh(wc, wz, z_agl, dz, 0., 2000.)
   uh_2to5[f,:,:] = calc_uh(wc, wz, z_agl, dz, 2000., 5000.)


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

uh_2to5_var = fout.createVariable('uh_2to5', 'f4', ('NE','NY','NX',))
uh_2to5_var.long_name = "2-5 km updraft helicity"
uh_2to5_var.units = "m^2/s^2"

uh_0to2_var = fout.createVariable('uh_0to2', 'f4', ('NE','NY','NX',))
uh_0to2_var.long_name = "0-2 km updraft helicity"
uh_0to2_var.units = "m^2/s^2"

wz_0to2_var = fout.createVariable('wz_0to2', 'f4', ('NE','NY','NX',))
wz_0to2_var.long_name = "Average 0-2 km vertical vorticity"
wz_0to2_var.units = "s^-1"

rain_var = fout.createVariable('rain', 'f4', ('NE','NY','NX',))
rain_var.long_name = "5-minute Accumulated rainfall"
rain_var.units = "in"

soil_moisture_var = fout.createVariable('soil_moisture', 'f4', ('NE','NY','NX',))
soil_moisture_var.long_name = "Top layer soil_moisture"
soil_moisture_var.units = "%"

comp_dz_var = fout.createVariable('comp_dz', 'f4', ('NE','NY','NX',))
comp_dz_var.long_name = "Composite reflectivity"
comp_dz_var.units = "dBZ"

w_up_var = fout.createVariable('w_up', 'f4', ('NE','NY','NX',))
w_up_var.long_name = "Max updraft velocity"
w_up_var.units = "m/s"

ws_80_var = fout.createVariable('ws_80', 'f4', ('NE','NY','NX',))
ws_80_var.long_name = "80-m wind speed"
ws_80_var.units = "kts"

hail_var = fout.createVariable('hail', 'f4', ('NE','NY','NX',))
hail_var.long_name = "Max hail size at the surface (NSSL 2-moment)"
hail_var.units = "in"

hailcast_var = fout.createVariable('hailcast', 'f4', ('NE','NY','NX',))
hailcast_var.long_name = "Max hail size at the surface (Hailcast)"
hailcast_var.units = "in"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['uh_2to5'][:] = uh_2to5
fout.variables['uh_0to2'][:] = uh_0to2
fout.variables['wz_0to2'][:] = wz_0to2
fout.variables['comp_dz'][:] = comp_dz
fout.variables['rain'][:] = rain
fout.variables['soil_moisture'][:] = soil_moisture
fout.variables['w_up'][:] = w_up
fout.variables['ws_80'][:] = ws_80
fout.variables['hail'][:] = hail
fout.variables['hailcast'][:] = hailcast

### Close output file ### 

fout.close()
del fout




