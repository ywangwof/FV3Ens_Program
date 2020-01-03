###################################################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
import pylab as P
import numpy as np
import sys
import os
import netCDF4
import pickle as pl
from optparse import OptionParser
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *

parser = OptionParser()
parser.add_option("-d", dest="indir", type="string", default= None, help="Input directory of member directories containing WRFOUT files to process")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for summary file)")

(options, args) = parser.parse_args()

if ((options.indir == None) or (options.outdir == None)):
    print
    parser.print_help()
    print
    sys.exit(1)
else:
    indir = options.indir
    outdir = options.outdir

############### Get grid info from summary file: #############################

edge = 7
resolution = 'h'
area_thresh = 1000. 
damage_files = [] 

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
#      if (file[0:6] == 'wrfwof'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
      if (file[0:6] == 'wrfout'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
         member_files.append(file)
   member_files.sort()

   files.append(os.path.join(temp_dir, member_files[0]))

##################### Process WRFOUT Files ########################

infile = files[0]

try:                                                 #open WRFOUT file
   fin = netCDF4.Dataset(infile, "r")
   print "Opening %s \n" % infile
except:
   print "%s does not exist! \n" %infile
   sys.exit(1)


##################### Get dimensions and attributes using first WRFOUT file ########################

      ### Get init and valid times ###

start_date = fin.START_DATE                               #Read initilization time string
init_year = start_date[0:4]
init_mon = start_date[5:7]
init_day = start_date[8:10]
init_hr = start_date[11:13]
init_min = start_date[14:16]

init_date = init_year + init_mon + init_day               #YYYYMMDD string for output file
init_time = init_hr + init_min                            #HHMM string for initialization time

valid_hr = infile[-8:-6]                          #Parse valid time from WRFOUT filename
valid_min = infile[-5:-3]

valid_time = valid_hr + valid_min                 #HHMM string for valid time

      ### Set output path ###
outname = "news-e_MAP_" + init_date + "_" + init_time + ".pickle"         #output file
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

######################### Set domain: ####################################

xlat = xlat[edge:-edge,edge:-edge]
xlon = xlon[edge:-edge,edge:-edge]

sw_lat = xlat[0,0]
sw_lon = xlon[0,0]
ne_lat = xlat[-1,-1]
ne_lon = xlon[-1,-1]

ny = xlat.shape[0]
nx = xlat.shape[1]

fin.close()
del fin

print 'BASEMAP part'

m = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, urcrnrlon=ne_lon, urcrnrlat=ne_lat, projection='lcc', lat_1=true_lat_1, lat_2=true_lat_2, lat_0=cen_lat, lon_0=stand_lon, resolution = resolution, area_thresh = area_thresh)

print output_path
pl.dump(m,open(output_path,'wb'),-1)

