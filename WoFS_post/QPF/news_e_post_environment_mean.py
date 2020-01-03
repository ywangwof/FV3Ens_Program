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

neighborhood               = 15                 #15 gridpoint radius for probability matched mean 

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
      if (file[0:6] == 'wrfwof'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
         member_files.append(file)
   member_files.sort()

   files.append(os.path.join(temp_dir, member_files[t]))

files.sort()  #should have sorted directory paths to each ensemble file to be processed


##################### Process WRFOUT Files ########################

for f, infile in enumerate(files):

   try:                                                 #open WRFOUT file
      fin = netCDF4.Dataset(infile, "r")
      print "Opening %s \n" % infile
   except:
      print "%s does not exist! \n" %infile
      sys.exit(1)

   if (f == 0):

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
      timestep = str(t)
      if (len(timestep) == 1):
         timestep = '0' + timestep
      outname = "news-e_ENV_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + ".nc"         #output file
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

      if (valid_time_seconds < 25000.):         #Convert values past 0000 UTC, assumes forecast not run past ~7 UTC 
         valid_time_seconds = valid_time_seconds + 86400.

################################ Initialize output variables ########################################

      u_10 = np.zeros((ne,ny,nx))                                #10m U component of wind (m/s)
      v_10 = np.zeros((ne,ny,nx))                                #10m V component of wind (m/s)

      p_sfc = np.zeros((ne,ny,nx))                               #surface pressure (Pa)
      t_2 = np.zeros((ne,ny,nx))                                 #2m temp (K)
      td_2 = np.zeros((ne,ny,nx))                                 #2m temp (K)
      qv_2 = np.zeros((ne,ny,nx))                                #2m water vapor mixing ratio (g/kg)

      dbz_col_max = np.zeros((ne,ny,nx))                         #Simulated column max reflectivity (dBZ)     

      th_v_75mb = np.zeros((ne,ny,nx))                           #lowest 75mb mixed layer virtual potential temperature (K)
      th_e_75mb = np.zeros((ne,ny,nx))                           #lowest 75mb mixed layer equivalent potential temperature (K)

      lcl_75mb = np.zeros((ne,ny,nx))                            #lowest 75mb mixed layer lifted condensation level (m AGL)
      lfc_75mb = np.zeros((ne,ny,nx))                            #lowest 75mb mixed layer level of free convection (m AGL)
      cape_75mb = np.zeros((ne,ny,nx))                           #lowest 75mb mixed layer CAPE (J/kg)
      cin_75mb = np.zeros((ne,ny,nx))                            #lowest 75mb mixed layer CIN (J/kg)
      cape_0to3_75mb = np.zeros((ne,ny,nx))                      #lowest 75mb mixed layer 0-3 km AGL CAPE (J/kg)

      shear_u_0to1 = np.zeros((ne,ny,nx))                        #u-component of the 0-1 km AGL wind shear (m/s)
      shear_v_0to1 = np.zeros((ne,ny,nx))                        #v-component of the 0-1 km AGL wind shear (m/s)
      shear_u_0to6 = np.zeros((ne,ny,nx))                        #u-component of the 0-6 km AGL wind shear (m/s)
      shear_v_0to6 = np.zeros((ne,ny,nx))                        #v-component of the 0-6 km AGL wind shear (m/s)
      bunk_r_u = np.zeros((ne,ny,nx))                            #u-component of the Bunkers storm motion (right mover) (m/s)
      bunk_r_v = np.zeros((ne,ny,nx))                            #v-component of the Bunkers storm motion (right mover) (m/s)
      srh_0to1 = np.zeros((ne,ny,nx))                            #0-1 km AGL storm-relative helicity (m^2/s^2)
      srh_0to3 = np.zeros((ne,ny,nx))                            #0-3 km AGL storm-relative helicity (m^2/s^2)
      u_500 = np.zeros((ne,ny,nx))                               #u-component of the wind at 500 m AGL (m/s)      
      v_500 = np.zeros((ne,ny,nx))                               #v-component of the wind at 500 m AGL (m/s)      

      stp_75mb = np.zeros((ne,ny,nx))                            #Significant Tornado Parameter for 75mb mixed layer

      sw_down = np.zeros((ne,ny,nx))	        	      #Downward flux of shortwave radiation at the ground (W m^-2) 


##################### Process WRFOUT file: ########################

############################### Read WRFOUT variables: ##################################

   u_10[f,:,:] = fin.variables["U10"][0,:,:]              	#expects var dimensions of (nt, ny, nx) with nt = 1
   v_10[f,:,:] = fin.variables["V10"][0,:,:]

   p_sfc_temp = fin.variables["PSFC"][0,:,:]
   p_sfc_temp = p_sfc_temp / 100.
   p_sfc[f,:,:] = p_sfc_temp                                    #convert to hPa
  
   t_2[f,:,:] = fin.variables["T2"][0,:,:]
   th_2 = fin.variables["TH2"][0,:,:]
   qv_2_temp = fin.variables["Q2"][0,:,:]

   qv_2_temp = np.where(qv_2_temp < 0., 0.0001, qv_2_temp)                #force qv to be positive
   qv_2[f,:,:] = qv_2_temp
   td_2[f,:,:] = calc_td(t_2[f,:,:], p_sfc_temp, qv_2_temp)

   u = fin.variables["U"][0,:,:,:]		                #expects var dimensions of (nt, nz, ny, nx) with nt = 1
   v = fin.variables["V"][0,:,:,:]

   nz = u.shape[0]                         	        #get number of vertical grid levels

   ph = fin.variables["PH"][0,:,:,:]
   phb = fin.variables["PHB"][0,:,:,:]
   p = fin.variables["P"][0,:,:,:]
   pb = fin.variables["PB"][0,:,:,:]

   uc = (u[:,:,:-1]+u[:,:,1:])/2.                          #convert staggered grids to centered
   vc = (v[:,:-1,:]+v[:,1:,:])/2.

   qv = fin.variables["QVAPOR"][0,:,:,:]
   qc = fin.variables["QCLOUD"][0,:,:,:]
   qr = fin.variables["QRAIN"][0,:,:,:]
   qi = fin.variables["QICE"][0,:,:,:]
   qs = fin.variables["QSNOW"][0,:,:,:]
   qh = fin.variables["QGRAUP"][0,:,:,:]

   qt = qc + qr + qi + qs + qh
   qv = np.where(qv < 0., 0.0001, qv)                      #force qv to be positive definite

   th = fin.variables["T"][0,:,:,:]
   th = th + 300.                                  	#add base state temp (300 K) to potential temp

   dbz = fin.variables["REFL_10CM"][0,:,:,:]

   sw_down[f,:,:] = fin.variables["SWDOWN"][0,:,:]

### Close WRFOUT file ###

   fin.close()
   del fin

########################## Calculate derived output variables (using news_e_post_cbook.py): ##############################

######### Calculate vertical grid values #########

   z, dz = calc_height(ph, phb)                            #height and layer thickness (m)
   z_agl = z - hgt

   p = (p + pb) / 100.                                     #pressure (hPa)
   temp = calc_t(th, p)

######### Sfc/2m layer values #########

   t_v = calc_thv(temp, qv, qt)
   td = calc_td(temp, p, qv)
   th_v = calc_thv(th, qv, qt)
   temp_0to3 = np.ma.masked_where((z_agl > 3000.), (t_v))  #Set values above 3 km AGL to zero for calculating 0-3 km CAPE

######### 75mb mixed-layer values #########

   masked_temp = np.ma.masked_where((p_sfc_temp - p) > 75., (temp))
   masked_th = np.ma.masked_where((p_sfc_temp - p) > 75., (th))
   masked_td = np.ma.masked_where((p_sfc_temp - p) > 75., (td))
   masked_th_v = np.ma.masked_where((p_sfc_temp - p) > 75., (th_v))
   masked_t_v = np.ma.masked_where((p_sfc_temp - p) > 75., (t_v))
   masked_qv = np.ma.masked_where((p_sfc_temp - p) > 75., (qv))
   masked_p = np.ma.masked_where((p_sfc_temp - p) > 75., (p))

   t_75mb = np.ma.average(masked_temp, axis=0, weights=dz)
   th_75mb = np.ma.average(masked_th, axis=0, weights=dz)
   td_75mb = np.ma.average(masked_td, axis=0, weights=dz)
   th_v_75mb[f,:,:] = np.ma.average(masked_th_v, axis=0, weights=dz)
   t_v_75mb = np.ma.average(masked_t_v, axis=0, weights=dz)
   qv_75mb = np.ma.average(masked_qv, axis=0, weights=dz)
   p_75mb = np.ma.average(masked_p, axis=0, weights=dz)

   th_e_75mb[f,:,:], lcl_t_75mb = calc_the_bolt(p_75mb, t_75mb, qv_75mb)
   lcl_p_75mb =  1000. / (th_v_75mb[f,:,:] / lcl_t_75mb)**(1004. / 287.)

   lcl_up_index, lcl_low_index, lcl_interp = find_upper_lower(lcl_p_75mb, p)
   lcl_75mb[f,:,:] = calc_interp(z_agl, lcl_up_index, lcl_low_index, lcl_interp)

   t_75mb_parcel = calc_parcel_dj(p, th_e_75mb[f,:,:], t_v_75mb, p_75mb)

   lfc_p_75mb, el_p_75mb = calc_lfc_el(t_v, t_75mb_parcel, p, lcl_t_75mb, lcl_p_75mb)
   lfc_up_index, lfc_low_index, lfc_interp = find_upper_lower(lfc_p_75mb, p)
   lfc_75mb_temp = calc_interp(z_agl, lfc_up_index, lfc_low_index, lfc_interp)

   lfc_p_75mb = np.where(lfc_75mb_temp > 1000000., p_sfc_temp, lfc_p_75mb)
   lfc_75mb[f,:,:] = np.where(lfc_75mb_temp > 1000000., 0., lfc_75mb_temp)

   cape_75mb[f,:,:] = calc_cape(t_v, t_75mb_parcel, p, lcl_p_75mb, dz)
   cin_75mb[f,:,:] = calc_cin(t_v, t_75mb_parcel, p, lfc_p_75mb, dz)

   t_75mb_parcel_0to3 = np.ma.masked_where((z_agl > 3000.), (t_75mb_parcel))

   cape_0to3_75mb[f,:,:] = calc_cape(temp_0to3, t_75mb_parcel_0to3, p, lcl_p_75mb, dz)

######### Wind values #########

   temp_layer = np.zeros((z_agl.shape[1], z_agl.shape[2])) + 500.
   agl500_upper, agl500_lower, agl500_interp = find_upper_lower(temp_layer, z_agl)
   u_500[f,:,:] = calc_interp(uc, agl500_upper, agl500_lower, agl500_interp)
   v_500[f,:,:] = calc_interp(vc, agl500_upper, agl500_lower, agl500_interp)

   shear_u_0to1[f,:,:], shear_v_0to1[f,:,:] = calc_wind_shear(z_agl, uc, vc, 0., 1000.)
   shear_u_0to6[f,:,:], shear_v_0to6[f,:,:] = calc_wind_shear(z_agl, uc, vc, 0., 6000.)

   bunk_r_u[f,:,:], bunk_r_v[f,:,:], bunk_l_u, bunk_l_v = calc_bunkers(p, z_agl, dz, uc, vc)

   srh_0to1[f,:,:] = calc_srh(z_agl, uc, vc, dz, 0., 1000., bunk_r_u[f,:,:], bunk_r_v[f,:,:])
   srh_0to3[f,:,:] = calc_srh(z_agl, uc, vc, dz, 0., 3000., bunk_r_u[f,:,:], bunk_r_v[f,:,:])

   stp_75mb[f,:,:] = calc_stp(cape_75mb[f,:,:], lcl_75mb[f,:,:], srh_0to1[f,:,:], shear_u_0to6[f,:,:], shear_v_0to6[f,:,:])

######### Convert wind speeds to kts ##########

   u_10[f,:,:] = u_10[f,:,:] *  1.943844
   v_10[f,:,:] = v_10[f,:,:] *  1.943844
   u_500[f,:,:] = u_500[f,:,:] *  1.943844
   v_500[f,:,:] = v_500[f,:,:] *  1.943844
   shear_u_0to1[f,:,:] = shear_u_0to1[f,:,:] *  1.943844
   shear_v_0to1[f,:,:] = shear_v_0to1[f,:,:] *  1.943844
   shear_u_0to6[f,:,:] = shear_u_0to6[f,:,:] *  1.943844
   shear_v_0to6[f,:,:] = shear_v_0to6[f,:,:] *  1.943844
   bunk_r_u[f,:,:] = bunk_r_u[f,:,:] * 1.943844
   bunk_r_v[f,:,:] = bunk_r_v[f,:,:] * 1.943844

######### Storm-scale values #########

   dbz_col_max[f,:,:] = np.max(dbz, axis=0)


##################### Calc Ens. Mean Values: ########################

### Calc probability matched mean for reflectivity only ###

pmm_dz = dbz_col_max[0,:,:] * 0.
temp_dz = np.where(dbz_col_max > 100000., 0., dbz_col_max)
temp_mean_dz = np.mean(temp_dz, axis=0)

pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)


mean_cape = np.mean(cape_75mb, axis=0)
mean_cape_0to3 = np.mean(cape_0to3_75mb, axis=0)
mean_cin = np.mean(cin_75mb, axis=0)
mean_lcl = np.mean(lcl_75mb, axis=0)
mean_lfc = np.mean(lfc_75mb, axis=0)
mean_stp = np.mean(stp_75mb, axis=0)

mean_bunk_r_u = np.mean(bunk_r_u, axis=0)
mean_bunk_r_v = np.mean(bunk_r_v, axis=0)
mean_srh_0to1 = np.mean(srh_0to1, axis=0)
mean_srh_0to3 = np.mean(srh_0to3, axis=0)
mean_shear_u_0to6 = np.mean(shear_u_0to6, axis=0)
mean_shear_v_0to6 = np.mean(shear_v_0to6, axis=0)
mean_shear_u_0to1 = np.mean(shear_u_0to1, axis=0)
mean_shear_v_0to1 = np.mean(shear_v_0to1, axis=0)

mean_u_10 = np.mean(u_10, axis=0)
mean_v_10 = np.mean(v_10, axis=0)
mean_u_500 = np.mean(u_500, axis=0)
mean_v_500 = np.mean(v_500, axis=0)

mean_p_sfc = np.mean(p_sfc, axis=0)
mean_t_2 = np.mean(t_2, axis=0)
mean_td_2 = np.mean(td_2, axis=0)
mean_the_ml = np.mean(th_e_75mb, axis=0)
mean_thv_ml = np.mean(th_v_75mb, axis=0)
mean_qv_2 = np.mean(qv_2, axis=0)

mean_tf_2 = (mean_t_2 - 273.15) * 1.8 + 32.     #convert to deg. F
mean_tdf_2 = (mean_td_2 - 273.15) * 1.8 + 32.     #convert to deg. F

mean_sw_down = np.mean(sw_down, axis=0)

mean_comp_dz = np.mean(dbz_col_max, axis=0)

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

### Ensemble mean variables ###

u_10 = fout.createVariable('u_10', 'f4', ('NY','NX',))
u_10.long_name = "U-component of 10-m wind"
u_10.units = "kts"

v_10 = fout.createVariable('v_10', 'f4', ('NY','NX',))
v_10.long_name = "V-component of 10-m wind"
v_10.units = "kts"

p_sfc = fout.createVariable('p_sfc', 'f4', ('NY','NX',))
p_sfc.long_name = "Surface Pressure"
p_sfc.units = "hPa"

t_2 = fout.createVariable('t_2', 'f4', ('NY','NX',))
t_2.long_name = "2-m temperature"
t_2.units = "K"

td_2 = fout.createVariable('td_2', 'f4', ('NY','NX',))
td_2.long_name = "2-m dewpoint temperature"
td_2.units = "K"

qv_2 = fout.createVariable('qv_2', 'f4', ('NY','NX',))
qv_2.long_name = "2-m water vapor mixing ratio"
qv_2.units = "Kg/Kg"

th_v_2 = fout.createVariable('th_v_ml', 'f4', ('NY','NX',))
th_v_2.long_name = "75 hPa mixed layer virtual potential temperature"
th_v_2.units = "K"

th_e_2 = fout.createVariable('th_e_ml', 'f4', ('NY','NX',))
th_e_2.long_name = "75 hPa mixed layer equivalent potential temperature"
th_e_2.units = "K"

lcl_ml = fout.createVariable('lcl_ml', 'f4', ('NY','NX',))
lcl_ml.long_name = "75 hPa mixed layer lifted condensation level height"
lcl_ml.units = "m"

lcl_ml = fout.createVariable('lfc_ml', 'f4', ('NY','NX',))
lcl_ml.long_name = "75 hPa mixed layer level of free convection"
lcl_ml.units = "m"

cape_ml = fout.createVariable('cape_ml', 'f4', ('NY','NX',))
cape_ml.long_name = "75 hPa mixed layer CAPE"
cape_ml.units = "J/Kg"

cin_ml = fout.createVariable('cin_ml', 'f4', ('NY','NX',))
cin_ml.long_name = "75 hPa mixed layer convective inhibition"
cin_ml.units = "J/Kg"

cape_0to3_ml = fout.createVariable('cape_0to3_ml', 'f4', ('NY','NX',))
cape_0to3_ml.long_name = "75 hPa mixed layer CAPE below 3 km"
cape_0to3_ml.units = "J/Kg"

u_500 = fout.createVariable('u_500', 'f4', ('NY','NX',))
u_500.long_name = "U-component of 500-m wind"
u_500.units = "kts"

v_500 = fout.createVariable('v_500', 'f4', ('NY','NX',))
v_500.long_name = "V-component of 500-m wind"
v_500.units = "kts"

shear_u_0to1 = fout.createVariable('shear_u_0to1', 'f4', ('NY','NX',))
shear_u_0to1.long_name = "U-component of 0-1 km vertical wind shear"
shear_u_0to1.units = "kts"

shear_v_0to1 = fout.createVariable('shear_v_0to1', 'f4', ('NY','NX',))
shear_v_0to1.long_name = "V-component of 0-1 km vertical wind shear"
shear_v_0to1.units = "kts"

shear_u_0to6 = fout.createVariable('shear_u_0to6', 'f4', ('NY','NX',))
shear_u_0to6.long_name = "U-component of 0-6 km vertical wind shear"
shear_u_0to6.units = "kts"

shear_v_0to6 = fout.createVariable('shear_v_0to6', 'f4', ('NY','NX',))
shear_v_0to6.long_name = "V-component of 0-6 km vertical wind shear"
shear_v_0to6.units = "kts"

bunk_r_u = fout.createVariable('bunk_r_u', 'f4', ('NY','NX',))
bunk_r_u.long_name = "U-component Bunkers right-mover storm motion"
bunk_r_u.units = "kts"

bunk_r_v = fout.createVariable('bunk_r_v', 'f4', ('NY','NX',))
bunk_r_v.long_name = "V-component of Bunkers right-mover storm motion"
bunk_r_v.units = "kts"

srh_0to1 = fout.createVariable('srh_0to1', 'f4', ('NY','NX',))
srh_0to1.long_name = "0-1 km storm relative helicity"
srh_0to1.units = "m^2/s^2"

srh_0to3 = fout.createVariable('srh_0to3', 'f4', ('NY','NX',))
srh_0to3.long_name = "0-3 km storm relative helicity"
srh_0to3.units = "m^2/s^2"

stp_ml = fout.createVariable('stp_ml', 'f4', ('NY','NX',))
stp_ml.long_name = "75-hPa mixed layer significant tornado parameter"
stp_ml.units = "unitless"

sw_down = fout.createVariable('sw_down', 'f4', ('NY','NX',))
sw_down.long_name = "Downward flux of shortwave radiation at the ground"
sw_down.units = "W/m^2"

ccomp_dz = fout.createVariable('comp_dz', 'f4', ('NY','NX',))
ccomp_dz.long_name = "Mean composite reflectivity (93 km neighborhood)"
ccomp_dz.units = "dBZ"

ccomp_dz_pmm = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm.units = "dBZ"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['u_10'][:,:] = mean_u_10
fout.variables['v_10'][:,:] = mean_v_10

fout.variables['p_sfc'][:,:] = mean_p_sfc
fout.variables['t_2'][:,:] = mean_tf_2
fout.variables['td_2'][:,:] = mean_tdf_2
fout.variables['qv_2'][:,:] = mean_qv_2

fout.variables['th_v_ml'][:,:] = mean_thv_ml
fout.variables['th_e_ml'][:,:] = mean_the_ml

fout.variables['lcl_ml'][:,:] = mean_lcl
fout.variables['lfc_ml'][:,:] = mean_lfc
fout.variables['cape_ml'][:,:] = mean_cape
fout.variables['cin_ml'][:,:] = mean_cin
fout.variables['cape_0to3_ml'][:,:] = mean_cape_0to3

fout.variables['u_500'][:,:] = mean_u_500
fout.variables['v_500'][:,:] = mean_v_500
fout.variables['shear_u_0to1'][:,:] = mean_shear_u_0to1
fout.variables['shear_v_0to1'][:,:] = mean_shear_v_0to1
fout.variables['shear_u_0to6'][:,:] = mean_shear_u_0to6
fout.variables['shear_v_0to6'][:,:] = mean_shear_v_0to6
fout.variables['bunk_r_u'][:,:] = mean_bunk_r_u
fout.variables['bunk_r_v'][:,:] = mean_bunk_r_v
fout.variables['srh_0to1'][:,:] = mean_srh_0to1
fout.variables['srh_0to3'][:,:] = mean_srh_0to3
fout.variables['stp_ml'][:,:] = mean_stp

fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['comp_dz'][:] = mean_comp_dz

fout.variables['sw_down'][:,:] = mean_sw_down

### Close output file ### 

fout.close()
del fout




