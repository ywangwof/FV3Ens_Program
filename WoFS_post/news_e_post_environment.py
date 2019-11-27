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
radius_max                 = 5
radius_gauss               = 2                  #grid point radius of convolution operator

kernel                     = gauss_kern(radius_gauss)   #convolution operator kernel

############################ Find WRFOUT files to process: #################################

### Find member dirs ###

ne = 40
member_dirs = []

member_dirs_temp = os.listdir(indir)

for d, dir in enumerate(member_dirs_temp):
   if (dir[0:4] == 'mem_'):
      member_dirs.append(dir)

member_dirs.sort() #sorts as members [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]

files = []

for n in range(0, len(member_dirs)):
   temp_dir = os.path.join(indir, member_dirs[n],'summary')

   member_files = []
   temp_files = os.listdir(temp_dir)

   for f, file in enumerate(temp_files):
      if (file[0:7] == 'FV3_ENS'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
#      if (file[0:6] == 'wrfwof'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
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
      init_mon = start_date[4:6]
      init_day = start_date[6:8]
      init_hr = start_date[8:10]
      init_min = start_date[10:12]

      init_date = init_year + init_mon + init_day               #YYYYMMDD string for output file
      init_time = init_hr + init_min                            #HHMM string for initialization time

      valid_hr = infile[-7:-5]				#Parse valid time from WRFOUT filename
      valid_min = infile[-5:-3]

      valid_time = valid_hr + valid_min                 #HHMM string for valid time

      ### Set output path ###
      timestep = str(t)
      if (len(timestep) == 1):
         timestep = '0' + timestep
      outname = "fv3sar_ENV_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + ".nc"         #output file
      output_path = os.path.join(outdir,'enspost')
      if not os.path.lexists(output_path): os.mkdir(output_path)
      output_path =  os.path.join(outdir,'enspost',outname)

      ### Get grid/projection info ###

      dx = fin.DX                                             #east-west grid spacing (m)
      dy = fin.DY                                             #north-south grid spacing (m)
      cen_lat = fin.CEN_LAT                                   #center of domain latitude (dec deg)
      cen_lon = fin.CEN_LON                                   #center of domain longitude (dec deg)
      stand_lon = fin.STAND_LON                               #center lon of Lambert conformal projection
      true_lat_1 = fin.TRUELAT1                               #true lat value 1 for Lambert conformal conversion (dec deg)
      true_lat_2 = fin.TRUELAT2                               #true lat value 2 for Lambert conformal conversion (dec deg)

      xlat = fin.variables["xlat"][:,:]                     #latitude (dec deg; Lambert conformal)
      xlon = fin.variables["xlon"][:,:]                    #longitude (dec deg; Lambert conformal)
      hgt = fin.variables["hgt"][:,:]                       #terrain height above MSL (m)

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
      u_bl = np.zeros((ne,ny,nx))                                #U component of wind at boundary layer (1000-925 hPa)
      v_bl = np.zeros((ne,ny,nx))                                #V component of wind at boundary layer (1000-925 hPa)
      u_850p = np.zeros((ne,ny,nx))                              #850-hPa U component of wind (m/s)
      u_700p = np.zeros((ne,ny,nx))                              #700-hPa U component of wind (m/s)
      v_850p = np.zeros((ne,ny,nx))                              #850-hPa V component of wind (m/s)
      v_700p = np.zeros((ne,ny,nx))                              #700-hPa V component of wind (m/s)

      mslp = np.zeros((ne,ny,nx))                                #2m water vapor mixing ratio (g/kg)
      p_sfc = np.zeros((ne,ny,nx))                               #surface pressure (Pa)
      t_2 = np.zeros((ne,ny,nx))                                 #2m temp (K)
      td_2 = np.zeros((ne,ny,nx))                                 #2m temp (K)
      t_850p = np.zeros((ne,ny,nx))                              #850-hPa temp (K)
      t_700p = np.zeros((ne,ny,nx))                              #700-hPa temp (K)
      td_850p = np.zeros((ne,ny,nx))                             #850-hPa dewpt temp (K)
      td_700p = np.zeros((ne,ny,nx))                             #700-hPa dewpt temp (K)

      qv_2 = np.zeros((ne,ny,nx))                                #2m water vapor mixing ratio (g/kg)

#      pbl_hgt = np.zeros((ne,ny,nx))                             #PBL height (m)  #not in WRFWOF files yet

      dbz_col_max = np.zeros((ne,ny,nx))                         #Simulated column max reflectivity (dBZ)

      th_v_100mb = np.zeros((ne,ny,nx))                           #lowest 100mb mixed layer virtual potential temperature (K)
      th_e_100mb = np.zeros((ne,ny,nx))                           #lowest 100mb mixed layer equivalent potential temperature (K)

      lcl_100mb = np.zeros((ne,ny,nx))                            #lowest 100mb mixed layer lifted condensation level (m AGL)
      lfc_100mb = np.zeros((ne,ny,nx))                            #lowest 100mb mixed layer level of free convection (m AGL)
      cape_100mb = np.zeros((ne,ny,nx))                           #lowest 100mb mixed layer CAPE (J/kg)
      cin_100mb = np.zeros((ne,ny,nx))                            #lowest 100mb mixed layer CIN (J/kg)
      cape_0to3_100mb = np.zeros((ne,ny,nx))                      #lowest 100mb mixed layer 0-3 km AGL CAPE (J/kg)
      pbl_mfc = np.zeros((ne,ny,nx))                             # Moisture flux convergence (g/(kg*s))

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

      stp_100mb = np.zeros((ne,ny,nx))                            #Significant Tornado Parameter for 100mb mixed layer

      sw_down = np.zeros((ne,ny,nx))	          	         #Downward flux of shortwave radiation at the ground (W m^-2)
      iwp = np.zeros((ne,ny,nx))	          	         #Cloud ice water path (ice water path units)
      lwp = np.zeros((ne,ny,nx))	          	         #Cloud liquid water path (liquid water path units)
      pw = np.zeros((ne,ny,nx))	        		         #Precipitable water (in)
      ctp = np.zeros((ne,ny,nx))	        	         #Cloud top pressure (hPa)

      u_dvg = np.zeros((ne,ny,nx))                               #u-component of upper-level (400-250 hPa) divergence (m/s) --- BCM
      v_dvg = np.zeros((ne,ny,nx))                               #v-component of upper-level (400-250 hPa) divergence (m/s) --- BCM
      corf_us_u = np.zeros((ne,ny,nx))                           #u-component of the Corfidi upshear vector (m/s) --- BCM
      corf_us_v = np.zeros((ne,ny,nx))                           #v-component of the Corfidi upshear vector (m/s) --- BCM
      corf_ds_u = np.zeros((ne,ny,nx))                           #u-component of the Corfidi downshear vector (m/s) --- BCM
      corf_ds_v = np.zeros((ne,ny,nx))                           #v-component of the Corfidi downshear vector (m/s) --- BCM

      u_850_300p = np.zeros((ne,ny,nx))                          #u-component of the 850-300 hPa mean wind (m/s) --- BCM
      v_850_300p = np.zeros((ne,ny,nx))                          #v-component of the 850-300 hPa mean wind (m/s) --- BCM
      w_700p = np.zeros((ne,ny,nx))                              #700-hPa vertical velocity (omega) --- BCM
      u_500p = np.zeros((ne,ny,nx))                              #u-component of the 500-hPa wind (m/s) --- BCM
      v_500p = np.zeros((ne,ny,nx))                              #v-component of the 500-hPa wind (m/s) --- BCM

##################### Process WRFOUT file: ########################

############################### Read WRFOUT variables: ##################################

   u_10[f,:,:] = fin.variables["u10"][:,:]              	#expects var dimensions of (nt, ny, nx) with nt = 1
   v_10[f,:,:] = fin.variables["v10"][:,:]

   temps =  fin.variables['t2'][:,:]
   stemps = temps[:,:]+6.5*fin.variables['hgt'][:,:]/1000.
   mslp[f,:,:] = fin.variables['psfc'][:,:]*np.exp(9.81/(287.0*stemps)*fin.variables['hgt'][:,:])*0.01 + (6.7 * fin.variables['hgt'][:,:] / 1000.)

   p_sfc_temp = fin.variables["psfc"][:,:]
   p_sfc_temp = p_sfc_temp / 100.
   p_sfc[f,:,:] = p_sfc_temp                                    #convert to hpa

   t_2[f,:,:] = fin.variables["t2"][:,:]
   th_2 = fin.variables["th2"][:,:]
   qv_2_temp = fin.variables["q2"][:,:]

   qv_2_temp = np.where(qv_2_temp < 0., 0.0001, qv_2_temp)                #force qv to be positive
   qv_2[f,:,:] = qv_2_temp
   td_2[f,:,:] = calc_td(t_2[f,:,:], p_sfc_temp, qv_2_temp)

   uc = fin.variables["u"][:,:,:]		                #expects var dimensions of (nt, nz, ny, nx) with nt = 1
   vc = fin.variables["v"][:,:,:]
   wc = fin.variables["w"][:,:,:]

   nz = uc.shape[0]                         	        #get number of vertical grid levels

   dz = fin.variables["delz"][:,:,:]
   z = np.cumsum(dz, axis=0)
   #phb = np.zeros((nz,ny,nz))
   p = fin.variables["p"][:,:,:]
   pb = fin.variables["pb"][:,:,:]

   #uc = (u[:,:,:-1]+u[:,:,1:])/2.                          #convert staggered grids to centered
   #vc = (v[:,:-1,:]+v[:,1:,:])/2.
   #wc = (w[:-1,:,:]+w[1:,:,:])/2.

   qv = fin.variables["qvapor"][:,:,:]
   qc = fin.variables["qcloud"][:,:,:]
   qr = fin.variables["qrain"][:,:,:]
   qi = fin.variables["qice"][:,:,:]
   qs = fin.variables["qsnow"][:,:,:]
   qh = fin.variables["qgraupel"][:,:,:]

   qt = qc + qr + qi + qs + qh
   qv = np.where(qv < 0., 0.0001, qv)                      #force qv to be positive definite

   th = fin.variables["t"][:,:,:]
   th = th + 300.                                  	#add base state temp (300 k) to potential temp
   print 'Min max th = %f  %f' % (np.min(th), np.max(th))
   dbz = fin.variables["refl_10cm"][:,:,:]

#   pbl_hgt[f,:,:] = fin.variables["pblh"][:,:]       #not in wrfwof files yet
   sw_down[f,:,:] = fin.variables["swdown"][:,:]
   iwp[f,:,:] = fin.variables["iwp"][:,:]
   lwp[f,:,:] = fin.variables["lwp"][:,:]

### Close WRFOUT file ###

   fin.close()
   del fin

########################## Calculate derived output variables (using news_e_post_cbook.py): ##############################

######### Calculate vertical grid values #########

   #z, dz = calc_height(ph, phb)                            #height and layer thickness (m)
   z_agl = z - hgt
   p = (p + pb) / 100.                                     #pressure (hPa)
   temp = calc_t(th, p)

######### Sfc/2m layer values #########

   t_v = calc_thv(temp, qv, qt)
   td = calc_td(temp, p, qv)
   th_v = calc_thv(th, qv, qt)
   print 'Min max th_v = %f  %f' % (np.min(th_v), np.max(th_v))
   temp_0to3 = np.ma.masked_where((z_agl > 3000.), (t_v))  #Set values above 3 km AGL to zero for calculating 0-3 km CAPE

######### Column integrated values #########

   pw[f,:,:] = calc_pw(qv, p)
   ctp[f,:,:] = calc_ctp(qt, p)
#   print 'sadlfk;jas', pw.shape, qv.shape, dz.shape

######### 100mb mixed-layer values #########

   masked_temp = np.ma.masked_where((p_sfc_temp - p) > 100., (temp))
   masked_th = np.ma.masked_where((p_sfc_temp - p) > 100., (th))
   masked_td = np.ma.masked_where((p_sfc_temp - p) > 100., (td))
   masked_th_v = np.ma.masked_where((p_sfc_temp - p) > 100., (th_v))
   masked_t_v = np.ma.masked_where((p_sfc_temp - p) > 100., (t_v))
   masked_qv = np.ma.masked_where((p_sfc_temp - p) > 100., (qv))
   masked_p = np.ma.masked_where((p_sfc_temp - p) > 100., (p))

   t_100mb = np.ma.average(masked_temp, axis=0, weights=dz)
   th_100mb = np.ma.average(masked_th, axis=0, weights=dz)
   td_100mb = np.ma.average(masked_td, axis=0, weights=dz)
   th_v_100mb[f,:,:] = np.ma.average(masked_th_v, axis=0, weights=dz)
   t_v_100mb = np.ma.average(masked_t_v, axis=0, weights=dz)
   qv_100mb = np.ma.average(masked_qv, axis=0, weights=dz)
   p_100mb = np.ma.average(masked_p, axis=0, weights=dz)

   th_e_100mb[f,:,:], lcl_t_100mb = calc_the_bolt(p_100mb, t_100mb, qv_100mb)
   lcl_p_100mb =  1000. / (th_v_100mb[f,:,:] / lcl_t_100mb)**(1004. / 287.)

   lcl_up_index, lcl_low_index, lcl_interp = find_upper_lower(lcl_p_100mb, p)
   lcl_100mb[f,:,:] = calc_interp(z_agl, lcl_up_index, lcl_low_index, lcl_interp)

   t_100mb_parcel = calc_parcel_dj(p, th_e_100mb[f,:,:], t_v_100mb, p_100mb)

   lfc_p_100mb, el_p_100mb = calc_lfc_el(t_v, t_100mb_parcel, p, lcl_t_100mb, lcl_p_100mb)
   lfc_up_index, lfc_low_index, lfc_interp = find_upper_lower(lfc_p_100mb, p)
   lfc_100mb_temp = calc_interp(z_agl, lfc_up_index, lfc_low_index, lfc_interp)

   lfc_p_100mb = np.where(lfc_100mb_temp > 1000000., p_sfc_temp, lfc_p_100mb)
   lfc_100mb[f,:,:] = np.where(lfc_100mb_temp > 1000000., 0., lfc_100mb_temp)

   cape_100mb[f,:,:] = calc_cape(t_v, t_100mb_parcel, p, lcl_p_100mb, dz)
   cin_100mb[f,:,:] = calc_cin(t_v, t_100mb_parcel, p, lfc_p_100mb, dz)

   t_100mb_parcel_0to3 = np.ma.masked_where((z_agl > 3000.), (t_100mb_parcel))

   cape_0to3_100mb[f,:,:] = calc_cape(temp_0to3, t_100mb_parcel_0to3, p, lcl_p_100mb, dz)

###### Boundary-layer values ######

   pbl_mfc[f,:,:] = calc_bl_moisture_conv(p, uc, vc, qv, p_sfc_temp)
   neigh_pbl_mfc = pbl_mfc * 0.

   for n in range(0, pbl_mfc.shape[0]):
      pbl_mfc_temp = get_local_maxima2d(pbl_mfc[n,:,:], radius_max)
      neigh_pbl_mfc[n,:,:] = signal.convolve2d(pbl_mfc_temp, kernel, 'same')

   masked_u_bl = np.ma.masked_where((p_sfc_temp - p) > 75, (uc)) # Used for MFC
   masked_v_bl = np.ma.masked_where((p_sfc_temp - p) > 75, (vc)) # Used for MFC

   u_bl[f,:,:] = np.ma.masked_invalid(np.ma.average(masked_u_bl, axis=0, weights=dz))
   v_bl[f,:,:] = np.ma.masked_invalid(np.ma.average(masked_v_bl, axis=0, weights=dz))

###### Pressure-level temp and dewpoint
   t_850p[f,:,:] = calc_var_at_plev(temp, pb, p, 850.)
   t_700p[f,:,:] = calc_var_at_plev(temp, pb, p, 700.)
   td_850p[f,:,:] = calc_var_at_plev(td, pb, p, 850.)
   td_700p[f,:,:] = calc_var_at_plev(td, pb, p, 700.)

   t_850p = np.ma.masked_where(t_850p < 50, t_850p)
   t_700p = np.ma.masked_where(t_700p < 50, t_700p)
   td_850p = np.ma.masked_where(td_850p < 50, td_850p)
   td_700p = np.ma.masked_where(td_700p < 50, td_700p)


######### Wind values #########

   temp_layer = np.zeros((z_agl.shape[1], z_agl.shape[2])) + 500.
   agl500_upper, agl500_lower, agl500_interp = find_upper_lower(temp_layer, z_agl)
   u_500[f,:,:] = calc_interp(uc, agl500_upper, agl500_lower, agl500_interp)
   v_500[f,:,:] = calc_interp(vc, agl500_upper, agl500_lower, agl500_interp)

   u_850p[f,:,:] = calc_var_at_plev(uc, p, pb, 850.)
   u_700p[f,:,:] = calc_var_at_plev(uc, p, pb, 700.)
   u_500p[f,:,:] = calc_var_at_plev(uc, p, pb, 500.)
   v_850p[f,:,:] = calc_var_at_plev(vc, p, pb, 850.)
   v_700p[f,:,:] = calc_var_at_plev(vc, p, pb, 700.)
   v_500p[f,:,:] = calc_var_at_plev(vc, p, pb, 500.)

   u_850_300p[f,:,:], v_850_300p[f,:,:] = calc_mean_wind_pres_nowgt(p, uc, vc, 850., 300.)

   u_dvg[f,:,:] = calc_ul_divergence(p, uc, vc, 400., 250.)
   v_dvg[f,:,:] = calc_ul_divergence(p, uc, vc, 400., 250.)

   shear_u_0to1[f,:,:], shear_v_0to1[f,:,:] = calc_wind_shear(z_agl, uc, vc, 0., 1000.)
   shear_u_0to6[f,:,:], shear_v_0to6[f,:,:] = calc_wind_shear(z_agl, uc, vc, 0., 6000.)

   bunk_r_u[f,:,:], bunk_r_v[f,:,:], bunk_l_u, bunk_l_v = calc_bunkers(p, z_agl, dz, uc, vc)

   corf_us_u[f,:,:], corf_us_v[f,:,:], corf_ds_u[f,:,:], corf_ds_v[f,:,:] = calc_corfidi(p, z_agl, uc, vc, 850., 300., 0., 1500.)

   srh_0to1[f,:,:] = calc_srh(z_agl, uc, vc, dz, 0., 1000., bunk_r_u[f,:,:], bunk_r_v[f,:,:])
   srh_0to3[f,:,:] = calc_srh(z_agl, uc, vc, dz, 0., 3000., bunk_r_u[f,:,:], bunk_r_v[f,:,:])

   stp_100mb[f,:,:] = calc_stp(cape_100mb[f,:,:], lcl_100mb[f,:,:], srh_0to1[f,:,:], shear_u_0to6[f,:,:], shear_v_0to6[f,:,:], cin_100mb[f,:,:])

######### Vertical velocity values #########

   w_700p[f,:,:] = calc_omega(p, wc, pb, temp, 700.)

######### Convert wind speeds to kts ##########

   u_10[f,:,:] = u_10[f,:,:] *  1.943844
   v_10[f,:,:] = v_10[f,:,:] *  1.943844
   u_500[f,:,:] = u_500[f,:,:] *  1.943844
   v_500[f,:,:] = v_500[f,:,:] *  1.943844
   u_850p[f,:,:] = u_850p[f,:,:] * 1.943844
   u_700p[f,:,:] = u_700p[f,:,:] * 1.943844
   u_500p[f,:,:] = u_500p[f,:,:] * 1.943844
   v_850p[f,:,:] = v_850p[f,:,:] * 1.943844
   v_700p[f,:,:] = v_700p[f,:,:] * 1.943844
   v_500p[f,:,:] = v_500p[f,:,:] * 1.943844
   u_850_300p[f,:,:] = u_850_300p[f,:,:] * 1.943844 # mean U 850-300 hPa flow[BCM]
   v_850_300p[f,:,:] = v_850_300p[f,:,:] * 1.943844 # mean V 850-300 hPa flow[BCM]
   u_bl[f,:,:] = u_bl[f,:,:] * 1.943844
   v_bl[f,:,:] = v_bl[f,:,:] * 1.943844
   u_dvg[f,:,:] = u_dvg[f,:,:] * 1.943844         #--- ULU divergence vectors[BCM]
   v_dvg[f,:,:] = v_dvg[f,:,:] * 1.943844         #--- ULV divergence vectors[BCM]
   shear_u_0to1[f,:,:] = shear_u_0to1[f,:,:] *  1.943844
   shear_v_0to1[f,:,:] = shear_v_0to1[f,:,:] *  1.943844
   shear_u_0to6[f,:,:] = shear_u_0to6[f,:,:] *  1.943844
   shear_v_0to6[f,:,:] = shear_v_0to6[f,:,:] *  1.943844
   bunk_r_u[f,:,:] = bunk_r_u[f,:,:] * 1.943844
   bunk_r_v[f,:,:] = bunk_r_v[f,:,:] * 1.943844
   corf_us_u[f,:,:] = corf_us_u[f,:,:] * 1.943844 #--- Corfidi USU vectors[BCM]
   corf_us_v[f,:,:] = corf_us_v[f,:,:] * 1.943844 #--- Corfidi USV vectors[BCM]
   corf_ds_u[f,:,:] = corf_ds_u[f,:,:] * 1.943844 #--- Corfidi DSU vectors[BCM]
   corf_ds_v[f,:,:] = corf_ds_v[f,:,:] * 1.943844 #--- Corfidi DSV vectors[BCM]

######### Storm-scale values #########

   dbz_col_max[f,:,:] = np.max(dbz, axis=0)


##################### Calc Ens. Mean Values: ########################

### Calc probability matched mean for reflectivity only ###

pmm_dz = dbz_col_max[0,:,:] * 0.
temp_dz = np.where(dbz_col_max > 100000., 0., dbz_col_max)
temp_mean_dz = np.mean(temp_dz, axis=0)

pmm_dz[:,:] = prob_match_mean(temp_dz[:,:,:], temp_mean_dz[:,:], neighborhood)

tf_2 = (t_2 - 273.15) * 1.8 + 32.     #convert to deg. F
tdf_2 = (td_2 - 273.15) * 1.8 + 32.     #convert to deg. F
tf_850p = (t_850p - 273.15) * 1.8 + 32.     #convert to deg. F
tdf_850p = (td_850p - 273.15) * 1.8 + 32.     #convert to deg. F
tf_700p = (t_700p - 273.15) * 1.8 + 32.     #convert to deg. F
tdf_700p = (td_700p - 273.15) * 1.8 + 32.     #convert to deg. F

#mean_cape = np.mean(cape_75mb, axis=0)
#mean_cape_0to3 = np.mean(cape_0to3_75mb, axis=0)
#mean_cin = np.mean(cin_75mb, axis=0)
#mean_lcl = np.mean(lcl_75mb, axis=0)
#mean_lfc = np.mean(lfc_75mb, axis=0)
#mean_stp = np.mean(stp_75mb, axis=0)

#mean_bunk_r_u = np.mean(bunk_r_u, axis=0)
#mean_bunk_r_v = np.mean(bunk_r_v, axis=0)
#mean_srh_0to1 = np.mean(srh_0to1, axis=0)
#mean_srh_0to3 = np.mean(srh_0to3, axis=0)
#mean_shear_u_0to6 = np.mean(shear_u_0to6, axis=0)
#mean_shear_v_0to6 = np.mean(shear_v_0to6, axis=0)
#mean_shear_u_0to1 = np.mean(shear_u_0to1, axis=0)
#mean_shear_v_0to1 = np.mean(shear_v_0to1, axis=0)

#mean_u_10 = np.mean(u_10, axis=0)
#mean_v_10 = np.mean(v_10, axis=0)
#mean_u_500 = np.mean(u_500, axis=0)
#mean_v_500 = np.mean(v_500, axis=0)

#mean_p_sfc = np.mean(p_sfc, axis=0)
#mean_t_2 = np.mean(t_2, axis=0)
#mean_td_2 = np.mean(td_2, axis=0)
#mean_the_ml = np.mean(th_e_75mb, axis=0)
#mean_thv_ml = np.mean(th_v_75mb, axis=0)
#mean_qv_2 = np.mean(qv_2, axis=0)

#mean_tf_2 = (mean_t_2 - 273.15) * 1.8 + 32.     #convert to deg. F
#mean_tdf_2 = (mean_td_2 - 273.15) * 1.8 + 32.     #convert to deg. F

#mean_sw_down = np.mean(sw_down, axis=0)

#mean_comp_dz = np.mean(dbz_col_max, axis=0)

##################### Write Summary File: ########################

##################### Get dimensions and attributes of WRFOUT file ########################

### Create file and dimensions: ###

try:
   print "Creating %s" % output_path
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)
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

### Ensemble mean variables ###

u_10_var = fout.createVariable('u_10', 'f4', ('NE','NY','NX',))
u_10_var.long_name = "U-component of 10-m wind"
u_10_var.units = "kts"

v_10_var = fout.createVariable('v_10', 'f4', ('NE','NY','NX',))
v_10_var.long_name = "V-component of 10-m wind"
v_10_var.units = "kts"

u_850p_var = fout.createVariable('u_850p', 'f4', ('NE','NY','NX',))
u_850p_var.long_name = "U-component of 850-hPa wind"
u_850p_var.units = "kts"

u_700p_var = fout.createVariable('u_700p', 'f4', ('NE','NY','NX',))
u_700p_var.long_name = "U-component of 700-hPa wind"
u_700p_var.units = "kts"

u_500p_var = fout.createVariable('u_500p', 'f4', ('NE','NY','NX',))
u_500p_var.long_name = "U-component of 500-hPa wind"
u_500p_var.units = "kts"

v_850p_var = fout.createVariable('v_850p', 'f4', ('NE','NY','NX',))
v_850p_var.long_name = "V-component of 850-hPa wind"
v_850p_var.units = "kts"

v_700p_var = fout.createVariable('v_700p', 'f4', ('NE','NY','NX',))
v_700p_var.long_name = "V-component of 700-hPa wind"
v_700p_var.units = "kts"

v_500p_var = fout.createVariable('v_500p', 'f4', ('NE','NY','NX',))
v_500p_var.long_name = "V-component of 500-hPa wind"
v_500p_var.units = "kts"

w_700p_var = fout.createVariable('w_700p', 'f4', ('NE','NY','NX',))
w_700p_var.long_name = "700-hPa omega vertical velocity"
w_700p_var.units = "pa/s"

u_850_300p_var = fout.createVariable('u_850_300p', 'f4', ('NE','NY','NX',))
u_850_300p_var.long_name = "U-component of 850-300 hPa wind"
u_850_300p_var.units = "kts"

v_850_300p_var = fout.createVariable('v_850_300p', 'f4', ('NE','NY','NX',))
v_850_300p_var.long_name = "V-component of 850-300 hPa wind"
v_850_300p_var.units = "kts"

u_bl_var = fout.createVariable('u_bl', 'f4', ('NE', 'NY', 'NX',))
u_bl_var.long_name = "U-component of 1000-925 hPa wind"
u_bl_var.units = "kts"

v_bl_var = fout.createVariable('v_bl', 'f4', ('NE', 'NY', 'NX',))
v_bl_var.long_name = "V-component of 1000-925 hPa wind"
v_bl_var.units = "kts"

p_sfc_var = fout.createVariable('p_sfc', 'f4', ('NE','NY','NX',))
p_sfc_var.long_name = "Surface Pressure"
p_sfc_var.units = "hPa"

#no pbl height in WRFWOF yet
#pbl_hgt_var = fout.createVariable('pbl_hgt', 'f4', ('NE','NY','NX',))
#pbl_hgt_var.long_name = "PBL Height"
#pbl_hgt_var.units = "m"

mslp_var = fout.createVariable('mslp', 'f4', ('NE', 'NY','NX',))
mslp_var.long_name = "Mean Sea-Level Pressure"
mslp_var.units = "hPa"

t_2_var = fout.createVariable('t_2', 'f4', ('NE','NY','NX',))
t_2_var.long_name = "2-m temperature"
t_2_var.units = "K"

t_850p_var = fout.createVariable('t_850p', 'f4', ('NE','NY','NX',))
t_850p_var.long_name = "850-hPa temperature"
t_850p_var.units = "K"

t_700p_var = fout.createVariable('t_700p', 'f4', ('NE','NY','NX',))
t_700p_var.long_name = "700-hPa temperature"
t_700p_var.units = "K"

td_2_var = fout.createVariable('td_2', 'f4', ('NE','NY','NX',))
td_2_var.long_name = "2-m dewpoint temperature"
td_2_var.units = "K"

td_850p_var = fout.createVariable('td_850p', 'f4', ('NE','NY','NX',))
td_850p_var.long_name = "850-hPa dewpoint temperature"
td_850p_var.units = "K"

td_700p_var = fout.createVariable('td_700p', 'f4', ('NE','NY','NX',))
td_700p_var.long_name = "700-hPa dewpoint temperature"
td_700p_var.units = "K"

qv_2_var = fout.createVariable('qv_2', 'f4', ('NE','NY','NX',))
qv_2_var.long_name = "2-m water vapor mixing ratio"
qv_2_var.units = "Kg/Kg"

pw_var = fout.createVariable('pw', 'f4', ('NE','NY','NX',))
pw_var.long_name = "Precipitable water"
pw_var.units = "in"

ctp_var = fout.createVariable('ctp', 'f4', ('NE','NY','NX',))
ctp_var.long_name = "Cloud top pressure"
ctp_var.units = "hPa"

th_v_2_var = fout.createVariable('th_v_ml', 'f4', ('NE','NY','NX',))
th_v_2_var.long_name = "100 hPa mixed layer virtual potential temperature"
th_v_2_var.units = "K"

th_e_2_var = fout.createVariable('th_e_ml', 'f4', ('NE','NY','NX',))
th_e_2_var.long_name = "100 hPa mixed layer equivalent potential temperature"
th_e_2_var.units = "K"

lcl_ml_var = fout.createVariable('lcl_ml', 'f4', ('NE','NY','NX',))
lcl_ml_var.long_name = "100 hPa mixed layer lifted condensation level height"
lcl_ml_var.units = "m"

lfc_ml_var = fout.createVariable('lfc_ml', 'f4', ('NE','NY','NX',))
lfc_ml_var.long_name = "100 hPa mixed layer level of free convection"
lfc_ml_var.units = "m"

cape_ml_var = fout.createVariable('cape_ml', 'f4', ('NE','NY','NX',))
cape_ml_var.long_name = "100 hPa mixed layer CAPE"
cape_ml_var.units = "J/Kg"

cin_ml_var = fout.createVariable('cin_ml', 'f4', ('NE','NY','NX',))
cin_ml_var.long_name = "100 hPa mixed layer convective inhibition"
cin_ml_var.units = "J/Kg"

cape_0to3_ml_var = fout.createVariable('cape_0to3_ml', 'f4', ('NE','NY','NX',))
cape_0to3_ml_var.long_name = "100 hPa mixed layer CAPE below 3 km"
cape_0to3_ml_var.units = "J/Kg"

pbl_mfc_var = fout.createVariable('pbl_mfc', 'f4', ('NE','NY','NX',))
pbl_mfc_var.long_name = "Boundary-layer moisture convergence"
pbl_mfc_var.units = "g Kg^-1 s^-1"

u_500_var = fout.createVariable('u_500', 'f4', ('NE','NY','NX',))
u_500_var.long_name = "U-component of 500-m wind"
u_500_var.units = "kts"

v_500_var = fout.createVariable('v_500', 'f4', ('NE','NY','NX',))
v_500_var.long_name = "V-component of 500-m wind"
v_500_var.units = "kts"

shear_u_0to1_var = fout.createVariable('shear_u_0to1', 'f4', ('NE','NY','NX',))
shear_u_0to1_var.long_name = "U-component of 0-1 km vertical wind shear"
shear_u_0to1_var.units = "kts"

shear_v_0to1_var = fout.createVariable('shear_v_0to1', 'f4', ('NE','NY','NX',))
shear_v_0to1_var.long_name = "V-component of 0-1 km vertical wind shear"
shear_v_0to1_var.units = "kts"

shear_u_0to6_var = fout.createVariable('shear_u_0to6', 'f4', ('NE','NY','NX',))
shear_u_0to6_var.long_name = "U-component of 0-6 km vertical wind shear"
shear_u_0to6_var.units = "kts"

shear_v_0to6_var = fout.createVariable('shear_v_0to6', 'f4', ('NE','NY','NX',))
shear_v_0to6_var.long_name = "V-component of 0-6 km vertical wind shear"
shear_v_0to6_var.units = "kts"

bunk_r_u_var = fout.createVariable('bunk_r_u', 'f4', ('NE','NY','NX',))
bunk_r_u_var.long_name = "U-component Bunkers right-mover storm motion"
bunk_r_u_var.units = "kts"

bunk_r_v_var = fout.createVariable('bunk_r_v', 'f4', ('NE','NY','NX',))
bunk_r_v_var.long_name = "V-component of Bunkers right-mover storm motion"
bunk_r_v_var.units = "kts"

corf_us_u_var = fout.createVariable('corf_us_u', 'f4', ('NE', 'NY', 'NX',))
corf_us_u_var.long_name = "U-component of Corfidi upshear vector"
corf_us_u_var.units = "kts"

corf_us_v_var = fout.createVariable('corf_us_v', 'f4', ('NE', 'NY', 'NX',))
corf_us_v_var.long_name = "V-component of Corfidi upshear vector"
corf_us_v_var.units = "kts"

corf_ds_u_var = fout.createVariable('corf_ds_u', 'f4', ('NE', 'NY', 'NX',))
corf_ds_u_var.long_name = "U-component of Corfidi downshear vector"
corf_ds_u_var.units = "kts"

corf_ds_v_var = fout.createVariable('corf_ds_v', 'f4', ('NE', 'NY', 'NX',))
corf_ds_v_var.long_name = "V-component of Corfidi downshear vector"
corf_ds_v_var.units = "kts"

u_dvg_var = fout.createVariable('u_dvg', 'f4', ('NE', 'NY', 'NX',))
u_dvg_var.long_name = "U-component of 400-250 hPa divergence vector"
u_dvg_var.units = "kts"

v_dvg_var = fout.createVariable('v_dvg', 'f4', ('NE', 'NY', 'NX',))
v_dvg_var.long_name = "V-component of 400-250 hPa divergence vector"
v_dvg_var.units = "kts"

srh_0to1_var = fout.createVariable('srh_0to1', 'f4', ('NE','NY','NX',))
srh_0to1_var.long_name = "0-1 km storm relative helicity"
srh_0to1_var.units = "m^2/s^2"

srh_0to3_var = fout.createVariable('srh_0to3', 'f4', ('NE','NY','NX',))
srh_0to3_var.long_name = "0-3 km storm relative helicity"
srh_0to3_var.units = "m^2/s^2"

stp_ml = fout.createVariable('stp_ml', 'f4', ('NE','NY','NX',))
stp_ml.long_name = "100-hPa mixed layer significant tornado parameter"
stp_ml.units = "unitless"

sw_down_var = fout.createVariable('sw_down', 'f4', ('NE','NY','NX',))
sw_down_var.long_name = "Downward flux of shortwave radiation at the ground"
sw_down_var.units = "W/m^2"

iwp_var = fout.createVariable('iwp', 'f4', ('NE','NY','NX',))
iwp_var.long_name = "Cloud ice water path"
iwp_var.units = ""

lwp_var = fout.createVariable('lwp', 'f4', ('NE','NY','NX',))
lwp_var.long_name = "Cloud liquid water path"
lwp_var.units = ""

ccomp_dz_var = fout.createVariable('comp_dz', 'f4', ('NE','NY','NX',))
ccomp_dz_var.long_name = "Mean composite reflectivity (93 km neighborhood)"
ccomp_dz_var.units = "dBZ"

ccomp_dz_pmm_var = fout.createVariable('comp_dz_pmm', 'f4', ('NY','NX',))
ccomp_dz_pmm_var.long_name = "Probability matched mean composite reflectivity (93 km neighborhood)"
ccomp_dz_pmm_var.units = "dBZ"

### Write variables ###

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:] = hgt

fout.variables['u_10'][:] = u_10
fout.variables['v_10'][:] = v_10

fout.variables['u_850p'][:] = u_850p    #---BCM
fout.variables['v_850p'][:] = v_850p    #---BCM
fout.variables['u_700p'][:] = u_700p    #---BCM
fout.variables['v_700p'][:] = v_700p    #---BCM
fout.variables['u_500p'][:] = u_500p    #---BCM
fout.variables['v_500p'][:] = v_500p    #---BCM
fout.variables['w_700p'][:] = w_700p    #---BCM

fout.variables['u_850_300p'][:] = u_850_300p    #---BCM
fout.variables['v_850_300p'][:] = v_850_300p    #---BCM

fout.variables['u_bl'][:] = u_bl        #---BCM
fout.variables['v_bl'][:] = v_bl        #---BCM

fout.variables['t_850p'][:] = tf_850p   #---BCM
fout.variables['td_850p'][:] = tdf_850p #---BCM

fout.variables['t_700p'][:] = tf_700p   #---BCM
fout.variables['td_700p'][:] = tdf_700p #---BCM

fout.variables['p_sfc'][:] = p_sfc
#fout.variables['pbl_hgt'][:] = pbl_hgt  #---BCM  #not yet
fout.variables['mslp'][:,:] = mslp
fout.variables['t_2'][:] = tf_2
fout.variables['td_2'][:] = tdf_2
fout.variables['qv_2'][:] = qv_2

fout.variables['pw'][:] = pw
fout.variables['ctp'][:] = ctp

fout.variables['th_v_ml'][:] = th_v_100mb
fout.variables['th_e_ml'][:] = th_e_100mb

fout.variables['lcl_ml'][:] = lcl_100mb
fout.variables['lfc_ml'][:] = lfc_100mb
fout.variables['cape_ml'][:] = cape_100mb
fout.variables['cin_ml'][:] = cin_100mb
fout.variables['cape_0to3_ml'][:] = cape_0to3_100mb
fout.variables['pbl_mfc'][:] = neigh_pbl_mfc  #---BCM

fout.variables['u_500'][:] = u_500
fout.variables['v_500'][:] = v_500
fout.variables['shear_u_0to1'][:] = shear_u_0to1
fout.variables['shear_v_0to1'][:] = shear_v_0to1
fout.variables['shear_u_0to6'][:] = shear_u_0to6
fout.variables['shear_v_0to6'][:] = shear_v_0to6
fout.variables['bunk_r_u'][:] = bunk_r_u
fout.variables['bunk_r_v'][:] = bunk_r_v
fout.variables['corf_us_u'][:] = corf_us_u #--- BCM
fout.variables['corf_us_v'][:] = corf_us_v #--- BCM
fout.variables['corf_ds_u'][:] = corf_ds_u #--- BCM
fout.variables['corf_ds_v'][:] = corf_ds_v #--- BCM
fout.variables['u_dvg'][:] = u_dvg         #--- BCM
fout.variables['v_dvg'][:] = v_dvg         #--- BCM
fout.variables['srh_0to1'][:] = srh_0to1
fout.variables['srh_0to3'][:] = srh_0to3
fout.variables['stp_ml'][:] = stp_100mb

fout.variables['comp_dz_pmm'][:] = pmm_dz

fout.variables['comp_dz'][:] = dbz_col_max

fout.variables['sw_down'][:] = sw_down
fout.variables['iwp'][:] = iwp
fout.variables['lwp'][:] = lwp

### Close output file ###

fout.close()
del fout




