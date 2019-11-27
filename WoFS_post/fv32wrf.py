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
from datetime import datetime, timedelta, date, time
from optparse import OptionParser
from news_e_post_cbook import *

############################### File Variables: ##########################################

parser = OptionParser()
parser.add_option("-d", dest="indir", type="string", default= None, help="Input directory containing FV3 file to process")
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

print "TIME : %d \n" % t

dt = t * 10

fcsttime = timedelta(minutes=dt)
total_seconds = int(fcsttime.total_seconds())
hours, remainder = divmod(total_seconds,60*60)
minutes, seconds = divmod(remainder,60)
fcststr='%03d:%02d:%02d'%(hours,minutes,seconds)


#gridfile = indir + "/grid_spec.nc"
dynfile = indir + "/dynf%s.nc"%fcststr
phyfile = indir + "/phyf%s.nc"%fcststr

##################### Process FV3HISTORY File ########################

#try:                                                 #open grid file
#  fgrid = netCDF4.Dataset(gridfile, "r")
#  print "Opening %s \n" % gridfile
#except:
#  print "%s does not exist! \n" % gridfile
#  sys.exit(1)

try:                                                 #open 3d file
  fdyn = netCDF4.Dataset(dynfile, "r")
  print "Opening %s \n" % dynfile
except:
  print "%s does not exist! \n" % dynfile
  sys.exit(1)

try:                                                 #open vertical coordinate file
  fphy = netCDF4.Dataset(phyfile, "r")
  print "Opening %s \n" % phyfile
except:
  print "%s does not exist! \n" % phyfile
  sys.exit(1)

################# Get dimensions and attributes using first file ####################

### Get init and valid times ###

time_var = fdyn.variables["time"]                             #Read title time string
start_date=time_var.getncattr("units")
init_year = start_date[-19:-15]
init_mon = start_date[-14:-12]
init_day = start_date[-11:-9]
init_hr = start_date[-8:-6]
init_min = start_date[-5:-3]
start_date_true = init_year + init_mon + init_day + init_hr + init_min

init_date = init_year + init_mon + init_day               #YYYYMMDD string for output file
init_time = init_hr + init_min                            #HHMM string for initialization time

valid_time_dt = datetime(int(init_year),int(init_mon), int(init_day), int(init_hr), int(init_min)) + timedelta(minutes=dt)
valid_time_tt = valid_time_dt.timetuple()
valid_time = '{:>02d}'.format(valid_time_tt[3]) + '{:>02d}'.format(valid_time_tt[4])              #HHMM string for valid time

print "valid_time = %s" % valid_time

valid_hr = str(valid_time_tt[3])
valid_min = str(valid_time_tt[4])

### Set output path ###
nt = t
timestep = "%02d"%nt

outname = "FV3_ENS_" + timestep + "_" + init_date + "_" + init_time + "_" + valid_time + ".nc"         #output file
output_path = os.path.join(outdir, outname)

### Get grid/projection info ###

dx = fdyn.dx                       #east-west grid spacing (m)
dy = fdyn.dy                       #north-south grid spacing (m)
cen_lat    = fdyn.cen_lat   #39.71457                                   #center of domain latitude (dec deg)
cen_lon    = fdyn.cen_lon   #-96.44037                                   #center of domain longitude (dec deg)
stand_lon  = fdyn.cen_lon   #-96.44036                                 #center lon of Lambert conformal projection
true_lat_1 = fdyn.stdlat1   #30.0                                #true lat value 1 for Lambert conformal conversion (dec deg)
true_lat_2 = fdyn.stdlat2   #60.0                                #true lat value 2 for Lambert conformal conversion (dec deg)


xlat = np.degrees(fdyn.variables["grid_yt"][:,:])       #latitude (dec deg; Lambert conformal)
xlon = np.degrees(fdyn.variables["grid_xt"][:,:])       #longitude (dec deg; Lambert conformal)

ny = xlat.shape[0]
nx = xlat.shape[1]
nz = fdyn.dimensions['pfull'].size

### Calculate initial and valid time in seconds ###

init_time_seconds = int(init_hr) * 3600. + int(init_min) * 60.
valid_time_seconds = int(valid_hr) * 3600. + int(valid_min) * 60.

if (valid_time_seconds < 43200.):         #Convert values past 0000 UTC, assumes forecast not run past 12 UTC
  valid_time_seconds = valid_time_seconds + 86400.

############################### Read WRFOUT variables: ##################################

#wspd80= np.hypot(fin.variables["ugrd100m"][t,:,:], fin.variables["vgrd100m"][t,:,:])
t = 0

w_up = fdyn.variables["upvvelmax"][t,:,:]

rain_temp = fphy.variables["prateb_ave"][t,:,:]                  #5-min average precip rate (kg/m^2/s)
rain = rain_temp / 1000.0 * 300.0 * 1000.0                     #convert to mm

#   soil_moisture[f,:,:] = fin.variables["SMOIS"][0,0,:,:]

#hail_temp = fin.variables["HAIL_MAXK1"][t,:,:]
hail = np.zeros((ny,nx))

#hailcast_temp = fin.variables["HAILCAST_DIAM_MAX"][t,:,:]
hailcast = np.zeros((ny,nx))                      #convert to inches

hgt = fdyn.variables["hgtsfc"][0,:,:]          #terrain height above MSL (m)
print "hgt shape = ", np.shape(hgt)

u = np.flip(fdyn.variables["ugrd"][t,:,:,:],0)                       #expects var dimensions of (nt, nz, ny, nx) with nt = 1
v = np.flip(fdyn.variables["vgrd"][t,:,:,:],0)
w = np.flip(fdyn.variables["dzdt"][t,:,:,:],0)
#rel_vort = np.flip(fin.variables["rel_vort"][t,:,:,:],0)
tk = np.flip(fdyn.variables["tmp"][t,:,:,:],0)
print 'Min max tk = %f  %f' % (np.min(tk), np.max(tk))
qv = np.flip(fdyn.variables["spfh"][t,:,:,:],0)
qc = np.flip(fdyn.variables["clwmr"][t,:,:,:],0)
qr = np.flip(fdyn.variables["rwmr"][t,:,:,:],0)
qi = np.flip(fdyn.variables["icmr"][t,:,:,:],0)
qs = np.flip(fdyn.variables["snmr"][t,:,:,:],0)
qg = np.flip(fdyn.variables["grle"][t,:,:,:],0)

po = 1.0E5
To = 300.0
A = 50.0
Rd = 287.0

ak = np.flip(fdyn.ak)   #np.flip(fcoord.variables["hyam"],0)/1.0E-5
bk = np.flip(fdyn.bk)   #np.flip(fcoord.variables["hybm"],0)

pdhs = po*np.exp((-To/A)+np.sqrt((To/A)**2-2*hgt*9.81/(A*Rd)))

pb = np.zeros((nz,ny,nx))

for i in arange(nz):
  pb[i,...] = ak[i]+bk[i]*pdhs

#pb = np.flip(fin.variables["pres_hyd"][t,:,:,:],0)
pfull = np.flip(fdyn.variables["pres"][t,:,:,:],0)
p = pfull-pb
th = (100000.00/pfull) ** (287/1004.5) * tk - 300.0

print 'Min max th = %f  %f' % (np.min(th), np.max(th))
if (np.max(tk) > 350.0) :
  sys.exit(1)


#nz = pb.shape[0]
#p = np.zeros((nz,ny,nx))

delz = np.flip(fdyn.variables["delz"][t,:,:,:],0)

phb = np.cumsum(delz[:,:,:], axis=0)
ph = np.zeros((nz,ny,nx))

dbz = np.flip(fphy.variables["refl_10cm"][t,:,:,:],0)
#wz = fin.variables["REL_VORT"][t,:,:,:],0)

u10 = fphy.variables["ugrd10m"][t,:,:]
v10 =fphy.variables["vgrd10m"][t,:,:]

q2 = fphy.variables["spfh2m"][t,:,:]
t2 = fphy.variables["tmp2m"][t,:,:]
psfc = fphy.variables["pressfc"][t,:,:]
th2 = (100000.00/psfc) ** (287/1004.5) * t2

#uh25 = fin.variables["uhmax25"][t,:,:]
#uh02 = fin.variables["uhmax03"][t,:,:]
#wz02 = fin.variables["maxvort02"][t,:,:]
uh02 = fdyn.variables["uhmax03"][t,:,:]   #calc_uh(w,rel_vort,phb,delz,0.,2000.)
uh25 = fdyn.variables["uhmax25"][t,:,:]   #calc_uh(w,rel_vort,phb,delz,2000.,5000.)
wz02 = fdyn.variables["maxvort02"][t,:,:] #calc_avg_vort(rel_vort,phb,delz,0.,2000.)


swd = fphy.variables["dswrf_ave"][t,:,:]
iwp = np.zeros((ny,nx))
lwp = np.zeros((ny,nx))

### Close WRFOUT file ###

fdyn.close()
fphy.close()
del fdyn, fphy

### Create file and dimensions: ###

try:
   print "Creating %s ...." % output_path
   fout = netCDF4.Dataset(output_path, "w")
except:
   print "Could not create %s!\n" % output_path

fout.createDimension('NX', nx)
fout.createDimension('NY', ny)
fout.createDimension('NZ', nz)

### Set Attributes: ###

setattr(fout,'DX',dx)
setattr(fout,'DY',dy)
setattr(fout,'CEN_LAT',cen_lat)
setattr(fout,'CEN_LON',cen_lon)
setattr(fout,'STAND_LON',stand_lon)
setattr(fout,'TRUELAT1',true_lat_1)
setattr(fout,'TRUELAT2',true_lat_2)
setattr(fout,'PROJECTION','Lambert Conformal')
setattr(fout,'START_DATE',start_date_true)
setattr(fout,'INIT_TIME_SECONDS',init_time_seconds)
setattr(fout,'VALID_TIME_SECONDS',valid_time_seconds)
setattr(fout,'FORECAST_TIME_STEP',nt)
print 'Create variables'

xlat1 = fout.createVariable('xlat', 'f4', ('NY','NX',))
xlat1.long_name = "Latitude"
xlat1.units = "degrees North"

xlon1 = fout.createVariable('xlon', 'f4', ('NY','NX',))
xlon1.long_name = "Longitude"
xlon1.units = "degrees West"

hgt1 = fout.createVariable('hgt', 'f4', ('NY','NX',))
hgt1.long_name = "Height AGL"
hgt1.units = "m"

ppsfc = fout.createVariable('psfc', 'f4', ('NY','NX',))
ppsfc.long_name = "Surface Pressure"
ppsfc.units = "Pa"

uu10 = fout.createVariable('u10', 'f4', ('NY','NX',))
uu10.long_name = "10-m Zonal Wind"
uu10.units = "m/s"

vv10 = fout.createVariable('v10', 'f4', ('NY','NX',))
vv10.long_name = "10-m Meridional Wind"
vv10.units = "m/s"

tt2 = fout.createVariable('t2', 'f4', ('NY','NX',))
tt2.long_name = "2-m Temperature"
tt2.units = "K"

tth2 = fout.createVariable('th2', 'f4', ('NY','NX',))
tth2.long_name = "2-m Potential Temperature"
tth2.units = "K"

qq2 = fout.createVariable('q2', 'f4', ('NY','NX',))
qq2.long_name = "2-m Water Vapor Mixing Ratio"
qq2.units = "kg/kg"

#wwspd80 = fout.createVariable('wspd80', 'f4', ('NY','NX',))
#wwspd80.long_name = "Maximum 100-m Wind Speed"
#wwspd80.units = "m/s"

wwup = fout.createVariable('w_up_max', 'f4', ('NY','NX',))
wwup.long_name = "Maximum Upwards Vertical Wind Speed"
wwup.units = "m/s"

uuh02 = fout.createVariable('uh_02_max', 'f4', ('NY','NX',))
uuh02.long_name = "Maximum 0-3 km Updraft Helicity"
uuh02.units = "m^2 s^-2"

uuh25 = fout.createVariable('uh_25_max', 'f4', ('NY','NX',))
uuh25.long_name = "Maximum 2-5 km Updraft Helicity"
uuh25.units = "m^2 s^-2"

wwz02 = fout.createVariable('wz_02_max', 'f4', ('NY','NX',))
wwz02.long_name = "Maximum 0-2 km Vertical Vorticity"
wwz02.units = "s^-1"

pprcp = fout.createVariable('prec_acc_nc', 'f4', ('NY','NX',))
pprcp.long_name = "5-minutes Accumulated Precipitation"
pprcp.units = "mm"

hhailmax = fout.createVariable('hail_maxk1', 'f4', ('NY','NX',))
hhailmax.long_name = "Max Hail Diameter K=1 (not used)"
hhailmax.units = "m"

hhailcastmax = fout.createVariable('hailcast_diam_max', 'f4', ('NY','NX',))
hhailcastmax.long_name = "Hailcast Max Hail Diameter (not used)"
hhailcastmax.units = "mm"

sswdown = fout.createVariable('swdown', 'f4', ('NY','NX',))
sswdown.long_name = "Downward Short Wave Flux at Ground Surface"
sswdown.units = "W/m^2"

iiwpath = fout.createVariable('iwp', 'f4', ('NY','NX',))
iiwpath.long_name = "Ice Water Path (not used)"
iiwpath.units = "g/m^2"

llwpath = fout.createVariable('lwp', 'f4', ('NY','NX',))
llwpath.long_name = "Liquid Water Path (not used)"
llwpath.units = "g/m^2"

rradar = fout.createVariable('refl_10cm', 'f4', ('NZ','NY','NX',))
rradar.long_name = "Simulated Radar Reflectivity"
rradar.units = "dBz"

uu = fout.createVariable('u', 'f4', ('NZ','NY','NX',))
uu.long_name = "Zonal Wind Speed"
uu.units = "m/s"

vv = fout.createVariable('v', 'f4', ('NZ','NY','NX',))
vv.long_name = "Meridional Wind Speed"
vv.units = "m/s"

ww = fout.createVariable('w', 'f4', ('NZ','NY','NX',))
ww.long_name = "Vertical Velocity"
ww.units = "m/s"

tt = fout.createVariable('t', 'f4', ('NZ','NY','NX',))
tt.long_name = "Perturbation Potential Temperature"
tt.units = "K"

pp = fout.createVariable('p','f4', ('NZ','NY','NX'))
pp.long_name = "Perturbation Pressure"
pp.units = "Pa"

ppb = fout.createVariable('pb','f4', ('NZ','NY','NX'))
ppb.long_name = "Base State Pressure"
ppb.units = "Pa"

ddz = fout.createVariable('delz','f4', ('NZ','NY','NX'))
ddz.long_name = "Centered Level Thickness"
ddz.units = "m"

ddz = fout.createVariable('phb','f4', ('NZ','NY','NX'))
ddz.long_name = "Base State Potential Height"
ddz.units = "gpm"

ddz = fout.createVariable('ph','f4', ('NZ','NY','NX'))
ddz.long_name = "Pertubation Potential Height (=0)"
ddz.units = "gpm"

#rrv = fout.createVariable('rel_vort','f4', ('NZ','NY','NX'))
#rrv.log_name = "Relative Vorticity"
#rrv.units = "1/s"

qqv = fout.createVariable('qvapor','f4', ('NZ','NY','NX'))
qqv.log_name = "Water Vapor Mixing Ratio"
qqv.units = "kg/m^2"

qqc = fout.createVariable('qcloud','f4', ('NZ','NY','NX'))
qqc.log_name = "Cloud Water Mixing Ratio"
qqc.units = "kg/m^2"

qqr = fout.createVariable('qrain','f4', ('NZ','NY','NX'))
qqr.log_name = "Rain Water Mixing Ratio"
qqr.units = "kg/m^2"

qqi = fout.createVariable('qice','f4', ('NZ','NY','NX'))
qqi.log_name = "Ice Mixing Ratio"
qqi.units = "kg/m^2"

qqs = fout.createVariable('qsnow','f4', ('NZ','NY','NX'))
qqs.log_name = "Snow Mixing Ratio"
qqs.units = "kg/m^2"

qqg = fout.createVariable('qgraupel','f4', ('NZ','NY','NX'))
qqg.log_name = "Graupel Mixing Ratio"
qqg.units = "kg/m^2"

print 'Write variables'

fout.variables['xlat'][:] = xlat
fout.variables['xlon'][:] = xlon
fout.variables['hgt'][:]  = hgt
fout.variables['psfc'][:] = psfc
fout.variables['u10'][:]  = u10
fout.variables['v10'][:]  = v10
fout.variables['t2'][:]   = t2
fout.variables['th2'][:]  = th2
fout.variables['q2'][:]   = q2
#fout.variables['wspd80'][:] = wspd80
fout.variables['w_up_max'][:]  = w_up
fout.variables['uh_25_max'][:] = uh25
fout.variables['uh_02_max'][:] = uh02
fout.variables['wz_02_max'][:] = wz02
fout.variables['refl_10cm'][:] = dbz
fout.variables['prec_acc_nc'][:] = rain
fout.variables['hail_maxk1'][:] = hail
fout.variables['hailcast_diam_max'][:] = hailcast
fout.variables['swdown'][:] = swd
fout.variables['iwp'][:] = iwp
fout.variables['lwp'][:] = lwp

fout.variables['u'][:] = u
fout.variables['v'][:] = v
fout.variables['w'][:] = w
fout.variables['p'][:] = p
fout.variables['pb'][:] = pb
fout.variables['delz'][:] = delz
fout.variables['ph'][:] = ph
fout.variables['phb'][:] = phb
fout.variables['t'][:] = th
#fout.variables['rel_vort'][:] = rel_vort
fout.variables['qvapor'][:] = qv
fout.variables['qcloud'][:] = qc
fout.variables['qrain'][:] = qr
fout.variables['qice'][:]= qi
fout.variables['qsnow'][:] = qs
fout.variables['qgraupel'][:] = qg

### Close output file ###

fout.close()
del fout
