#!/usr/bin/env python
#
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import netCDF4 as ncdf
#from optparse import OptionParser

import strmrpt

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Main function defined to return correct sys.exit() calls
#

#parser = OptionParser()
#
#parser.add_option("-d", "--date",    dest="date" )
#(options, args) = parser.parse_args()
if len(sys.argv) != 2:
  print "Usage: %s YYYYMMDD"% sys.argv[0]
  sys.exit(0)

dtime = sys.argv[1]
filename  = os.path.join("/scratch/ywang/EPIC/test_runs","%s18/mem_001/grid_spec.nc"%dtime)
filename2 = os.path.join("/scratch/ywang/EPIC/test_runs","%s18/mem_001/INPUT/C3337_grid.tile7.nc"%dtime)

figure = plt.figure(figsize = (12,12) )

print "filename  = %s"% filename
print "filename2 = %s"% filename2

f         = ncdf.Dataset(filename, "r")
glat      = f['grid_lat'][...]
glon      = f['grid_lon'][...]
glatt     = f['grid_latt'][...]
glont     = f['grid_lont'][...]

g = ncdf.Dataset(filename2, "r")
ctrlat = g.plat
ctrlon = g.plon
g.close()

print "ctrlat = %f, ctrlon = %f"% (ctrlat, ctrlon)

sw_lon = glon.min() - 360.
sw_lat = glat.min()
ne_lon = glon.max() - 360.
ne_lat = glat.max()
lat_0  = (ne_lat+sw_lat)/2.0    #f.target_lat
lon_0  = (ne_lon+sw_lon)/2.0    #f.target_lon

print "sw_lat = %f, sw_lon = %f" %(sw_lat, sw_lon)
print "ne_lat = %f, ne_lon = %f" %(ne_lat, ne_lon)

ny, nx = glat.shape

#-----------------------------------------------------------------------
#
# gnomomic grid
#
#-----------------------------------------------------------------------

map0 = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, \
             urcrnrlon=ne_lon, urcrnrlat=ne_lat, \
             lat_0=lat_0, lon_0=lon_0, \
             projection = 'gnom',      \
             resolution='c',   \
             area_thresh=5000.,
             suppress_ticks=True)

x, y = map0(glon, glat)
xc, yc = map0(glont, glatt)

xwidth = x.max() - x.min()
ywidth = y.max() - y.min()

# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS84 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
nx1 = 650   #1301
ny1 = 550    #921
dx1 = 3000.
dy1 = 3000.
xsize=(nx1-1)*dx1
ysize=(ny1-1)*dy1

map1 = Basemap(width=xsize,height=ysize,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=30.,lat_2=60.,lat_0=ctrlat,lon_0=ctrlon)

y1 = np.arange(map1.ymin,map1.ymax+1000,dy1)
print 'y1 = ',y1[0],y1[1],y1[2],map1.ymax, y1.shape
x1 = np.arange(map1.xmin,map1.xmax+1000,dx1)
print 'x1 = ',x1[0],x1[1],x1[2], map1.xmax,x1.shape
x1_2, y1_2 = np.meshgrid(x1,y1)

print x1_2.shape,y1_2.shape

lon1, lat1 = map1(x1_2,y1_2,inverse=True)


#-----------------------------------------------------------------------
#
# base grid to be plotted
#
#-----------------------------------------------------------------------
xwidth = xwidth + 500000
ywidth = ywidth + 500000

#print "width = ",xwidth, ywidth

mymap = Basemap(width=xwidth, height=ywidth, \
            lat_0=lat_0, lon_0=lon_0, \
            projection = 'gnom',      \
            resolution='c',   \
            area_thresh=5000.,
            suppress_ticks=True)


mymap.drawcoastlines(linewidth=0.1)
mymap.drawcountries(linewidth=0.1)
mymap.drawstates(linewidth=0.1)

parallels = np.arange(10.,60, 10.)
# labels = [left,right,top,bottom]
mymap.drawparallels(parallels,color='b',labels=[True,True,False,False])
meridians = np.arange(-130.,-60.,10.)
mymap.drawmeridians(meridians,color='b',labels=[False,False,False,True])

#-----------------------------------------------------------------------
#
# compute aspect
#
#-----------------------------------------------------------------------

x0, y0 = mymap(glon, glat)

dx = x0[:,:-1] - x0[:,1:]
dy = y0[:-1,:] - y0[1:,:]

print "dx = ", dx.min(), dx.max(), dx.shape
print "dy = ", dy.min(), dy.max(), dy.shape
#print x.min(), x.max()
#print y.min(), y.max()

aspect = dx[1:,:] / dy[:,1:]

print "aspect_max = %f at "%aspect.max(), np.unravel_index(aspect.argmax(), dims = aspect.shape)
print "aspect_min = %f at "%aspect.min(), np.unravel_index(aspect.argmin(), dims = aspect.shape)

xc,yc = mymap(ctrlon,ctrlat)
plt.text(xc,yc,'o',color='k')
plt.annotate('Center', xy=(xc,yc),  xycoords='data',
                xytext=(-12, -24), textcoords='offset points',
                color='k',
                arrowprops=dict(arrowstyle="->")
                )

mymap.contourf(x0[1:,1:], y0[1:,1:], aspect)
plt.colorbar(shrink=0.5)

#
# Plot original gnomonic grid
#
mymap.plot(x0[:,    0], y0[:,   0], color='g', linewidth=0.5)
mymap.plot(x0[:, nx-1], y0[:,nx-1], color='g', linewidth=0.5)
mymap.plot(x0[0,    :], y0[0,   :], color='g', linewidth=0.5)
mymap.plot(x0[ny-1, :], y0[ny-1,:], color='g', linewidth=0.5)
plt.text(x0[0,0],y0[0,0]-24,'fv3 grid',color='g')

#
# plot lambert
#
x1,y1 = mymap(lon1,lat1)

mymap.plot(x1[:,     0], y1[:,    0], color='r', linewidth=0.5)
mymap.plot(x1[:, nx1-1], y1[:,nx1-1], color='r', linewidth=0.5)
mymap.plot(x1[0,     :], y1[0,    :], color='r', linewidth=0.5)
mymap.plot(x1[ny1-1, :], y1[ny1-1,:], color='r', linewidth=0.5)
plt.text(x1[ny1-20,nx1-90],y1[ny1-20,nx1-90],'Lambert grid',color='r')

print 'lat1 = %f, lon1 = %f' %(lat1[0,0],lon1[0,0])
plt.annotate('SW corner', xy=(x1[0,0], y1[0,0]),  xycoords='data',
                xytext=(12, 12), textcoords='offset points',
                color='k',
                arrowprops=dict(arrowstyle="fancy", color='r')
                )

strmrpt.plot_storms(mymap,dtime,'/oldscratch/ywang/EPIC/Program/scripts')

plt.title("FV3 GRID on %s with grid aspect ratio"%dtime)
figname = "grid_%s.png"%dtime
print "Saving figure to %s ..." % figname
figure.savefig(figname, format='png')

#os.system("open %s" % newfilename)
# End of file
