#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from scipy import signal
from scipy import *
from scipy import ndimage
from skimage.morphology import label
from skimage.measure import regionprops
import math
from math import radians, tan, sin, cos, pi, atan, sqrt, pow, asin, acos
import pylab as P
import numpy as np
from numpy import NAN
import sys
import netCDF4
from optparse import OptionParser
import pickle as pl
from netcdftime import utime
import os
import time as timeit
import ctables
import news_e_post_cbook
from news_e_post_cbook import *
import news_e_plotting_cbook_v2
from news_e_plotting_cbook_v2 import *


shapefile = '/work1/jessica.choate/HWT_LSR/lsr_201705161900_201705170600'
exp_file = '/work1/jessica.choate/summary_files/20170516_rlt_odin/2300/2017-05-16_23:00:00_01.nc'

edge = 7

try:
    fin = netCDF4.Dataset(exp_file, "r")
    print "Opening %s \n" % exp_file
except:
    print "%s does not exist! \n" % exp_file
    sys.exit(1)

xlat = fin.variables['XLAT'][:]
xlon = fin.variables['XLON'][:]

xlat = xlat[edge:-edge,edge:-edge]
xlon = xlon[edge:-edge,edge:-edge]

sw_lat = xlat[0,0]
sw_lon = xlon[0,0]
ne_lat = xlat[-1,-1]
ne_lon = xlon[-1,-1]

cen_lat = fin.CEN_LAT
cen_lon = fin.CEN_LON
stand_lon = fin.STAND_LON
true_lat1 = fin.TRUE_LAT1
true_lat2 = fin.TRUE_LAT2

fin.close()
del fin

map = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, urcrnrlon=ne_lon, urcrnrlat=ne_lat, projection='lcc', lat_1=true_lat1, lat_2=true_lat2, lat_0=cen_lat, lon_0=cen_lon, resolution = 'c', area_thresh = 10000.)

map.readshapefile(shapefile, 'lsr', drawbounds = False) #read shapefile

hail_21 = []
wind_21 = []
tornado_21 = []

hail_22 = []
wind_22 = []
tornado_22 = []

t_1 = 21 * 3600.
t_2 = 22 * 3600. 
t_3 = 23 * 3600. 

for info, shape in zip(map.lsr_info, map.lsr):
    temp_inithr = info['VALID'][8:10]
    temp_initmin = info['VALID'][10:12]
    init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.
    
    temp_lat = double(info['LAT'])
    temp_lon = double(info['LON'])
    
    if ((init_sec >= t_1) and (init_sec <= t_2) and (temp_lat >= sw_lat) and (temp_lat <= ne_lat) and (temp_lon >= sw_lon) and (temp_lon <= ne_lon)): 
        if (info['TYPECODE'] == 'H'):
            hail_21.append(info['TYPECODE'])
        elif (info['TYPECODE'] == 'D'):
            wind_21.append(info['TYPECODE'])
        elif (info['TYPECODE'] == 'T'):
            tornado_21.append(info['TYPECODE'])

    if ((init_sec >= t_2) and (init_sec <= t_3) and (temp_lat >= sw_lat) and (temp_lat <= ne_lat) and (temp_lon >= sw_lon) and (temp_lon <= ne_lon)): 
        if (info['TYPECODE'] == 'H'):
            hail_22.append(info['TYPECODE'])
        elif (info['TYPECODE'] == 'D'):
            wind_22.append(info['TYPECODE'])
        elif (info['TYPECODE'] == 'T'):
            tornado_22.append(info['TYPECODE'])

print 'Hail 21-22: ', len(hail_21), 'Wind 21-22: ', len(wind_21), 'Tor 21-22: ', len(tornado_21)
print 'Hail 22-23: ', len(hail_22), 'Wind 22-23: ', len(wind_22), 'Tor 22-23: ', len(tornado_22)




