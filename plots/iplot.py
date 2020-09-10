#!/usr/bin/env python3

import math
import os
import sys
from argparse import Namespace
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from netCDF4 import Dataset
from matplotlib.colors import LinearSegmentedColormap
from datetime import date, datetime, timedelta
import glob
from pathlib import Path
import fnmatch

initdate='2019052018'
inithr=initdate[8:10]
fhr=6

# Provide a file path to a forecast directory.
# The example below creates a dictionary containing 2 experiments, expt, through the first 4 forecast hours (including 0)
exptlist=['mem_001', 'mem_002']

file_path=os.path.join('/scratch/ywang/EPIC/test_runs',f'{initdate}')

# whether the two runs are using the same vertical coordinates,
# this decides whether to plot fields at different levels where they have max values
#diff_vert_lev=True
diff_vert_lev=False

folder = {expt: os.path.join(file_path,expt) for expt in exptlist }
filecount1 = min(len(fnmatch.filter(os.listdir(folder[exptlist[0]]), 'phyf*.nc')),len(fnmatch.filter(os.listdir(folder[exptlist[0]]), 'dynf*.nc')))
filecount2 = min(len(fnmatch.filter(os.listdir(folder[exptlist[1]]), 'phyf*.nc')),len(fnmatch.filter(os.listdir(folder[exptlist[1]]), 'dynf*.nc')))
last_fhr = min(filecount1, filecount2) - 1
print('available forecast hours=', last_fhr)
print('plotting forecast hour=', fhr)

# check if fhr is out of range of the forecast files
if (fhr > last_fhr):
	print("forecast hour out of range, use the last fhr")
	fhr = last_fhr

files = {expt: {x: [os.path.join(file_path, expt, x + f'001:00:00.nc')
         for i in range(last_fhr + 1)]  for x in ['dynf', 'phyf']} for expt in exptlist }

start_lat, start_lon=(34,-98)
end_lat, end_lon=(40, -88)
numpoints=100

def get_validdate(initdate, fhr):
    # In this section: initdate + fhr = validdate
    validdate=(datetime.strptime(initdate,"%Y%m%d%H")+timedelta(hours=fhr)).strftime("%Y%m%d%H")
    validday=validdate[0:8]
    validhr=validdate[8:10]
    return validdate, validday, validhr

validdate, validday, validhr=get_validdate(initdate, fhr)

# Load each file in the files dict into a NetCDF Dataset

def load_Dataset(files: dict, ds: dict):
    assert(isinstance(files, dict))
    for k, v in files.items():
        if isinstance(v, dict):
            ds[k] = {}
            load_Dataset(v, ds[k])
        else:
            ds[k] = [Dataset(f, 'r') for f in list(v)]

ds_dict = {}
load_Dataset(files, ds_dict)

# Load the dictionary into a Namespace data structure.
# This step is not necessary, but cuts down the syntax needed to reference each item in the dict.
#
# Example: Retrieve the 0 hr forecast Dataset from GFS Dynamics
#            dict: ds_dict['GFS']['dynf'][0]
#       Namespace: datasets.GFS.dynf[0]


def make_namespace(ns: Namespace(), d: dict):
    assert(isinstance(d, dict))
    for k, v in d.items():
        if isinstance(v, dict):
            leaf_ns = Namespace()
            ns.__dict__[k] = leaf_ns
            make_namespace(leaf_ns, v)
        else:
            ns.__dict__[k] = v


datasets = Namespace()

make_namespace(datasets, ds_dict)

#print(ds_dict[exptlist[0]]['dynf'][0])

# print('~~~~~~~~~ DYNAM FILE from GFS ~~~~~~~~~~~~~~~')
# # for v, info in datasets.GFS.dynf[0].variables.items():
# for v, info in ds_dict[exptlist[0]]['dynf'][0].variables.items():

def plot_data_zoom(dataL, dataR, start_lat, start_lon, end_lat, end_lon, title, expt, zmax=None):

    '''
    Input parameters:

        data: 2D Numpy array to be plotted in Left column
        lat: 2D Numpy array of latitude
        lon: 2D Numpy array of longitude
        title: String describing the variable being plotted.

    Draws a Basemap representation with the contoured data overlayed,
    with a colorbar for each experiment, and the difference between the two.

    '''

    def trim_grid():
        '''
        The u, v, and H data from analysis are all on grids either one column, or one row smaller than lat/lon.
        Return the smaller lat, lon grids, given the shape of the data to be plotted.
        Has no effect when all grids are the same size.
        '''
        y, x = np.shape(dataL)
        return lat[:y, :x], lon[:y, :x]

    def eq_contours(indata, indata2=None):

        '''
        Returns a balanced set of contours for data that has negative values.
        Also returns default colorbar to use for balanced, vs all positive values.
        '''
        cmap=plt.cm.get_cmap('bwr',20)
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmaplist[9] = [1,1,1,1]
        cmaplist[10] = [1,1,1,1]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

        minval = np.amin(indata) if indata2 is None else min(np.amin(indata), np.amin(indata2))
        maxval = np.amax(indata) if indata2 is None else max(np.amax(indata), np.amax(indata2))

        if minval == maxval:
            return np.linspace(-1, 1, 11), 'seismic'
        if np.amin(indata) < 0:
        # Set balanced contours. Choose an odd number in linspace below
            maxval = max(abs(minval), abs(maxval))
            return np.linspace(-maxval, maxval, 21), cmap
        else:
            c = np.linspace(minval, maxval, 21)
            if c[0] == 0:
                c[0] = c[1]/10
            return c, 'jet'

    lat_trim, lon_trim = trim_grid()

    #  add two more argumeents (fig, ax) to def plot_data(...)
    #     fig, ax = plt.subplots(figsize=(24, 12))
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    for i, data in enumerate([dataL, dataR]):
        m = Basemap(projection='mill',
                    llcrnrlon=360+start_lon-5,
                    urcrnrlon=360+end_lon+5 ,
                    llcrnrlat=start_lat-5,
                    urcrnrlat=end_lat+5,
                    resolution='c',
                    ax=ax[i],
                   )
        x, y = m(lon_trim, lat_trim)

        xlon=[360+start_lon, 360+end_lon]
        xlat=[start_lat, end_lat]
        xx, yy =m(xlon, xlat)

        # Use the same contour values for both experiments.
        if i < 2:
#         if i < 0:
            contours, cm = eq_contours(dataL, dataR)
        else:
            contours, cm = eq_contours(data)

        # Draw the contoured data over the map
        cs = m.contourf(x, y, data, contours, cmap=cm, ax=ax[i])
        m.plot(xx, yy, color='black', linewidth=3, linestyle='--', ax=ax[i])
        m.drawstates();
        m.drawcoastlines();
        m.drawmapboundary();
        m.drawparallels(np.arange(-90.,120.,5),labels=[1,0,0,0]);
        m.drawmeridians(np.arange(-180.,180.,10),labels=[0,0,0,1]);
        fig.colorbar(cs, ax=ax[i], orientation='vertical', shrink=0.7);
        if zmax is None:
            ax[i].set_title(f"{expt[i]}:\n {title}")
        else:
            ax[i].set_title(f"{expt[i]}:\n {title} at level {zmax[i]}")
    plt.show()

def latlon_to_xy(target_lat, target_lon):
    y, x = np.unravel_index((np.abs(lat - target_lat) + np.abs(lon - target_lon-360)).argmin(), lat.shape)
    return y,x

# function to derive values on the cross plane (3D to 2D) and make plots
def plot_vertical_cross(dataL, dataR, start_lat, start_lon, end_lat, end_lon, numpoints, title, expt):

# get list of x,y points from start and end lat/lon
    def get_xy_list(start_lat, start_lon, end_lat, end_lon, numpoints):
        lat_dist=end_lat - start_lat
        lon_dist=end_lon - start_lon
        lat_list=np.zeros(numpoints)
        lon_list=np.zeros(numpoints)
        y_list=np.zeros(numpoints)
        x_list=np.zeros(numpoints)
        for i in range(numpoints):
            lat_list[i]=start_lat + i * lat_dist / numpoints
            lon_list[i]=start_lon + i * lon_dist / numpoints
            y_list[i], x_list[i]=latlon_to_xy(lat_list[i], lon_list[i])
        return y_list, x_list

    y_list, x_list=get_xy_list(start_lat, start_lon, end_lat, end_lon, numpoints)

    def eq_contours(indata, indata2=None):

        '''
        Returns a balanced set of contours for data that has negative values.
        Also returns default colorbar to use for balanced, vs all positive values.
        '''
        cmap=plt.cm.get_cmap('bwr',20)
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmaplist[9] = [1,1,1,1]
        cmaplist[10] = [1,1,1,1]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

        minval = np.amin(indata) if indata2 is None  else min(np.amin(indata), np.amin(indata2))
        maxval = np.amax(indata) if indata2 is None  else max(np.amax(indata), np.amax(indata2))

        if minval == maxval:
            return np.linspace(-1, 1, 11), 'seismic'
        if np.amin(indata) < 0:
        # Set balanced contours. Choose an odd number in linspace below
            maxval = max(abs(minval), abs(maxval))
            return np.linspace(-maxval, maxval, 21), cmap
        else:
            c = np.linspace(minval, maxval, 21)
            if c[0] == 0:
                c[0] = c[1]/10
            return c, 'jet'

    # derive data on the cross section points
    def var3dto2d(data, y_list, x_list, numpoints):
        nlev=np.shape(data)[0]
        data2d=np.zeros((nlev,numpoints))
        for i in range(0, numpoints-1):
            iy=int(y_list[i])
            ix=int(x_list[i])
            data2d[:, i]=data[:, iy, ix]
        return data2d

    data2dL=var3dto2d(dataL, y_list, x_list, numpoints)
    data2dR=var3dto2d(dataR, y_list, x_list, numpoints)

# making vertical cross section plots
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    for i, data in enumerate([data2dL, data2dR]):

        x=np.shape(data)[1]
        y=np.shape(data)[0]

        if i < 2:
#         if i < 0:
             contours, cm = eq_contours(data2dL, data2dR)
        else:
             contours, cm = eq_contours(data)

        cs=ax[i].contourf(range(x),range(y)[::-1], data, contours, cmap=cm)
        ax[i].set_title(f"{expt[i]}: Cross Section \n {title}")
        plt.colorbar(cs, ax=ax[i], orientation='vertical', shrink=0.9)
        ax[i].grid()

    for ax in ax.flat:
        ax.set(xlabel='Grid Points', ylabel='Vertical Levels')

print(fhr)

lat = ds_dict[exptlist[0]]['dynf'][fhr]['grid_yt'][::] * 180 / math.pi
lon = ds_dict[exptlist[0]]['dynf'][fhr]['grid_xt'][::] * 180 / math.pi

print(np.shape(lat))
print(np.shape(lon))
data=ds_dict[exptlist[0]]['dynf'][fhr]['ugrd'][::]
print(np.shape(data))

sys.exit(0)
# Variables to plot from dynf and phyf files
vars_ = {
#     'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
    'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
#     'phyf': ['ugrd10m', 'vgrd10m', 'f10m', 'tmpsfc', 'tmp2m', 'spfh2m', 'pwatclm', 'tprcp', 'prate_ave', 'soilm', 'soilt1', 'pressfc', 'albdo_ave',  'lhtfl_ave', 'shtfl_ave', 'hpbl' ],

}
### find the level where variable field has max values
def find_lev_max(var3d):
    nlev=np.shape(var3d)[0]
    varmax=max(np.amax(var3d[i,:,:]) for i  in range(nlev))
#     print(varmax)
    for i in range(nlev):
        if np.amax(var3d[i,:,:])==varmax :
            zmax=i
#     print(zmax)
    return zmax

file = 'dynf'
for var in vars_[file]:
    dataL = np.squeeze(ds_dict[exptlist[0]][file][fhr][var][::])
    dataR = np.squeeze(ds_dict[exptlist[1]][file][fhr][var][::])
    if var == 'spfh':
        dataL = dataL * 1000.
        dataR = dataR * 1000.
    title = f'{var}: {fhr} hr fcst from {initdate} '
    #print(title)
    #plot_vertical_cross(dataL, dataR, start_lat, start_lon, end_lat, end_lon, numpoints, title, exptlist)

    # for 2 model runs with same verticall levels, use zmax from one of the two
    if diff_vert_lev:
        zmaxL=find_lev_max(dataL)
        zmaxR=find_lev_max(dataR)
        zmax=[zmaxL, zmaxR]
        data2dL=dataL[zmaxL,:,:]
        data2dR=dataR[zmaxR,:,:]
        title = f'{var}: {fhr} hr fcst from {initdate} '
        plot_data_zoom(data2dL, data2dR, start_lat, start_lon, end_lat, end_lon, title, exptlist, zmax)
    else:
        zmax=find_lev_max(dataL)
        data2dL=dataL[zmax,:,:]
        data2dR=dataR[zmax,:,:]
        title = f'{var}: {fhr} hr fcst from {initdate} at lev {zmax}'
        print(title)
        plot_data_zoom(data2dL, data2dR, start_lat, start_lon, end_lat, end_lon, title, exptlist)


# y1,x1=latlon_to_xy(start_lat, start_lon)
# y2,x2=latlon_to_xy(end_lat, end_lon)
# print(y1,x1)
# print(y2,x2)
# print(lat[y1,x1], lon[y1,x1])
# print(lat[y2,x2], lon[y2,x2])
