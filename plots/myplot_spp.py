#!/usr/bin/env python3

import os, sys, glob, math
#from argparse import Namespace

import numpy as np
from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

from datetime import datetime, timedelta

########################################################################
#
# Load the dictionary into a Namespace data structure.
# This step is not necessary, but cuts down the syntax needed to reference each item in the dict.
#
# Example: Retrieve the 0 hr forecast Dataset from GFS Dynamics
#            dict: ds_dict['GFS']['dynf'][0]
#       Namespace: datasets.GFS.dynf[0]

def make_namespace(d: dict):
    assert(isinstance(d, dict))
    ns =  Namespace()
    for k, v in d.items():
        if isinstance(v, dict):
            leaf_ns = make_namespace(v)
            ns.__dict__[k] = leaf_ns
        else:
            ns.__dict__[k] = v

    return ns

########################################################################

# Load each file in the files dict into a NetCDF Dataset

def load_Dataset(files: dict):
    assert(isinstance(files, dict))
    ds = dict()
    for k, v in files.items():
        if isinstance(v, dict):
            ds[k] = load_Dataset(v)
        else:
            ds[k] = [Dataset(f, 'r') for f in list(v)]

    return ds

########################################################################

def eq_contours(indata, fieldname=""):

    '''
    Returns a balanced set of contours for data that has negative values.
    Also returns default colorbar to use for balanced, vs all positive values.
    '''

    if fieldname.startswith('refl'):
        c5 =  (0.0,                 0.9254901960784314, 0.9254901960784314)
        c10 = (0.00392156862745098, 0.6274509803921569, 0.9647058823529412)
        c15 = (0.0,                 0.0,                0.9647058823529412)
        c20 = (0.0,                 1.0,                0.0)
        c25 = (0.0,                 0.7843137254901961, 0.0)
        c30 = (0.0,                 0.5647058823529412, 0.0)
        c35 = (1.0,                 1.0,                0.0)
        c40 = (0.9058823529411765,  0.7529411764705882, 0.0)
        c45 = (1.0,                 0.5647058823529412, 0.0)
        c50 = (1.0,                 0.0,                0.0)
        c55 = (0.8392156862745098,  0.0,                0.0)
        c60 = (0.7529411764705882,  0.0,                0.0)
        c65 = (1.0,                 0.0,                1.0)
        c70 = (0.6,                 0.3333333333333333, 0.788235294117647)
        c75 = (0.0,                 0.0,                0.0)

        nws_dz_cmap = mpl.colors.ListedColormap([c20, c25, c30, c35, c40, c45, c50, c55, c60, c65, c70])

        return np.arange(20.0,80.,5.), nws_dz_cmap

    else:

      cmap=plt.cm.get_cmap('bwr',20)
      cmaplist = [cmap(i) for i in range(cmap.N)]
      cmaplist[9] = [1,1,1,1]
      cmaplist[10] = [1,1,1,1]
      cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

      indata2 = None
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
          if c[0] == 0:  c[0] = c[1]/10
          return c, 'jet'


########################################################################

def plot_files(fhs,field,fmin, initdate, case, expt):

    '''
    Input parameters:

        fhs:   list of file handles
        field: string of field to be plotted
        expt:  list of description of each files
        fmin:  Valid forecast in minutes

    Draws a Basemap representation with the contoured data overlayed,
    with a colorbar for each experiment, and the difference between the two.

    '''
    fh = fhs[0]

    xsize = (fh.nx-1)*fh.dx
    ysize = (fh.ny-1)*fh.dy

    lat = np.degrees(fh['grid_yt'][::])
    lon = np.degrees(fh['grid_xt'][::])

    title = f'{field}: {fmin//60} hr fcst from {initdate}'

    ntotal = len(fhs)
    ncols  = 8
    nrows  = ntotal//ncols

    fig, ax = plt.subplots(nrows, ncols, figsize=(11, 8.5), dpi=600)
    axs = ax.flat
    #print(ax.shape,axs.shape)

    for n, fh in enumerate(fhs):
        data3d = np.squeeze(fh[field][::])
        #data2d = np.amax(data3d,axis=0)    # for relfectivity
        data2d = data3d[0,:,:]    # bottom layer

        m = Basemap(width=xsize,height=ysize,
                   rsphere=(6378137.00,6356752.3142),\
                   resolution='l',area_thresh=1000.,projection='lcc',\
                   lat_1=fh.stdlat1,lat_2=fh.stdlat2,lat_0=fh.cen_lat,lon_0=fh.cen_lon,
                   suppress_ticks=True,  ax=axs[n] )

        m.drawstates(linewidth=0.1)
        m.drawcoastlines(linewidth=0.1)
        m.drawmapboundary(linewidth=0.1)
        m.drawcountries(linewidth=0.1)
        #m.drawparallels(np.arange(-90.,120.,5),labels=[1,0,0,0]);
        #m.drawmeridians(np.arange(-180.,180.,10),labels=[0,0,0,1]);
        x, y = m(lon, lat)

        #contours, cm = eq_contours(data2d,'refl')
        contours, cm = eq_contours(data2d,field)

        # Draw the contoured data over the map
        cs = m.contourf(x, y, data2d, contours, cmap=cm, ax=axs[n])
        #m.plot(xx, yy, color='black', linewidth=3, linestyle='--', ax=ax[i])
        #fig.colorbar(cs, ax=axs[n], orientation='vertical', shrink=0.7);
        axs[n].set_title(expt[n],fontsize=8)
        handles, labels = axs[n].get_legend_handles_labels()

    plt.suptitle(f"{title} {case}", fontsize=10)

    fig.subplots_adjust(right=0.90)
    cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.5])
    fig.colorbar(cs, cax=cbar_ax, orientation='vertical', shrink=0.7)

    #plt.show()
    figname = f"{field}_{fmin//60:02d}h_spp.{case}.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')
    plt.close()

########################################################################

def plot_series(name,var1ds,times,field,initdate, cases):

    '''
    Input parameters:

        fhs:   list of file handles
        field: string of field to be plotted
        expt:  list of description of each files
        fmin:  Valid forecast in minutes

    Draws a Basemap representation with the contoured data overlayed,
    with a colorbar for each experiment, and the difference between the two.

    '''

    title = f'{field}: domain average {name} from {initdate}'


    fig = plt.figure(figsize=(11, 8.5))
    for var1d in var1ds:
        plt.plot(times,var1d)
    plt.xlabel("forecast hours")
    plt.ylabel(f"{name} for {field}")
    plt.legend(cases)

    plt.title(title, fontsize=10)

    #plt.show()
    figname = f"{field}_{name}.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

if __name__ == "__main__":

    initdate='2019052018'

    imin = 60
    fmin = 360
    ntimes = fmin//imin+1
    nens   = 40

    # Provide a file path to a forecast directory.
    # The example below creates a dictionary containing 2 experiments, expt, through the first 4 forecast hours (including 0)
    #casenames = {'SPH': 'test_runs','SPP': 'test_spp', 'MPH': 'test_mp'}
    casenames = {'v0': 'test_spp','v1': 'test_spp','v4': 'test_spp','v5': 'test_spp'}
    memlist=[f'mem_{mid:03d}' for mid in range(1,nens+1)]

    #print('available forecast files =', fmin//imin+1)

    folders = { case: os.path.join(f'/scratch/ywang/EPIC/{casenames[case]}',f'{initdate}.{case}')
                for case in casenames
               }
    filenms = [ f'{i//60:03d}:{i%60:02d}:00' for i in range(0,fmin+imin,imin)]

    files = {case:
                {mem:
                    {
                        x: [os.path.join(folders[case], mem, f'{x}{filenms[i]}.nc') for i in range(ntimes) ]
                        for x in ['phyf',]      #['dynf', 'phyf']
                    } for mem in memlist
                } for case in casenames
            }

    #print(files['SPP'][memlist[0]]['dynf'])
    #print(files['SPH'][memlist[0]]['phyf'])
    #sys.exit(0)

    # Variables to plot from dynf and phyf files
    varplots = {
         #'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
         #'phyf': ['refl_10cm', ],
         'phyf': ['shum_wts', 'skebu_wts', 'skebv_wts'],

    }

    varens = {
         #'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
         'phyf': ['tmp2m', 'spfh2m', 'ugrd10m', 'vgrd10m' ],
    }

    #
    # Plot stamp plot of the ensemble members
    #
    fileindx = 'phyf'
    for case in casenames:
        fmins =[60,180,360]
        for fmin in fmins:
            print(f"Ploting {case} at {fmin} minutes forecast ...")
            filenames = []
            i = fmin//60
            for mem in memlist:
                filenames.append(files[case][mem][fileindx][i])

            fhs = [Dataset(f, 'r') for f in filenames]

            for var in varplots[fileindx]:
                plot_files(fhs, var, fmin, initdate, case, memlist)

    #
    # Compute time series of ensemble properties
    #
    #fileindx = 'phyf'
    #
    #varstds  = {var: [] for var in varens[fileindx]}
    #for case in casenames:
    #    ds_dict = load_Dataset(files[case])
    #
    #    fh = ds_dict[memlist[0]][fileindx][0]
    #    nlat = fh.ny
    #    nlon = fh.nx
    #    vardata = np.zeros((nens,ntimes,nlat,nlon))
    #    var1d   = np.zeros((nens,ntimes,nlat*nlon))
    #
    #    for v, var in enumerate(varens[fileindx]):
    #        for n,mem in enumerate(memlist):
    #            for t in range(ntimes):
    #                fh = ds_dict[mem][fileindx][t]
    #                vardata[n,t,:,:] = np.squeeze(fh[var])
    #                var1d[n,t,:] = vardata[n,t,:,:].ravel()
    #                #print(f"{mem}-{t}: {var} -> {fh[var].shape}")
    #        #varmean = var1d.mean(axis=0)
    #        varstd  = var1d.std(axis=0)
    #        varstd1d = varstd.mean(axis=1)
    #        varstds[var].append(varstd1d)
    #        #print(vardata.shape)
    #        #print(var1d.shape)
    #        #print(varmean.shape)
    #        #print(varstd.shape)
    #        #print(varstd1d.shape)
    #        #sys.exit(0)
    #
    #for var in varens[fileindx]:
    #    plot_series('STD',varstds[var],range(fmin//60+1),var,initdate,casenames)


