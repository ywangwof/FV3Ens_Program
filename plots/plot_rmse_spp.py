#!/usr/bin/env python3
import os, sys
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

########################################################################

def plot_series(name,var1ds,times,field,initdate, cases):

    '''
    Input parameters:

        name:   Plot name
        var1ds: list of variables to be plot, 2D (n,t)
        time:   List of times
        field:  Filed name
        initdate:
        cases:  List of cases, 1D (n,)

    Plot time series of RMSE or BIAS

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

########################################################################

def plot_all(var_stds,var_rmse,var_bias,times,varens,initdate, cases):

    '''
    Input parameters:

        name:   Plot name
        var1ds: list of variables to be plot, 2D (n,t)
        time:   List of times
        field:  Filed name
        initdate:
        cases:  List of cases, 1D (n,)

    Plot time series of RMSE or BIAS

    '''

    nrows = 2
    ncols = 2

    fig, ax = plt.subplots(nrows, ncols, figsize=(11, 8))
    axs = ax.flat
    #print(ax.shape,axs.shape)

    colors=['red','green','blue','cyan']
    unitss=['K','K','m/s','m/s']

    for n,var in enumerate(varens):
        std2d = var_stds[n]
        rms2d = var_rmse[n]
        bia2d = var_bias[n]
        if var == 'spfh2m':
            std2d = 1000.*std2d
            rms2d = 1000.*rms2d
            bia2d = 1000.*bia2d

        print(f"{var}: RMS ---")
        for var1d, c, case in zip(rms2d,colors,cases):
            print(case,var1d)
            axs[n].plot(times,var1d,'.--',color=c,label=f"{case} rmse")
        print(f"{var}: STD ---")
        for var1d,c,case in zip(std2d,colors,cases):
            print(case,var1d)
            axs[n].plot(times,var1d,'o-',color=c,label=f"{case} std")
        print(f"{var}: BIAS ---")
        for var1d, c,case in zip(bia2d,colors,cases):
            print(case,var1d)
            axs[n].plot(times,var1d,'*:',color=c,label=f"{case} bias")
        print("=== === ===")

        units = unitss[n]
        if n > 1: axs[n].set(xlabel="forecast hours",ylabel=f"{var} ({units})")
        else:     axs[n].set(                        ylabel=f"{var} ({units})")
        #axs[n].label_outer()
        #ax[n].legend(cases)

        title = f'{var}: domain average at {initdate}'
        axs[n].set_title(title,fontsize=8)
        handles, labels = axs[n].get_legend_handles_labels()

    #fig.subplots_adjust(bottom=0.90)
    #cbar_ax = fig.add_axes([0.92, 0.92, 0.2, 0.8])
    #fig.legend(handles,labels,loc="best")
    fig.subplots_adjust(right=0.85)
    #cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.5])
    fig.legend(handles,labels, loc="center right")

    #plt.show()
    figname = f"stat2m_{initdate}_spp.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')

########################################################################

def plot_indiv(varname,var_in,times,varens,initdate, cases):

    '''
    Input parameters:

        varname:   Plot name, rmse, std, bias
        var_in:    variables to be plot, 3D (n,c,t)
        times:     List of times, (t)
        varens:    list of variables, (n)
        initdate:
        cases:     List of cases, 1D (c,)

    Plot time series of RMSE, std or BIAS

    '''

    nrows = 2
    ncols = 2

    fig, ax = plt.subplots(nrows, ncols, figsize=(11, 8))
    axs = ax.flat
    #print(ax.shape,axs.shape)

    colors=['red','green','blue','cyan']
    unitss=['K','K','m/s','m/s']

    for n,var in enumerate(varens):
        var2d = var_in[n]
        if var == 'spfh2m':
            var2d = 1000.*var2d

        print(f"{var}: {varname} ------------------------------------")
        for var1d, c, case in zip(var2d,colors,cases):
            print(case,var1d)
            axs[n].plot(times,var1d,'.--',color=c,label=f"{case} {varname}")
        print("=======================================================")

        units = unitss[n]
        if n > 1: axs[n].set(xlabel="forecast hours",ylabel=f"{var} ({units})")
        else:     axs[n].set(                        ylabel=f"{var} ({units})")
        #axs[n].label_outer()
        #ax[n].legend(cases)

        title = f'{varname}: domain average at {initdate}'
        axs[n].set_title(title,fontsize=8)
        handles, labels = axs[n].get_legend_handles_labels()

    #fig.subplots_adjust(bottom=0.90)
    #cbar_ax = fig.add_axes([0.92, 0.92, 0.2, 0.8])
    #fig.legend(handles,labels,loc="best")
    fig.subplots_adjust(right=0.85)
    #cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.5])
    fig.legend(handles,labels, loc="center right")

    #plt.show()
    figname = f"stat2m_{varname}_{initdate}_spp.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')

########################################################################

def plot_irregular(var,lat,lon,field,fmin,fh):

    '''
    Input parameters:

        fhs:   list of file handles
        field: string of field to be plotted
        expt:  list of description of each files
        fmin:  Valid forecast in minutes

    Draws a Basemap representation with the contoured data overlayed,
    with a colorbar for each experiment, and the difference between the two.

    '''
    import matplotlib.tri as tri

    xsize = (fh.nx-1)*fh.dx
    ysize = (fh.ny-1)*fh.dy

    xlat = np.degrees(fh['grid_yt'][::])
    xlon = np.degrees(fh['grid_xt'][::])


    title = f'{field}: {fmin//60} hr fcst from {initdate}'

    fig, axs = plt.subplots(1, 1, figsize=(11, 8))
    #axs = ax.flatten()
    #print(ax.shape,axs.shape)

    m = Basemap(width=xsize,height=ysize,
                   rsphere=(6378137.00,6356752.3142),\
                   resolution='l',area_thresh=1000.,projection='lcc',\
                   lat_1=fh.stdlat1,lat_2=fh.stdlat2,lat_0=fh.cen_lat,lon_0=fh.cen_lon,
                   suppress_ticks=True,  ax=axs )

    m.drawstates(linewidth=0.1)
    m.drawcoastlines(linewidth=0.1)
    m.drawmapboundary(linewidth=0.1)
    m.drawcountries(linewidth=0.1)
    #m.drawparallels(np.arange(-90.,120.,5),labels=[1,0,0,0]);
    #m.drawmeridians(np.arange(-180.,180.,10),labels=[0,0,0,1]);
    x, y = m(lon, lat)

    xgrid, ygrid = m(xlon,xlat)

    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, var)
    data2d = interpolator(xgrid, ygrid)

    #contours, cm = eq_contours(var)
    contours = np.arange(250.0,350.,10.)

    # Draw the contoured data over the map
    cs = m.contourf(xgrid, ygrid, data2d, contours, ax=axs)
    #m.plot(xx, yy, color='black', linewidth=3, linestyle='--', ax=ax[i])
    #fig.colorbar(cs, ax=axs[n], orientation='vertical', shrink=0.7);
    axs.set_title(field,fontsize=8)
    handles, labels = axs.get_legend_handles_labels()

    plt.suptitle(f"{title}", fontsize=10)

    fig.subplots_adjust(right=0.90)
    cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.5])
    fig.colorbar(cs, cax=cbar_ax, orientation='vertical', shrink=0.7)

    #plt.show()
    figname = f"{field}_{fmin//60:02d}h.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')

########################################################################

def plot_regular(var,field,fmin,fh):

    '''
    Input parameters:

        fhs:   list of file handles
        field: string of field to be plotted
        expt:  list of description of each files
        fmin:  Valid forecast in minutes

    Draws a Basemap representation with the contoured data overlayed,
    with a colorbar for each experiment, and the difference between the two.

    '''
    import matplotlib.tri as tri

    xsize = (fh.nx-1)*fh.dx
    ysize = (fh.ny-1)*fh.dy

    xlat = np.degrees(fh['grid_yt'][::])
    xlon = np.degrees(fh['grid_xt'][::])


    title = f'{field}: {fmin//60} hr fcst from {initdate}'

    fig, axs = plt.subplots(1, 1, figsize=(11, 8))
    #axs = ax.flatten()
    #print(ax.shape,axs.shape)

    m = Basemap(width=xsize,height=ysize,
                   rsphere=(6378137.00,6356752.3142),\
                   resolution='l',area_thresh=1000.,projection='lcc',\
                   lat_1=fh.stdlat1,lat_2=fh.stdlat2,lat_0=fh.cen_lat,lon_0=fh.cen_lon,
                   suppress_ticks=True,  ax=axs )

    m.drawstates(linewidth=0.1)
    m.drawcoastlines(linewidth=0.1)
    m.drawmapboundary(linewidth=0.1)
    m.drawcountries(linewidth=0.1)
    #m.drawparallels(np.arange(-90.,120.,5),labels=[1,0,0,0]);
    #m.drawmeridians(np.arange(-180.,180.,10),labels=[0,0,0,1]);

    xgrid, ygrid = m(xlon,xlat)

    #contours, cm = eq_contours(var)
    contours = np.arange(250.0,350.,10.)

    # Draw the contoured data over the map
    cs = m.contourf(xgrid, ygrid, var, contours, ax=axs)
    #m.plot(xx, yy, color='black', linewidth=3, linestyle='--', ax=ax[i])
    #fig.colorbar(cs, ax=axs[n], orientation='vertical', shrink=0.7);
    axs.set_title(field,fontsize=8)
    handles, labels = axs.get_legend_handles_labels()

    plt.suptitle(f"{title}", fontsize=10)

    fig.subplots_adjust(right=0.90)
    cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.5])
    fig.colorbar(cs, cax=cbar_ax, orientation='vertical', shrink=0.7)

    #plt.show()
    figname = f"{field}_{fmin//60:02d}h.png"
    print(f"Saving figure to {figname} ...")
    fig.savefig(figname, format='png')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if __name__ == "__main__":

    initdate='2019052018'

    imin = 60
    fmin = 360
    ntimes = fmin//imin+1
    nens   = 40

    # Provide a file path to a forecast directory.
    # The example below creates a dictionary containing 2 experiments, expt, through the first 4 forecast hours (including 0)
    casenames = {'v0': 'test_spp','v1': 'test_spp','v4': 'test_spp','v5': 'test_spp'}

    #print('available forecast files =', fmin//imin+1)
    fileindx = 'phyf'

    # Variables to plot from dynf and phyf files
    varens = {
         #'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
         'phyf': ['tmp2m', 'td2m', 'ugrd10m', 'vgrd10m' ],
    }

    outfile = f"/oldscratch/ywang/EPIC/Program/plots/stat2m.{initdate}_spp_v0145.npy"
    #
    # Compute time series of ensemble properties
    #
    with open(outfile,'rb') as f:
        var_rmse = np.load(f)
        var_bias = np.load(f)
        var_stds = np.load(f)

    #for v, var in enumerate(varens[fileindx]):
    #    print(f"Plotting {var} ..." )
    #    plot_series('RMSE',var_rmse[v],range(fmin//60+1),var,initdate,casenames)
    #    plot_series('BIAS',var_bias[v],range(fmin//60+1),var,initdate,casenames)
    #

    #
    # plot rmse, std, bias in one figure
    #
    #plot_all(var_stds,var_rmse,var_bias,range(fmin//60+1),varens[fileindx],initdate,casenames)

    #
    # plot rmse, std, bias in its own figure individually
    #
    plot_indiv('rmse',var_rmse,range(fmin//60+1),varens[fileindx],initdate,casenames)
    plot_indiv('std_',var_stds,range(fmin//60+1),varens[fileindx],initdate,casenames)
    plot_indiv('bias',var_bias,range(fmin//60+1),varens[fileindx],initdate,casenames)
