#!/usr/bin/env python

import cape as capecalc
import numpy as np
import os

__all__ = ['mucape','mlcape1','mlcape2','sfcape']

pwd = os.getcwd()
lookup_file = os.path.join(pwd,'psadilookup.dat')

def mucape(prs_mb,tmp,mixr,hgt,ter,psfc_mb,ter_follow,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension, but code will check for this **

    output
    ------
    mucape: 2D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 2D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 2D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    cape_type = 0
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    mucape,mucin,mulcl,mulfc = capecalc.dcapecalc3d(prs_mb[:,:,:],tmp[:,:,:],mixr[:,:,:],hgt[:,:,:],ter[:,:],psfc_mb[:,:],cape_type,ter_follow,lookup_file)

    return mucape,mucin,mulcl,mulfc

def mlcape1(prs_mb,tmp,mixr,hgt,ter,psfc_mb,ter_follow,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *Lifts mixed-layer parcel quantity from surface

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    cape_type = 2
    mlcape,mlcin,mllcl,mllfc = capecalc.dcapecalc3d(prs_mb[:,:,:],tmp[:,:,:],mixr[:,:,:],hgt[:,:,:],ter[:,:],psfc_mb[:,:],cape_type,ter_follow,lookup_file)
    return mlcape,mlcin,mllcl,mllfc


def mlcape2(prs_mb,tmp,mixr,hgt,ter,psfc_mb,ter_follow,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *UPP version of MLCAPE (lifts from mixed-layer level)

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    cape_type = 3
    mlcape,mlcin,mllcl,mllfc = capecalc.dcapecalc3d(prs_mb[:,:,:],tmp[:,:,:],mixr[:,:,:],hgt[:,:,:],ter[:,:],psfc_mb[:,:],cape_type,ter_follow,lookup_file)
    return mlcape,mlcin,mllcl,mllfc

def sfcape(prs_mb,tmp,mixr,hgt,ter,psfc_mb,ter_follow,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    sfcape: 2D array of surface based cape (J/kg)
    sfcin: 2D array of surface based cin (J/kg)
    sflcl: 2D array of surface based lcl height (m)
    sflfc: 2D array of surface based lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    cape_type = 1
    sfcape,sfcin,sflcl,sflfc = capecalc.dcapecalc3d(prs_mb[:,:,:],tmp[:,:,:],mixr[:,:,:],hgt[:,:,:],ter[:,:],psfc_mb[:,:],cape_type,ter_follow,lookup_file)
    return sfcape,sfcin,sflcl,sflfc

