#!/usr/bin/env python3
import os, sys
import struct
import numpy as np
import ncepbufr
from netCDF4 import Dataset
from datetime import datetime, timedelta

from invdisttree import Invdisttree as INVTree

########################################################################

# Load each file in the files dict into a NetCDF Dataset

def load_Dataset(files):
    if isinstance(files, dict):
        assert(isinstance(files, dict))
        ds = dict()
        for k, v in files.items():
            if isinstance(v, dict):
                ds[k] = load_Dataset(v)
            else:
                ds[k] = Dataset(v, 'r')

    elif isinstance(files, list):
        ds = []
        for v in files:
            ds.append(Dataset(v,'r'))
    else:
        print("Unsupported argument type in load_Dataset")
        raise ValueError

    return ds

########################################################################

def read_bufrsfc(fileaname):
    ''' Read specific SFC observations from one PrepBufr file

        filter out missing values
    '''

    exp_reports = ['ADPSFC',]
    exp_mnemonic = ['SID','XOB', 'YOB', 'DHR']
    exp_tdmnem = ['TDO']
    exp_tcmnem = ['TOB']
    exp_uvmnem = ['UOB', 'VOB']

    tcobs = []
    tdobs = []
    uvobs = []

    #outfile  = 'rap.2019052018.adpsfc.txt'
    #outtable = 'rap.table'

    bufr = ncepbufr.open(filename)
    #bufr.advance()
    #bufr.load_subset()
    #bufr.dump_subset(outfile)
    #temp = bufr.read_subset('SSTH').squeeze()-273.15  # convert from Kelvin to Celsius

    #bufr.dump_table(outtable)
    #bufr.print_table()
    #msgs=bufr.inventory()
    while bufr.advance() == 0:
        #print(bufr.msg_counter,bufr.msg_type,bufr.msg_date, bufr.receipt_time,bufr.subsets)
        if bufr.msg_type in exp_reports:
            while bufr.load_subset() == 0:
                tco = bufr.read_subset('TOB')
                tdo = bufr.read_subset('TDO')
                #tdo = bufr.read_subset('QOB')
                uvo = bufr.read_subset('UOB VOB')
                hdr = bufr.read_subset('SID XOB YOB DHR').squeeze()
                if not np.ma.is_masked(tco):
                    tcobs.append((struct.pack('d',hdr[0]),hdr[1],hdr[2],hdr[3],tco[0][0]))
                if not np.ma.is_masked(tdo):
                    tdobs.append((struct.pack('d',hdr[0]),hdr[1],hdr[2],hdr[3],tdo[0][0]))
                if not np.ma.is_masked(uvo):
                    uvobs.append((struct.pack('d',hdr[0]),hdr[1],hdr[2],hdr[3],uvo[0][0],uvo[1][0]))

                #print(f"hdr={hdr},tc={tco},td={tdo},uob={uob},vob={vob}")

    bufr.close()

    obs = {'TOB': tcobs, 'TDO': tdobs, 'UVO': uvobs}

    return obs

########################################################################

def calc_td(t, p, qv):
    """Calculates dewpoint temp for a np.ndarray

    Adapted from similar wrftools routine (https://github.com/keltonhalbert/wrftools/)
    Copyright (c) 2015, Kelton Halbert

    Input:
        t - np.ndarray of temperature (K)
        p - np.ndarray of pressure (Pa)
        qv - np.ndarray of water vapor mixing ratio (Kg/Kg)

    Returns:

        td - np.ndarray of dewpoint temperature (K)
    """

    ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #class Constants:
    #
    #    def __init__(self):
    #        self.e_0 = 6.1173        # std atm vapor pressure (hPa)
    #        self.t_0 = 273.16        # std atm surface temp (K)
    #        self.Rv  = 461.50        # gas constant for water vapor (J/K Kg)
    #        self.Lv  = 2501000.      # latent heat of vaporization (J/Kg)
    #
    ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #
    #def vapor_pressure(t, p, qv, const):
    #    """ Returns vapor pressure """
    #    e_s = const.e_0 * np.exp((const.Lv / const.Rv) * ((1. / const.t_0) - (1. / t)))
    #                  # sat vapor pressure via Clasius Clapeyron equation (hPa)
    #    w_s = (0.622 * e_s) / (p - e_s)              # sat mixing ratio (g/Kg)
    #    rh = (qv / w_s) * 100.                       # relative humidity
    #    e = (rh / 100.) * e_s                        # vapor pressure (hPa)
    #
    #    return e
    #
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #
    #const = Constants()
    #
    #pc = p/100.         # convert to hPa
    #qvc = np.where(qv < 0.0, np.finfo(float).eps, qv*1000.)      # convert to g/Kg
    #
    #with np.errstate(all='raise'):
    #    try:
    #        e = vapor_pressure(t,pc,qvc,const)
    #        td = 1. / ((1. / const.t_0) - ((const.Rv / const.Lv) * np.log(e / const.e_0)))
    #    except FloatingPointError:
    #        x,y = np.nonzero(e <= 0.0)
    #        noe = e[x,y]
    #        print(f"{x[0]}, {y[0]} => {noe[0]}")
    #        print(f"{t[x[0],y[0]]},{p[x[0],y[0]]},{qv[x[0],y[0]]}")
    #        sys.exit(0)
    #
    #td[np.isnan(td)] = 0.
    #
    #return td

    """ Calculate the dewpoint temperature in Kelvin given the
    Temperature (TEMP), Pressure (PRES), and Mixing Ratio (QVAPOR).
    Arrays must be the same shape.
    ---------------------
    TEMP (numpy.ndarray): ndarray of 'normal' temperature in degrees Kelvin
    ---------------------
    returns:
        numpy.ndarray of dewpoint values (c) same shape as TEMP, PRES, and QVAPOR
    """

    pres   = p/100.         # convert to hPa
    qvapor = np.where(qv < 0.0, np.finfo(float).eps, qv)      # filter unphysical values
    t_inv  = 1 / t

    ## constants for the Clausius Clapeyron Equation
    E_0    = 6.1173           ## mb
    T_0    = 273.16           ## K
    RV     = 461.50           ## J K-1 Kg-1
    LV_0   = 2.501 * 10**6    ## J Kg-1

    K1     = LV_0 / RV
    K2     = 1 / T_0
    K1_INV = RV / LV_0

    ## Compute portions of the equation, Clausius-Clapeyron Equation
    e_s   = E_0 * np.exp( K1 * ( K2 - t_inv ) )     ## mb
    ## get saturation mixing ratio for RH
    w_s = ( 0.622 * e_s ) / ( pres - e_s )
    rh  = ( qvapor / w_s ) * 100
    ## back out the vapor pressure
    e   = ( rh / 100 ) * e_s ## mb
    ## compute individual terms when solving the equation for Td
    k4  = np.log( e / E_0 )
    td  = 1./( K2 - (K1_INV * k4) )

    return td

########################################################################

def process_obs(obs,lat2d,lon2d):
    ''' process observation with units and domain check'''

    latmax = lat2d.max()
    latmin = lat2d.min()
    lonmax = lon2d.max()
    lonmin = lon2d.min()
    #print(f"{latmin} - {latmax}; {lonmin} -> {lonmax}")

    validobs = {'tmp2m'  : {'sid': [], 'lat': [], 'lon': [], 'time': [], 'OBS': []},
                #'spfh2m' : {'sid': [], 'lat': [], 'lon': [], 'time': [], 'OBS': []},
                'td2m'   : {'sid': [], 'lat': [], 'lon': [], 'time': [], 'OBS': []},
                'ugrd10m': {'sid': [], 'lat': [], 'lon': [], 'time': [], 'OBS': []},
                'vgrd10m': {'sid': [], 'lat': [], 'lon': [], 'time': [], 'OBS': []}
               }

    for key,val in obs.items():
        for record in val:
            lat = record[2]
            lon = record[1]
            #print(f"lat={lat},lon={lon}")
            if lat > latmin and lat < latmax and lon > lonmin and lon < lonmax:
                if key == 'TOB':
                    validobs['tmp2m']['OBS'].append(record[4] + 273.16 )   # covert to Kelvin
                    validobs['tmp2m']['lat'].append(lat)
                    validobs['tmp2m']['lon'].append(lon)
                    validobs['tmp2m']['sid'].append(record[0])
                    validobs['tmp2m']['time'].append(record[3])
                elif key == 'TDO':
                    #validobs['spfh2m']['OBS'].append(record[4]/1000000)   # mg/kg to kg/kg
                    validobs['td2m']['OBS'].append(record[4] + 273.16)     # convert to Kelvin
                    validobs['td2m']['lat'].append(lat)
                    validobs['td2m']['lon'].append(lon)
                    validobs['td2m']['sid'].append(record[0])
                    validobs['td2m']['time'].append(record[3])
                else:
                    validobs['ugrd10m']['lat'].append(lat)
                    validobs['ugrd10m']['lon'].append(lon)
                    validobs['ugrd10m']['sid'].append(record[0])
                    validobs['ugrd10m']['time'].append(record[3])
                    validobs['ugrd10m']['OBS'].append(record[4])
                    validobs['vgrd10m']['lat'].append(lat)
                    validobs['vgrd10m']['lon'].append(lon)
                    validobs['vgrd10m']['sid'].append(record[0])
                    validobs['vgrd10m']['time'].append(record[3])
                    validobs['vgrd10m']['OBS'].append(record[5])

    return validobs

########################################################################

def get_filename(case,initdate,mem,ftime,filebase):
    '''
    case name in ['test_runs','test_spp','test_mp']
    initdate  in YYYYMMDDHH
    mem       in mem_xxx or a list of 'mem_xxx'
    ftime     in minutes
    filebase  in ['phyf','dynf']
    '''
    rootdir  = '/scratch/ywang/EPIC'

    itime = f'{ftime//60:03d}:{ftime%60:02d}:00'

    if isinstance(mem,str):
        filename = os.path.join(rootdir,case,initdate,mem,f'{filebase}{itime}.nc')

        if not os.path.lexists(filename):
            print(f"{filename} not exist.")
            raise ValueError

    elif isinstance(mem,list):
        filename = [os.path.join(rootdir,case,initdate,memid,f'{filebase}{itime}.nc') for memid in mem]

        for filen in filename:
            if not os.path.lexists(filen):
                print(f"{filen} not exist.")
                raise ValueError

    return filename


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if __name__ == "__main__":

    initdates=['2019050718','2019051618','2019052018']

    imin = 60
    fmin = 360
    ntimes = fmin//imin+1
    nens   = 40

    # Provide a file path to a forecast directory.
    # The example below creates a dictionary containing 2 experiments, expt, through the first 4 forecast hours (including 0)
    casenames = {'SPH': 'test_runs','SPP': 'test_spp','MPH': 'test_mp', 'MSP' : 'test_mspp'}
    #casenames = {'MPH': 'test_mp'}
    memlist=[f'mem_{mid:03d}' for mid in range(1,nens+1)]

    #print('available forecast files =', fmin//imin+1)

    #folders = { case: os.path.join(f'/scratch/ywang/EPIC/{casenames[case]}',f'{initdate}')
    #            for case in casenames
    #          }
    #filenms = [ f'{i//60:03d}:{i%60:02d}:00' for i in range(0,fmin+imin,imin)]
    #          #{  #'SPP': [ f'{i//60:03d}' for i in range(0,fmin+imin,imin)] }
    #          #  'SPP': [ f'{i//60:03d}:{i%60:02d}:00' for i in range(0,fmin+imin,imin)]
    #          #  'MPH': [ f'{i//60:03d}:{i%60:02d}:00' for i in range(0,fmin+imin,imin)]
    #          #}
    #
    #files = [{case:
    #            {mem:
    #                {
    #                    x: os.path.join(folders[case], mem, f'{x}{filenms[i]}.nc')
    #                    for x in ['phyf',]      #['dynf', 'phyf']
    #                } for mem in memlist
    #            } for case in casenames
    #        } for i in range(ntimes) ]

    #print(files['SPP'][memlist[0]]['dynf'])
    #print(files['SPH'][memlist[0]]['phyf'])
    #sys.exit(0)

    # Variables to plot from dynf and phyf files
    varens = {
         #'dynf': ['ugrd', 'vgrd', 'dzdt', 'tmp', 'spfh', 'dpres','delz'],
         'phyf': ['tmp2m', 'td2m', 'ugrd10m', 'vgrd10m' ],
         #'phyf' : ['td2m']
    }

    #
    # Compute time series of ensemble properties
    #
    fileindx = 'phyf'

    nensvar = len(varens[fileindx])
    ncases  = len(casenames.keys())

    #-------------------------------------------------------------------
    #
    # Read all observations for date and time
    #
    #-------------------------------------------------------------------
    validobs = [[None for y in range(0,fmin+imin,imin)] for x in range(0,len(initdates))]
    for d,initdate in enumerate(initdates):
        #
        # Read FV3 grid
        #
        file0 = get_filename(casenames['MPH'],initdate,memlist[0],0,'phyf')
        fh = Dataset(file0,'r')
        lat2d = np.degrees(fh['grid_yt'])
        lon2d = np.degrees(fh['grid_xt'])
        fh.close()

        t = 0
        for ftime in range(0,fmin+imin,imin):
            inittime = datetime.strptime(initdate,'%Y%m%d%H')
            currtime = inittime+timedelta(minutes=ftime)
            #
            # Read OBS and process obs
            #
            filename = f'/scratch/wof/realtime/OBSGEN/SBUFR/{currtime:%Y/%m/%d/rap.%Y%m%d%H}.prepbufr.tm00'
            sfcobs = read_bufrsfc(filename)
            validobs[d][t] = process_obs(sfcobs,lat2d,lon2d)
            t = t+1

    #-------------------------------------------------------------------
    #
    # Accumulate error variables
    #
    #-------------------------------------------------------------------

    # invdist tree interpolation parameters
    Nnear    = 4        # 8 2d, 11 3d => 5 % chance one-sided -- Wendel, mathoverflow.com
    leafsize = 10
    eps      = .1       # approximate nearest, dist <= (1 + eps) * true nearest
    p        = 1        # weights ~ 1 / distance**p

    var_rmse = np.zeros((nensvar,ncases,ntimes))
    var_bias = np.zeros((nensvar,ncases,ntimes))
    var_stds = np.zeros((nensvar,ncases,ntimes))
    var_cunt = np.zeros((nensvar,ncases,ntimes))
    var_nobs = np.zeros((nensvar,ncases,ntimes))

    for n,case in enumerate(casenames):

        for d,initdate in enumerate(initdates):

            t = 0
            for ftime in range(0,fmin+imin,imin):

                filenames = get_filename(casenames[case],initdate,memlist,ftime,fileindx)
                ds_files  = load_Dataset(filenames)

                fh = ds_files[0]
                lat2d = np.degrees(fh['grid_yt'])
                lon2d = np.degrees(fh['grid_xt'])
                nlat  = fh.ny
                nlon  = fh.nx

                gridloc = np.vstack((np.ravel(lon2d), np.ravel(lat2d))).T

                for v, var in enumerate(varens[fileindx]):

                    obsloc = np.vstack((validobs[d][t][var]['lon'],validobs[d][t][var]['lat'])).T
                    nobs   = len(validobs[d][t][var]['OBS'])

                    print(f"{case}:{initdate} for {var} at {ftime} min - nobs = {nobs}, nx*ny = {nlat*nlon}" )

                    vardiff = 0.0
                    varbias = 0.0
                    varcnt  = 0
                    #var2d = np.zeros((nens,nlat*nlon))
                    var2d = np.zeros((nens,nobs))
                    for m,mem in enumerate(memlist):
                        print(f"Computing {case}:{mem} for {var} at {ftime} min" )
                        fh = ds_files[m]
                        if var == 'td2m':
                            sph2m = np.squeeze(fh['spfh2m'])
                            tmp2m = np.squeeze(fh['tmp2m'])
                            presfc = np.squeeze(fh['pressfc'])
                            vardata = calc_td(tmp2m, presfc, sph2m)
                            #print(vardata[0:2,0:2])
                            #print(vardata.max(),vardata.min(),vardata.mean())
                            #varobs = validobs[var]['OBS']
                            #print(np.max(varobs),np.min(varobs),np.mean(varobs))
                            #sys.exit(0)
                        else:
                            vardata = np.squeeze(fh[var])   # variable on model grid
                        invdisttree   = INVTree( gridloc, np.ravel(vardata), leafsize=leafsize, stat=1 )
                        varobs = invdisttree( obsloc, nnear=Nnear, eps=eps, p=p )  # dimensioned by number of observations
                        #plot_regular(vardata,f"{var}_grid",i,fh)
                        #plot_irregular(varobs,validobs[var]['lat'],validobs[var]['lon'],f"{var}_obs",i,fh)
                        vardiff += np.sum( (varobs-validobs[d][t][var]['OBS'])**2 )
                        varbias += np.sum( (varobs-validobs[d][t][var]['OBS'])    )
                        #var2d[m,:] = vardata.ravel()
                        var2d[m,:]  = varobs
                        varcnt  += len(varobs)
                        assert(nobs == len(varobs))

                    #varstd = var2d.std(axis=0)         # std over each grid point over ##ensemble members
                    #var_stds[v,n,t] += varstd.mean()    # std over each grid point over observation points

                    varstd = var2d.std(axis=0)         # std over each grid point over observation points
                    var_stds[v,n,t] += np.sum(varstd)
                    var_rmse[v,n,t] += vardiff
                    var_bias[v,n,t] += varbias
                    var_cunt[v,n,t] += varcnt
                    var_nobs[v,n,t] += nobs

                t = t+1

    #-------------------------------------------------------------------
    #
    # Compute RSMS, STD, Bias etc.
    #
    #-------------------------------------------------------------------
    for v in range(nensvar):
        for n in range(ncases):
            for t in range(ntimes):
                #var_stds[v,n,t] = var_stds[v,n,t]/len(initdates)
                assert(var_nobs[v,n,t]*len(memlist) == var_cunt[v,n,t])
                var_stds[v,n,t] = var_stds[v,n,t] / var_nobs[v,n,t]
                var_rmse[v,n,t] = np.sqrt( var_rmse[v,n,t] / var_cunt[v,n,t] )
                var_bias[v,n,t] = var_bias[v,n,t] / var_cunt[v,n,t]

    #-------------------------------------------------------------------
    #
    # Dump RSMS, STD, Bias etc.
    #
    #-------------------------------------------------------------------
    outfile = f"/oldscratch/ywang/EPIC/Program/plots/stat2m_norm.aggregation.npy"
    with open(outfile,'wb') as f:
        np.save(f,var_rmse)
        np.save(f,var_bias)
        np.save(f,var_stds)

