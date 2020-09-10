#!/usr/bin/env python3
import os, sys
import struct
from datetime import datetime, timedelta
import numpy as np
import ncepbufr

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

def read_snd(filename,outtable=None):

    bufr = ncepbufr.open(filename)
    #bufr.advance()
    #bufr.load_subset()
    #bufr.dump_subset(outfile)
    #temp = bufr.read_subset('SSTH').squeeze()-273.15  # convert from Kelvin to Celsius

    if outtable is not None: bufr.dump_table(outtable)
    #bufr.print_table()
    #msgs=bufr.inventory()

    stations = {}
    while bufr.advance() == 0:
        station = {}
        if bufr.msg_type == "ADPUPA":
            #print(bufr.msg_counter,bufr.msg_type,bufr.msg_date, bufr.receipt_time,bufr.subsets)
            msgtime = datetime.strptime(str(bufr.msg_date),'%Y%m%d%H')

            while bufr.load_subset() == 0:
                sco = bufr.read_subset('SID').squeeze()
                typ = bufr.read_subset('TYP').squeeze()
                xob = bufr.read_subset('XOB').squeeze()
                yob = bufr.read_subset('YOB').squeeze()
                dhr = bufr.read_subset('DHR').squeeze()

                tob = bufr.read_subset('TOB').squeeze()
                tdo = bufr.read_subset('TDO').squeeze()
                pob = bufr.read_subset('POB').squeeze()
                uob = bufr.read_subset('UOB').squeeze()
                vob = bufr.read_subset('VOB').squeeze()
                #ddo = bufr.read_subset('DDO').squeeze()
                #sso = bufr.read_subset('SOB').squeeze()

                sid = struct.pack('d',sco).decode("utf-8")
                obt = msgtime + timedelta(hours=dhr.compressed()[0])
                print(sid,typ,bufr.msg_date, dhr.compressed())
                if sid not in stations.keys():
                    stations[sid] = {'lat': yob, 'lon': xob if xob < 180 else xob-360.,
                                     'time': obt}

                if typ < 200:               # 100-199 for mass observations
                    #stations[sid]['ps']  = pso.tolist(-9999.) * units.hPa
                    #stations[sid]['T']   = tco.tolist(-9999.) * units.degC
                    #stations[sid]['Td']  = tdo.tolist(-9999.) * units.degC
                    mask = np.any([tob.mask, tdo.mask], axis=0)
                    pob.mask = mask; tob.mask = mask; tdo.mask = mask
                    pso = pob.compressed()       # remove missing tob or tdo levels
                    tds = tdo.compressed()
                    tso = tob.compressed()

                    stations[sid]['ps']  = pso * units.hPa
                    stations[sid]['T']   = tso * units.degC
                    stations[sid]['Td']  = tds * units.degC
                else:                       # 200-299 for wind observations
                    mask = np.any([uob.mask, vob.mask], axis=0)
                    pob.mask = mask; uob.mask = mask; vob.mask = mask
                    pso = pob.compressed()        # remove missing uob or vob levels
                    uso = uob.compressed()
                    vso = vob.compressed()

                    stations[sid]['pv']  = pso * units.hPa
                    stations[sid]['u']   = uso * units.meter/units.second
                    stations[sid]['v']   = vso * units.meter/units.second
                    #stations[sid]['dir'] = ddo.tolist(-9999.) * units.degree
                    #stations[sid]['spd'] = sso.tolist(-9999.) * units.meter/units.second

                #    stations.append(station)
                #if np.ma.is_masked(tco[0]): continue
                #else: break
                #print(tco)
                #print(tdo)
                #print(pso)
                #print(uvo)
                #sys.exit(0)
                #bufr.dump_subset(outfile,append=True)
    bufr.close()

    return stations


def plot_skewt(p,t,td,puv=None,u=None,v=None,title=None,outfile=None):


    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(9, 9))
    skew = SkewT(fig, rotation=30)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p, t, 'r', linewidth=2)
    skew.plot(p, td, 'g', linewidth=2)
    if u is not None and v is not None:
        skew.plot_barbs(puv, u, v)

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 60)

    # Calculate the LCL
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], t[0], td[0])

    # Calculate the parcel profile.
    parcel_prof = mpcalc.parcel_profile(p, t[0], td[0]).to('degC')

    # Plot LCL temperature as black dot
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

    # Plot the parcel profile as a black line
    skew.plot(p, parcel_prof, 'k--', linewidth=1)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p, t, parcel_prof)
    skew.shade_cape(p,t, parcel_prof)

    # Plot a zero degree isotherm
    #skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    if title is not None:
        plt.title(title)

    # Show the plot
    #plt.show()
    if outfile is None: outfile = 'skewt.png'
    fig.savefig(outfile, format='png')

if __name__ == "__main__":
    filename = '/scratch/wof/realtime/OBSGEN/SBUFR/2019/05/20/rap.2019052018.prepbufr.tm00'
    outtable = 'rap.table'

    stations = read_snd(filename)
    for sid,stn in stations.items():
        if 'pv' not in stn:
            stn['pv'] = None
            stn['u']  = None
            stn['v']  = None
            lenuv = 0
        else:
            lenuv = len(stn['pv'])

        lenmass = len(stn['ps'])

        if len(stn['ps']) > 10:
            outfile = f'2019052018_{sid.rstrip()}.png'
            title = f"{sid} {stn['time']:%Y%m%d %H:%M:%S} ({stn['lat']:.2f}\N{DEGREE SIGN}N, {stn['lon']:.2f}\N{DEGREE SIGN}E)"
            print(f"Ploting  {sid} {lenmass},{lenuv} to {outfile}...")
            plot_skewt(stn['ps'],stn['T'],stn['Td'],stn['pv'],stn['u'],stn['v'],title,outfile)
        else:
            print(f"Skipping {sid} {lenmass},{lenuv}")
