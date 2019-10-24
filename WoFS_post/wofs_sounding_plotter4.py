#!/usr/bin/env python
import glob
import os
import numpy as np
from netCDF4 import Dataset
import sharppy.sharptab.profile as profile
import datetime
from test_wof_sounding_plot import plot_wof
import sharppy.sharptab.utils as utils
import sharppy.sharptab.thermo as thermo
import argparse
import time
import multiprocessing
import sys
import warnings

warnings.filterwarnings("ignore",message="converting a masked element to nan",append=True)
warnings.filterwarnings("ignore",message="overflow encountered in multiply",append=True)
warnings.simplefilter(action="ignore",category=RuntimeWarning)
warnings.simplefilter(action="ignore",category=UserWarning)

############## Parse command line options #############

parser = argparse.ArgumentParser()
#parser.add_argument("indir", type=str, default= None, help="Input directory of member directories containing WRFOUT files to process")
parser.add_argument("summary_dir", type=str, default= None, help="Input directory for ENV summary files to process")
parser.add_argument("fcst_nt", type=int, help = "Number of timesteps to process")
parser.add_argument("outdir", type=str, help = "Output Directory (for summary file)")

args = parser.parse_args()
#indir = args.indir
summary_dir = args.summary_dir
outdir = args.outdir
fcst_nt = args.fcst_nt

############## Subroutine to call plot_wof subroutine in test_wof_sounding_plot.py #############

def make_plot(index):
    #index = (r,c)
    r,c = index
    print(r,c)
    xlat = xlats[r,c]
    xlon = xlons[r,c]

    mean_hgt = np.nanmean(z_agl[:,:,r,c], axis=0)
    mean_t = np.nanmean(t[:,:,r,c], axis=0)
    mean_q = np.nanmean(q[:,:,r,c], axis=0)
    mean_p = np.nanmean(p[:,:,r,c], axis=0)
    mean_td = thermo.temp_at_mixrat(mean_q*1000., mean_p)
    mean_u = np.nanmean(u[:,:,r,c], axis=0)
    mean_v = np.nanmean(v[:,:,r,c], axis=0)    
    mean_omega = np.nanmean(omega[:,:,r,c], axis=0)
    print("start prof")
    prof = profile.create_profile(profile='convective', pres=mean_p, hght=mean_hgt, tmpc=mean_t, \
                    dwpc=mean_td, u=mean_u, v=mean_v, omeg=mean_omega, missing=-9999, strictQC=False, latitude=xlat,date=vdateobj)

    print("end prof")
    member_profs = np.empty(len(z_agl[:,1]), dtype=object)
    print("start member profs")
    # Loop over all of the members and create profile objects for each of them.
    for m in range(len(z_agl[:,1,r,c])): 
        member_profs[m] =  profile.create_profile(profile='convective', pres=p[m,:,r,c], hght=z_agl[m,:,r,c], tmpc=t[m,:,r,c], \
                        dwpc=td[m,:,r,c], u=u[m,:,r,c], v=v[m,:,r,c], missing=-9999, strictQC=False, latitude=xlat,date=vdateobj)
        members = {'hght': z_agl[:,:,r,c], 'pres': p[:,:,r,c], 'tmpc': t[:,:,r,c], 'dwpc': td[:,:,r,c], 'u': u[:,:,r,c], 'v': v[:,:,r,c], 'member_profs': member_profs}

    print("end member profs")

    #find summary file lat/lons
    wherelatlon = np.where((sumlats == xlat)&(sumlons == xlon))
    latwhere = wherelatlon[0][0]
    lonwhere = wherelatlon[1][0]
    cape_ml = cape_ml_0[:,latwhere-3:latwhere+4,lonwhere-3:lonwhere+4]
    srh_0to1 = srh_0to1_0[:,latwhere-3:latwhere+4,lonwhere-3:lonwhere+4]

    figname = os.path.join(outdir,'wofs_snd_{:02d}_{:02d}_f{:03d}.png'.format(r+1,c+1,leadminute))
    # pass the data to the plotting script.
    print('rounded lat/lon',round(xlat,2), round(xlon,2))
    plot_wof(prof, members, figname, str(round(xlat,2)), str(round(xlon,2)), idateobj, vdateobj, x_pts=cape_ml, y_pts=srh_0to1)

############## While loop to plot every 3rd output time in forecast #############

timesteps = np.arange(0,fcst_nt+3,3)
process_flag = timesteps * 0 

while (np.min(process_flag)  == 0):
  
   for t_index, t in enumerate(timesteps): 
      if (process_flag[t_index] == 0): 
         temp_sounding = glob.glob(summary_dir+'/'+'wofs_SND_{:02d}_*.nc'.format(t))
         if (len(temp_sounding) > 0): 
            sounding_file = glob.glob(summary_dir+'/'+'wofs_SND_{:02d}_*.nc'.format(t))[0]
            print('sounding file',sounding_file)
            if sounding_file:

               temp_summary= glob.glob(summary_dir+'/'+'news-e_ENV_{:02d}_*.nc'.format(t))
               if (len(temp_summary) > 0): 
                  summary_file= glob.glob(summary_dir+'/'+'news-e_ENV_{:02d}_*.nc'.format(t))[0]

                  if summary_file:

                     summary_file_stats = os.stat(summary_file)
                     if (summary_file_stats.st_size > 390000000):  #hard-coded check on ENV summary file size

                        if (t_index == 0):
                           time.sleep(30)

                        try:                                                 #open WRFOUT file
                           fin = Dataset(sounding_file, "r")
                           print("Opening {} \n".format(sounding_file))
                        except:
                           print("{} does not exist! \n".format(sounding_file))
                           sys.exit(1)

                        process_flag[t_index] = 1

                        #start_date = fin.START_DATE
                        split = os.path.basename(sounding_file).split('_')
                        start_date = split[3] + split[4] 
                        print('start datestr',start_date)
                        idateobj = datetime.datetime.strptime(start_date,'%Y%m%d%H%M')
                        inittimesecs = fin.INIT_TIME_SECONDS
                        if (inittimesecs < 43200): 
                           inittimesecs = inittimesecs + 86400
                        validtimesecs = fin.VALID_TIME_SECONDS
                        diffseconds = validtimesecs - inittimesecs
                        vdateobj = idateobj + datetime.timedelta(seconds=diffseconds)
                        print('start/valid',idateobj,vdateobj)
                        leadminute = int(diffseconds / 60)
#               if (leadminute > 400): #account for day change
#                  leadminute = leadminute - 1440
                        print('lead minute', leadminute)
                        xlats = fin.variables["xlat"][:,:]
                        xlons = fin.variables["xlon"][:,:]
                        u = fin.variables["u"][:,:,:,:]*1.94384
                        v = fin.variables["v"][:,:,:,:]*1.94384
                        p = fin.variables["p"][:,:,:,:]
                        t = fin.variables["t"][:,:,:,:]
                        q = fin.variables["q"][:,:,:,:]
                        z_agl = fin.variables["z_agl"][:,:,:,:]
                        omega = fin.variables["omega"][:,:,:,:]
                        td = thermo.temp_at_mixrat(q*1000., p)

                        nc = Dataset(summary_file,'r')
                        sumlats = nc.variables['xlat'][:,:]
                        sumlons = nc.variables['xlon'][:,:]
                        cape_ml_0 = nc.variables['cape_ml'][:,:,:]
                        srh_0to1_0 = nc.variables['srh_0to1'][:,:,:]
                        nc.close()
                        del nc

                        pool = multiprocessing.Pool(processes=24)

                        indices = []
                        for index,x in np.ndenumerate(xlats[12:16,:]):#[0:1,0:1]):
                           new_index = (index[0] + 12, index[1])
                           indices.append(new_index)

                        pool.map(make_plot,indices)
                        pool.close()
                        pool.join()

   time.sleep(1)

