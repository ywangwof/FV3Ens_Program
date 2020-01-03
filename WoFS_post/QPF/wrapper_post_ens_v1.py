###################################################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
import pylab as P
import numpy as np
import sys, glob
import os
import time
from optparse import OptionParser
import netCDF4
from news_e_post_cbook import *

from multiprocessing import Pool

sys.path.append("/scratch/software/Anaconda2/bin")

###################################################################################################
# run_script is a function that runs a system command

def run_script(cmd):
    print "Executing command:  " + cmd
    os.system(cmd)
    print cmd + "  is finished...."
    return

###################################################################################################

parser = OptionParser()
parser.add_option("-d", dest="fcst_dir", type="string", default= None, help="Input Directory (of wrfouts)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for summary file)")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.fcst_dir == None) or (options.outdir == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   fcst_dir = options.fcst_dir
   outdir = options.outdir
   fcst_nt = options.fcst_nt

pool = Pool(processes=(24))              # set up a queue to run

######### Check to see if all forecasts from all members are present: #########

ne = 18
member_dirs = []

while (len(member_dirs) < ne):         #wait for all members to be running
   member_dirs = []
   member_dirs_temp = os.listdir(fcst_dir)
   for d, dir in enumerate(member_dirs_temp):
      if (dir[0:3] == 'ENS'):
         member_dirs.append(dir)
   print 'NUMBER OF MEMBER DIRS: ', len(member_dirs)
   if (len(member_dirs) == ne):
      print 'found all member dirs ... '
      time.sleep(1)
   else:
      time.sleep(10)

member_dirs.sort()

######### Post it. #########

current_t = 0
ens_t = 0
iteration = 0

while (current_t < fcst_nt): 
#for t in range(0, fcst_nt):
   process = 1

##### process swath files from completed ens files

   summary_files_temp = os.listdir(outdir)

   for f, file in enumerate(summary_files_temp):
      str_ens_t = str(ens_t) 
      if (len(str_ens_t) == 1): 
         str_ens_t = '0' + str_ens_t

      if ((file[-28:-25] == 'ENS') and (file[-24:-22] == str_ens_t)):
         cmd = "python news_e_post_swath.py -d %s -t %d -n %d " % (outdir, ens_t, fcst_nt)
         pool.apply_async(run_script, (cmd,))
         ens_t = ens_t + 1

   for n in range(0, len(member_dirs)):
      temp_dir = os.path.join(fcst_dir, member_dirs[n])

      member_files = []
      temp_files = os.listdir(temp_dir)

      for f, file in enumerate(temp_files):
         if (file[0:9] == 'wrfout_d0'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
            member_files.append(file)
      if (len(member_files) < (current_t+1)): 
         process = 0
     
   if (process == 1): 
      if (current_t >= 3): 
         cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-3))
         pool.apply_async(run_script, (cmd,))
      current_t = current_t + 1
      iteration = 0
   else: 
      iteration = iteration + 1
      iter_time = iteration * 5.
      time.sleep(5)
      print 'iterating ... ', iter_time, current_t
   print '1st while loop: ', file, current_t
   time.sleep(1)

time.sleep(10)

#get last fcst timestep

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-3))
pool.apply_async(run_script, (cmd,))

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-2))
pool.apply_async(run_script, (cmd,))

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-1))
pool.apply_async(run_script, (cmd,))

time.sleep(20)

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t))
print 'last post_ensemble command', current_t, cmd

pool.apply_async(run_script, (cmd,))

##### process swath files from completed ens files

while (ens_t < fcst_nt):
   time.sleep(1)
   summary_files_temp = os.listdir(outdir)
   summary_files_temp.sort()
   print summary_files_temp[-1] 
   for f, file in enumerate(summary_files_temp):
      str_ens_t = str(ens_t)
      if (len(str_ens_t) == 1):
         str_ens_t = '0' + str_ens_t

      print '2nd while loop: ', file, ens_t, str_ens_t
      if ((file[-28:-25] == 'ENS') and (file[-24:-22] == str_ens_t)):
         print 'ENS T: ', ens_t
         cmd = "python news_e_post_swath.py -d %s -t %d -n %d " % (outdir, ens_t, fcst_nt)
         pool.apply_async(run_script, (cmd,))
         ens_t = ens_t + 1

#finish = 0
print 'CURRENT T; ENS T: ', current_t, ens_t
#while (finish == 0):
#summary_files_temp = os.listdir(outdir)

#for f, file in enumerate(summary_files_temp):
#   print 'asdf: ', file
#   if ((file[-28:-25] == 'ENS') and (file[-24:-22] == '36')):

cmd = "python news_e_post_swath.py -d %s -t %d -n %d " % (outdir, ens_t, fcst_nt)
print 'final swath command: ', cmd
pool.apply_async(run_script, (cmd,))

cmd = "python news_e_post_summary.py -d %s -t %d -n %d " % (outdir, ens_t, fcst_nt)
print 'final summary command: ', cmd
pool.apply_async(run_script, (cmd,))
 
#      finish == 1
#   else: 
#         time.sleep(10)

pool.close()
pool.join()


