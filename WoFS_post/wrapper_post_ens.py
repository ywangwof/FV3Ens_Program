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
ne = 40
member_dirs = []

while (len(member_dirs) < ne):         #wait for all members to be running
   member_dirs = []
   member_dirs_temp = os.listdir(fcst_dir)
   for d, dir in enumerate(member_dirs_temp):
      if (dir[0:4] == 'mem_'):
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

while (current_t < fcst_nt-1):
#for t in range(0, fcst_nt):
   process = 1

##### process swath files from completed ens files

   for n in range(0, len(member_dirs)):
      temp_dir = os.path.join(fcst_dir, member_dirs[n],'summary')

      member_files = []
      temp_files = os.listdir(temp_dir)

      for f, file in enumerate(temp_files):
         if (file[0:7] == 'FV3_ENS'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
#         if (file[0:6] == 'wrfwof'):                               #assumes filename format of: "wrfout_d02_yyyy-mm-dd_hh:mm:ss
            member_files.append(file)
      #if (len(member_files) < (current_t+1)):
      #   process = 0

   if (process == 1):
      if (current_t >= 3):
         cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-3))
         #print cmd
         pool.apply_async(run_script, (cmd,))
      current_t = current_t + 1
      iteration = 0
   else:
      iteration = iteration + 1
      iter_time = iteration * 5.
      time.sleep(5)
      print 'iterating ... ', iter_time, current_t

time.sleep(10)

#get last fcst timestep

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-3))
pool.apply_async(run_script, (cmd,))

time.sleep(10)

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-2))
pool.apply_async(run_script, (cmd,))

time.sleep(10)

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t-1))
pool.apply_async(run_script, (cmd,))

time.sleep(30)

cmd = "python news_e_post_ensemble.py -d %s -o %s -t %d " % (fcst_dir, outdir, (current_t))
print 'last post_ensemble command', current_t, cmd

pool.apply_async(run_script, (cmd,))


pool.close()
pool.join()


