###################################################################################################

#from mpl_toolkits.basemap import Basemap
import matplotlib
import math
from scipy import *
import pylab as P
import numpy as np
import sys, glob
import os
import time
from optparse import OptionParser
#import netCDF4
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
parser.add_option("-d", dest="summary_dir", type="string", default= None, help="Input Directory (of summary files)")
parser.add_option("-o", dest="outdir", type="string", help = "Output Directory (for summary file)")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.outdir == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   outdir = options.outdir
   fcst_nt = options.fcst_nt

pool = Pool(processes=(24))              # set up a queue to run

##### process swath files from completed ens files

fcst_times = np.arange(fcst_nt)
process = fcst_times * 0
counter = 0

#time.sleep(120) # give post_ens_rain.py a head start

while (np.min(process) == 0):
   counter = counter + 1
   print 'iteration: ', counter, fcst_times
   for t in fcst_times:
      str_t = str(t)
      if (len(str_t) == 1):
         str_t = '0' + str_t

      summary_files_temp = os.listdir(summary_dir)
      summary_files_temp.sort()

      for f, filename in enumerate(summary_files_temp):
         if ((filename[-28:-25] == 'ENS') and (filename[-24:-22] == str_t) and (process[t] == 0)):
            process[t] = 1
            cmd = "python news_e_post_swath_rainfall_v1.py -d %s -t %d -n %d " % (summary_dir, t, fcst_nt)
            pool.apply_async(run_script, (cmd,))
            time.sleep(2)

      time.sleep(5)

#cmd = "python news_e_post_swath_rainfall_v1.py -d %s -t %d -n %d " % (outdir, fcst_nt, fcst_nt)
#pool.apply_async(run_script, (cmd,))

time.sleep(2)

cmd = "python news_e_post_summary_rainfall.py -d %s -t %d -n %d " % (summary_dir, fcst_nt-1, fcst_nt-1)
pool.apply_async(run_script, (cmd,))

print 'final summary command: ', cmd

pool.close()
pool.join()

########################################################
#Deprecated code below from 2018 version
########################################################


#ens_t = 0
#iteration = 0
#
#while (ens_t < fcst_nt):
#   summary_files_temp = os.listdir(summary_dir)
#   summary_files_temp.sort()
#
#   process = []
#   for i in range(0, ens_t): ## Make sure all times prior to the current one are available
#      str_i = str(i)
#      if (len(str_i) == 1):
#         str_i = '0' + str_i
#
#      for f, file in enumerate(summary_files_temp):
#         if ((file[-28:-25] == 'ESR') and (file[-24:-22] == str_i)):
#            process.append(i)
#
#   for f, file in enumerate(summary_files_temp):  ## If new time available, process.
#      str_ens_t = str(ens_t)
#      if (len(str_ens_t) == 1):
#         str_ens_t = '0' + str_ens_t
#
#      if ((file[-28:-25] == 'ESR') and (file[-24:-22] == str_ens_t) and (len(process) >= ens_t)):
#         time.sleep(1) #give file time to write ... just in case
#         cmd = "python news_e_post_swath_rainfall.py -d %s -t %d -n %d " % (outdir, ens_t, fcst_nt)
#         pool.apply_async(run_script, (cmd,))
#         ens_t = ens_t + 1
#
##get last fcst timestep
#
#time.sleep(90)
#
#cmd = "python news_e_post_swath_rainfall.py -d %s -t %d -n %d " % (outdir, fcst_nt, fcst_nt)
#print 'final swath rainfall command: ', cmd
#pool.apply_async(run_script, (cmd,))
#
#cmd = "python news_e_post_summary_rainfall.py -d %s -t %d -n %d " % (outdir, fcst_nt, fcst_nt)
#print 'final summary command: ', cmd
#pool.apply_async(run_script, (cmd,))
#
#pool.close()
#pool.join()


