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
parser.add_option("-m", dest="mapname", type="string", help = "Path to pickled basemap for plotting")
parser.add_option("-e", dest="fcst_nt", type="int", help = "Total number of timesteps in forecast")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.outdir == None) or (options.mapname == None) or (options.fcst_nt == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   outdir = options.outdir
   mapname = options.mapname
   fcst_nt = options.fcst_nt

pool = Pool(processes=(24))              # set up a queue to run

#variables = ['uh_2to5', 'uh_0to2', 'wz_0to2', 'comp_dz', 'rain', 'ws_80', 'w_up', 'hail', 'hailcast']
variables = ['comp_dz', 'rain']

######### Plot swath rainfall products #########
fcst_times = np.arange(fcst_nt)
process = fcst_times * 0
counter = 0

while (np.min(process) == 0):
   counter = counter + 1
   print 'iteration: ', counter
   for t in fcst_times:
      str_t = str(t)
      if (len(str_t) == 1):
         str_t = '0' + str_t

      summary_files_temp = os.listdir(summary_dir)
      summary_files_temp.sort()

      for f, file in enumerate(summary_files_temp):
         if ((file[-28:-25] == 'PCP') and (file[-24:-22] == str_t) and (process[t] == 0)):
            process[t] = 1
            for v, variable in enumerate(variables):
               cmd = "python news_e_swath_rainfall.py -d %s -o %s -m %s -n %s -t %d -f %d " % (summary_dir, outdir, mapname, variable, t, fcst_nt)
               pool.apply_async(run_script, (cmd,))

   time.sleep(1)

# Plot summary rain products

print("Now plotting summary rain products")

time.sleep(10)

for v, variable in enumerate(variables):
   cmd = "python news_e_summary_rainfall.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, fcst_nt)
   pool.apply_async(run_script, (cmd,))

time.sleep(5)

pool.close()
pool.join()


################################################
# Deprecated code from the 2018 Version
################################################

######### Plot swath products #########

#ens_t = 0
#iteration = 0
#
#while (ens_t < fcst_nt):
#   summary_files_temp = os.listdir(summary_dir)
#
#   for f, file in enumerate(summary_files_temp):
#      str_ens_t = str(ens_t)
#      if (len(str_ens_t) == 1):
#         str_ens_t = '0' + str_ens_t
#
#      if ((file[-28:-25] == 'PCP') and (file[-24:-22] == str_ens_t)):
#         time.sleep(2)
#
#         for v, variable in enumerate(variables):
#            cmd = "python news_e_swath_rainfall.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, ens_t)
#            pool.apply_async(run_script, (cmd,))
#         ens_t = ens_t + 1
##      else:
##         time.sleep(2)
#
##Plot summary products
#
#time.sleep(90)
#
#for v, variable in enumerate(variables):
#   cmd = "python news_e_summary_rainfall.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, fcst_nt)
#   pool.apply_async(run_script, (cmd,))
#
#for v, variable in enumerate(variables):
#   cmd = "python news_e_swath_rainfall.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, fcst_nt)
#   pool.apply_async(run_script, (cmd,))
#
#pool.close()
#pool.join()


