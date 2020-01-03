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

variables = ['uh_2to5', 'uh_0to2', 'wz_0to2', 'comp_dz', 'rain', 'ws_80', 'w_up', 'hail', 'hailcast']
#variables = ['multi', 'uh_2to5', 'uh_0to2', 'wz_0to2', 'comp_dz', 'rain', 'ws_80', 'w_up', 'hail', 'hailcast']

######### Plot swath products #########

ens_t = 0
iteration = 0

while (ens_t < (fcst_nt-1)):
   summary_files_temp = os.listdir(summary_dir)

   for f, file in enumerate(summary_files_temp):
      str_ens_t = str(ens_t)
      if (len(str_ens_t) == 1):
         str_ens_t = '0' + str_ens_t

      if ((file[-28:-25] == 'SWT') and (file[-24:-22] == str_ens_t)):
         time.sleep(2)

         for v, variable in enumerate(variables): 
            cmd = "python news_e_swath.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, ens_t)
            pool.apply_async(run_script, (cmd,))
         ens_t = ens_t + 1
#      else: 
#         time.sleep(2)

#Plot summary products

time.sleep(30)

for v, variable in enumerate(variables):
   cmd = "python news_e_swath.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, ens_t)
   pool.apply_async(run_script, (cmd,))

time.sleep(90)

for v, variable in enumerate(variables):
   cmd = "python news_e_summary.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, fcst_nt)
   pool.apply_async(run_script, (cmd,))

time.sleep(90)

for v, variable in enumerate(variables):
   cmd = "python news_e_swath.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, (fcst_nt-2))
   pool.apply_async(run_script, (cmd,))

for v, variable in enumerate(variables):
   cmd = "python news_e_swath.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, (fcst_nt-1))
   pool.apply_async(run_script, (cmd,))

for v, variable in enumerate(variables):
   cmd = "python news_e_swath.py -d %s -o %s -m %s -n %s -t %d " % (summary_dir, outdir, mapname, variable, fcst_nt)
   pool.apply_async(run_script, (cmd,))

pool.close()
pool.join()


