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
parser.add_option("-f", dest="summary_dir", type="string", help = "Summary file directory")
parser.add_option("-o", dest="image_dir", type="string", help = "Image directory")

(options, args) = parser.parse_args()

if ((options.summary_dir == None) or (options.image_dir == None)):
   print
   parser.print_help()
   print
   sys.exit(1)
else:
   summary_dir = options.summary_dir
   image_dir = options.image_dir

pool = Pool(processes=(2))              # set up a queue to run

######### Check to see if all forecasts from all members are present: #########

fcst_times = ['1900', '2000', '2030', '2100', '2130', '2200', '2230', '2300', '2330', '0000', '0030', '0100', '0130', '0200', '0230', '0300']
fcst_timesteps = [48, 36, 18, 36, 18, 36, 18, 36, 18, 36, 18, 36, 18, 36, 18, 36]

for f, fcst in enumerate(fcst_times): 
   temp_summary_dir = os.path.join(summary_dir, fcst)
   temp_summary_dir = temp_summary_dir + '/'

   temp_image_dir = os.path.join(image_dir, fcst)
   temp_image_dir = temp_image_dir + '/'

   cmd = "python wrapper_swath.py -d %s -o %s -e %d " % (temp_summary_dir, temp_image_dir, fcst_timesteps[f])
   pool.apply_async(run_script, (cmd,))

#   cmd = "python wrapper_timestep.py -d %s -o %s -e %d " % (temp_summary_dir, temp_image_dir, fcst_timesteps[f])
#   pool.apply_async(run_script, (cmd,))

pool.close()
pool.join()


