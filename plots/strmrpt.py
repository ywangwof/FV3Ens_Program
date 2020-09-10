#!/usr/bin/env python
## ---------------------------------------------------------------------
##
## Download SPC storm report and plot it in a map
##
## ---------------------------------------------------------------------
##
## HISTORY:
##
##   Yunheng Wang (04/17/2018)
##   Initial version.
##
##
########################################################################
##
## Requirements:
##
##   o Python 2.7 or above
##
########################################################################

import os, sys, re
from datetime import datetime, timedelta

import urllib2
import logging
import csv

########################################################################

def download_srpt(wrkdir,filename):
    try:
        logging.info('Starting downloading storm report file %s ...'%filename)
        url1="http://www.spc.noaa.gov/climo/reports"

        url="%s/%s" % (url1,filename)

        output=os.path.join(wrkdir,filename)

        tCHUNK = 512*1024
        downstr = ' '
        try :
          req = urllib2.urlopen(url, None, 20)
          with open(output,'wb') as fp :
            while True :
              fchunk = req.read(tCHUNK)
              if not fchunk : break
              fp.write(fchunk)
            downstr = ' ... '
        except urllib2.HTTPError as e:
          logging.error('HTTP Error with code: %d.' % e.code )
        except urllib2.URLError as e :
          logging.error('URL Error with reason: %s.' % e.reason )

        ##fsize = os.path.getsize(output)
        ##if fsize < expectfs :
        ##    logging.info('%s size (%d) less than %dK.'% (output,fsize,expectfs/1024))
        ##else :
        ##    logging.info('%s%s(%d)' % (output,downstr,fsize) )

    except Exception,x:
        logging.error('Error:download SPC storm report file: %s' % x)
        raise x
        #sys.exit(1)

########################################################################

def decode_srptfile(filename,starttime=None,endtime=None):
  '''
     starttime and endtime are integers
     if the time is before 10:00, add 1 before it. For examples
     21:00  -> 2100
     03:00  -> 10300
  '''
  srptdict = {'tornado': [], 'hail': [], 'wind': [] }

  #tmatch="Time,F_Scale,Location,County,State,Lat,Lon,Comments"
  #hmatch="Time,Size,Location,County,State,Lat,Lon,Comments"
  #wmatch="Time,Speed,Location,County,State,Lat,Lon,Comments"

  with open(filename, 'rb') as csvfile:
      reader = csv.reader(csvfile)
      for row in reader:
          if row[1] == "F_Scale":    # tornado lines begin
              stype = 'tornado'
          elif row[1] == 'Size':     # hail lines begin
              stype = 'hail'
          elif row[1] == 'Speed':    # wind lines begin
              stype = 'wind'
          else:                        # data
              if row[0] < "1000":
                timeval = 10000+int(row[0])
              else:
                timeval = int(row[0])

              if starttime is not None:
                  if timeval < starttime:  continue
              if endtime is not None:
                  if timeval > endtime:    continue
              #srptdict[stype].append((timeval,float(row[5]),float(row[6])))
              srptdict[stype].append((row[0],float(row[5]),float(row[6])))

  return srptdict

########################################################################

def plot_srpt(mymap,srpts):
  ''' mps is an instance of basemap
  '''

  markers = { 'tornado' : 'r^', 'hail': 'gD', 'wind': 'b>'}

  for key in srpts.keys():
    #print key, srpts[key]
    x, y = mymap([x[2] for x in srpts[key]],[x[1] for x in srpts[key]])
    mymap.plot(x, y, markers[key], markersize=4)

########################################################################

def plot_storms(mymap,eventdate,wrkdir):

  srptdtstr = datetime.strptime(eventdate, '%Y%m%d').strftime('%y%m%d')
  srptfile  = '%s_rpts.csv'%srptdtstr

  wrkrptfile = os.path.join(wrkdir,srptfile)
  if not os.path.lexists(wrkrptfile):
      try:
        download_srpt(wrkdir,srptfile)
      except:
        logging.error('ERROR: no storm reports will be plotted.')
        return

  #srpts = decode_srptfile(wrkrptfile,2100,10200)
  srpts = decode_srptfile(wrkrptfile)

  plot_srpt(mymap,srpts)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == "__main__":

  eventdate = "20170509"
  wrkdir = "/scratch/ywang/test_runs/HWT2017_11"

  srptdtstr = datetime.strptime(eventdate, '%Y%m%d').strftime('%y%m%d')
  srptfile  = '%s_rpts.csv'%srptdtstr

  wrkrptfile = os.path.join(wrkdir,srptfile)
  if not os.path.lexists(wrkrptfile):
      download_srpt(wrkdir,srptfile)

  srpts = decode_srptfile(wrkrptfile,2100,10200)

  for key,value in srpts.iteritems():
      print key, '->', len(value)
