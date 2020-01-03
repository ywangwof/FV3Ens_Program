#!/bin/csh
#==================================================================

source ~/.tcshrc
set echo

# cd to directory where job was submitted from
cd /scratch2/patrick.skinner/python_2017_plotting

# User defined variables: 

set dates = (20170516)
#set dates = (20170508 20170509 20170511 20170515 20170516 20170517 20170518 20170519 20170522 20170523 20170524 20170525 20170526 20170527 20170531 20170601 20170602 20170605 20170606 20170607 20170608)

#set times = (2300)
set times = (1900 2000 2030 2100 2130 2200 2230 2300 2330 0000 0030 0100 0130 0200 0230 0300)
set ne = 18

# Build directories needed: 

set asos = "/work1/skinnerp/2017_newse_post/asos/"
set rot_qc = "/work1/skinnerp/MRMS_verif/rot_qc_nsslmod3/"
# Run wrapper script for forecast: 

sleep 5

foreach date ($dates)

#   set summary = "/work1/Thomas.Jones/2017_reruns/summary_files/"$date"/NSSLMOD2/"
   set summary = "/work1/skinnerp/2017_newse_post/summary_files_test/"$date"_COMBO1/"
   set image = "/www/www.nssl.noaa.gov/projects/wof/news-e/newse_images/"$date"_COMBO1/"
   set mrms = "/work1/skinnerp/MRMS_verif/mrms_cressman/"$date"/"

   foreach dir ($times)
      set tempdir = $image$dir"/"
      echo $tempdir

      set tempsumdir = $summary$dir"/"
      echo $tempsumdir

      set rotfile = $rot_qc$date"_"$dir"_rotqc.nc"
      echo $rotfile

      set tempfiles=`ls -a ${tempsumdir} | wc | awk '{print $1}'`

      if ( "${tempfiles}" >= $ne ) then
         sleep 2
         python wrapper_plot.py -d $tempsumdir -i $tempdir -n $ne
         sleep 2
#         python wrapper_hwt.py -d $tempsumdir -i $tempdir -n $ne
#         sleep 2
         python wrapper_paintball.py -d $tempsumdir -i $tempdir -n $ne
         sleep 2
         python wrapper_dot.py -d $tempsumdir -a $asos -i $tempdir -n $ne
         sleep 2
         python wrapper_rot_object.py -d $tempsumdir -o $tempdir -i $tempdir -m $mrms -n $ne
         sleep 2
#         python wrapper_obj.py -r $rotfile -i $tempdir
      endif
   end
end

sleep 10
                                                                         


