#!/bin/csh
#==================================================================

source ~/.tcshrc
set echo

# cd to directory where job was submitted from
cd /scratch2/patrick.skinner/python_2017_post

# User defined variables: 

set dates = (20170910)

#set times = (2300)
set times = (1900 2000 2030 2100 2130 2200 2230 2300 2330 0000 0030 0100 0130 0200 0230 0300)
#set times = (2330 0000 0030 0100 0130 0200 0230 0300)
#set nt = (19 37 19 37 19 37 19 36)
set nt = (49 37 19 37 19 37 19 37 19 37 19 37 19 37 19 36)
set ne = 18

# Run wrapper script for forecast: 

sleep 5

foreach date ($dates)

   set fcst = "/scratch/tajones/gsi/"$date"/GSI_RLT/FCST/"
   set summary = "/work1/skinnerp/2017_newse_post/summary_files_test/"$date"_gsi/"

   set i = 1

   foreach dir ($times)
      set tempnt = $nt[$i]
      echo $tempnt

      set tempsumdir = $summary$dir"/"
      echo $tempsumdir

      set fcstdir = $fcst$dir"/"
      echo $fcstdir

      if ( -d $fcstdir ) then
         python wrapper_post.py -d $fcstdir -o $tempsumdir -e $tempnt
      endif

      @ i++
   end
end

sleep 10
                                                                         


