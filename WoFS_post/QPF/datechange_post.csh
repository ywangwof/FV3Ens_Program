#!/bin/csh -f
set echo

set base_dir = '/work/brian.matilla/python_realtime/'

setenv yyyy `date -u +%Y`
setenv mm `date -u +%m`
setenv dd `date -u +%d `

setenv ddy `date -ud 1-day-ago +%d`
setenv ddt `date -ud 1-day +%d`

set curr_date = $yyyy$mm$dd
set prev_date = $yyyy$mm$ddy
set next_date = $yyyy$mm$ddt

cd $base_dir

exec sed -r -s -i -e "s/$prev_date/$curr_date/g" run*
