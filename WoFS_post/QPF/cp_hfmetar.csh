#!/bin/csh

# User defined variables: 

set year = 2019
set month = 05 
set month2 = 05
set day = 17
set day2 = 18

set in_dir = "/work/LDM/MADIS/hfmetar/"$year"/"$month"/"$day"/"
set in_dir2 = "/work/LDM/MADIS/hfmetar/"$year"/"$month2"/"$day2"/"

set out_dir = "/oldscratch/skinnerp/2019_wofs_post/asos/"

# Run wrapper script for forecast: 

cp "${in_dir}"2* "${out_dir}"
cp "${in_dir2}"2* "${out_dir}"

sleep 5

gunzip "${out_dir}"*.gz
gunzip "${out_dir}"*.gz

sleep 10

find "${out_dir}" -name '20190[0-9][0-9][0-9]_[012][0-9][0-9][0-9]' -exec sh -c 'mv $0 $0.nc' {} \;
find "${out_dir}" -name '20190[0-9][0-9][0-9]_[012][0-9][0-9][0-9]' -exec sh -c 'mv $0 $0.nc' {} \;

