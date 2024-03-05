#!/bin/bash
# Sample download script obtained from https://eastcoast.coastwatch.noaa.gov/ecn_data_download_sample.txt
# ecn_data_download_sample.sh

# Sample script to download files from CoastWatch East Coast Node's http server
# (modify as needed for your data files of interest)

# - uses wget to download MODIS chlorophyll daily composite HDF file
#   for NE region from CoastWatch East Coast Node http server
# - downloads all files not already downloaded

# East Coast Node directory tree structure is of the form:
#     http://eastcoast.coastwatch.noaa.gov/data/SENSOR/PRODUCT/INTERVAL/REGION/index.php

# Ron Vogel, SMRC for NOAA CoastWatch East Coast Node
# May 31, 2013



# get index.php filelisting from East Coast Node http data directory
# wget http://eastcoast.coastwatch.noaa.gov/data/modis/chl/daily/ne/index.php
wget https://eastcoast.coastwatch.noaa.gov/data/avhrr-viirs/sst-ngt/daily/mr/index.php

# grep for needed filetype (in style of grep -e)
#grep -e 'MODSCW_[0-9]\{7\}_DAILY_AQUA_CHLORA_NE_1KM.hdf' index.php > filelist.txt
grep -e 'ACSPOCW_[0-9]\{7\}_DAILY_MULTISAT_SST-NGT_MR_750M.nc4' index.php > filelist.txt

# alternative grep if downloading swath data (e.g. MODCYCW for Chesapeake Bay)
#grep -e 'MODCYCW_[0-9]\{7\}_[0-9]\{4\}_AQUA_CHLORA_CD_1KM.hdf' index.php > filelist.txt

# remove extraneous html to generate clean list of filenames to download
#filenames="$(awk '{print substr($3, 12, index($3,"hdf") - 9)}' filelist.txt)"
filenames="$(awk '{print substr($3, 13, index($3,"nc4") - 10)}' filelist.txt)"
YOUR_DIRECTORY="../data/external/SST"

for file in $filenames ; do
    echo Starting $file
   # get all files not already downloaded
   if [[ !( -f $YOUR_DIRECTORY/$file) ]]; then
      echo Getting file:
      echo $file
#      wget http://eastcoast.coastwatch.noaa.gov/data/modis/chl/daily/ne/$file
      wget https://eastcoast.coastwatch.noaa.gov/data/avhrr-viirs/sst-ngt/daily/mr/$file
      wait

      mv $file $YOUR_DIRECTORY/$file
   fi
done

# change name of filelist.txt to include today's date
mv filelist.txt $YOUR_DIRECTORY/filelist_VIIRS_AVHRR_SST_$(date +%Y%m%d).txt

# rm filelist.txt
rm index.php
