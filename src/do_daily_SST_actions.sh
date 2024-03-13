#!/bin/bash
# This script runs the download and plot scripts for the WFIP3 SST data
# author: jtomfarrar

# Define the paths to the scripts and directories
download_script="/path/to/download_VIIRS_AVHRR_SST_data.sh"
plot_script="/path/to/plot_WFIP_SST_v0.py"
source_dir="/path/to/source/directory"
dest_dir_data="/path/to/destination/directory"
dest_dir_img="/path/to/destination/directory"

# Run the download script
bash $download_script

# Run the plot script
python $plot_script

# Copy the plots and nc4 files to the shared directory
cp $source_dir/*.png $dest_dir_img
cp $source_dir/*.nc4 $dest_dir_data

