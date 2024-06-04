#!/bin/bash

# This script combines the actions of downloading, plotting, and managing SST data

# Author: jtomfarrar, yugao

# Activate the Python environment
# /home/yugao/software/miniforge3/condabin/mamba activate coringa

# Define the paths to the scripts and directories
download_script="/home/yugao/projects/WFIP3_SST/src/download_VIIRS_AVHRR_SST_data.sh"
plot_script="/home/yugao/projects/WFIP3_SST/src/plot_WFIP_SST_v0.py"
img_source_dir="/home/yugao/projects/WFIP3_SST/img/SST_plots/"
dest_dir_img="/home/yugao/projects/WFIP3_SST/img/destination_dir/"
nc_source_dir="/home/yugao/projects/WFIP3_SST/data/processed/SST/good_data/"
dest_dir_data="/home/yugao/projects/WFIP3_SST/data/destination_dir/"

# Run the download script
bash "$download_script"

# Navigate to the directory containing the plot script and run it
cd "/home/yugao/projects/WFIP3_SST/src"
/home/yugao/software/miniforge3/envs/coringa/bin/python3 "plot_WFIP_SST_v0.py"
/home/yugao/software/miniforge3/envs/coringa/bin/python3 "plot_WFIP_HFR_v0.py"
# /home/yugao/software/miniforge3/envs/coringa/bin/python3 "plot_WFIP_HFR_SST_v0.py"
# need to download new files

echo "Script executed from: ${PWD}"

# Deactivate the Python environment
# /home/yugao/software/miniforge3/condabin/mamba deactivate

# Check if there are PNG files in the source directory
if ls "${img_source_dir}"*.png >/dev/null 2>&1; then
    echo "this one is getting called."
    echo "Source directory: $dest_dir_img"
    echo "Destination directory: $dest_dir_img"
    
    # Copy all PNG files to the destination directory
    cp "${img_source_dir}"*.png "$dest_dir_img"
    cp "$img_source_dir"/*.png "$dest_dir_img"
else
    echo "No PNG files found in $img_source_dir"
fi

# Copy the NC4 files to the shared directory (if they exist)
# YG: the files in $nc_source_dir is *.nc, not nc4
if ls "$nc_source_dir"/*.nc >/dev/null 2>&1; then
    echo "this one is getting called."
    echo "Source directory: $nc_source_dir"
    echo "Destination directory: $dest_dir_data"
    cp "$nc_source_dir"/*.nc "$dest_dir_data"

else
    echo "No NC4 files found in $nc_source_dir"
fi
