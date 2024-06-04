# Script to plot SST for the WFIP3 region
#
# @jtomfarrar, jfarrar@whoi.edu
# Started 2024-03-04

# %%
# make sure the working directory is the same as the directory as this script
import os
os.chdir('/home/yugao/projects/WFIP3_SST/src')

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import glob

import matplotlib
from matplotlib.patheffects import Stroke
import copy
import cartopy.crs as ccrs                   # import projections
import cartopy
# Data need to be downloaded using the bash script download_VIIRS_AVHRR_SST_data.sh

# This may be useful for later; GitHub Copilot suggested it
#data_dir = '/Users/jfarrar/Google Drive File Stream/My Drive/Projects/WFIP3/VIIRS_AVHRR_SST_data'

# %%
# %matplotlib inline
# %matplotlib widget
# %matplotlib qt5
plt.rcParams['figure.figsize'] = (5,4)
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 400

savefig = False # set to true to save plots as file

__figdir__ = '../img/SST_plots/'
os.system('mkdir  -p ' + __figdir__) #make directory if it doesn't exist
savefig_args = {'bbox_inches':'tight', 'pad_inches':0.2}
plotfiletype='png'

# %%
data_dir = '../data/external/SST'
output_dir = '../data/processed/SST/good_data'
os.system('mkdir  -p ' + output_dir) #make directory if it doesn't exist
# get a list of *.nc4 files in the data directory
nc_files = glob.glob(data_dir + '/*.nc4')
print('get a list of *.nc4 files in the data directory')
print(nc_files)

# %%
# See if xr.open_mfdataset can open all the files at once
ds = xr.open_mfdataset(nc_files)
ds

# %%
# I know that 2024-02-18 had valid data; try to select that time.  
# ds_sub = ds.sel(time='2024-04-28')
# 
ds_sub = ds.isel(time=1)
# Select data by index instead of date because the dates in the files may shift

# %%
# Plot the SST
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent([-72, -69, 40, 42.5])
plt.pcolormesh(ds_sub.lon,ds_sub.lat,ds_sub.sst.squeeze(),transform=ccrs.PlateCarree())
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.colorbar(label='SST (C)')
datestr = ds_sub.time.dt.strftime('%Y-%m-%d').values
plt.title('Sea Surface Temperature, ' + datestr)


# %%
V = [2, 8]
#Set color scale range (define V with input file)
levels = np.linspace(V[0],V[1],21)
sst = ds_sub

[xmin, xmax, ymin, ymax] = [-72, -69, 40, 42.5]


# %%

def plot_sst(sst, levels, xmin, xmax, ymin, ymax, savefig=False):
    plt.figure()
    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=290))  # Orthographic
    extent = [xmin, xmax, ymin, ymax]
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    day_str=sst.time.dt.strftime("%Y-%m-%d %H:%M").values
    day_str2=sst.time.dt.strftime("%Y-%m-%d").values
    #plt.set_cmap(cmap=plt.get_cmap('nipy_spectral'))
    plt.set_cmap(cmap=plt.get_cmap('turbo'))
    cmap = copy.copy(matplotlib.cm.get_cmap("turbo"))
    cmap.set_bad(color='white')
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    coast = cartopy.feature.GSHHSFeature(scale="full")
    ax.add_feature(coast, zorder=3, facecolor=[.6,.6,.6], edgecolor='black')
    # arr = np.ma.array(sst.sea_surface_temperature, mask=(sst.quality_level == 1))
    #cs = ax.pcolormesh(sst.lon,sst.lat,np.squeeze(arr)-273.15, vmin=levels[0], vmax=levels[-1], transform=ccrs.PlateCarree())
    cs = ax.pcolormesh(sst.lon,sst.lat,np.squeeze(sst.sea_surface_temperature)-273.15, vmin=levels[0], vmax=levels[-1], transform=ccrs.PlateCarree())
    cb = plt.colorbar(cs,fraction = 0.022,extend='both')
    cb.set_label('SST [$\circ$C]',fontsize = 10)

    ax.set_title('SST, ' + day_str , size = 10.)


    # Export figure
    if savefig:
        plt.savefig(__figdir__+'SST_map_'  + day_str2 + '.' +plotfiletype,**savefig_args)
    # plt.show(block=False)

# %%
plot_sst(ds_sub, V, xmin, xmax, ymin, ymax, savefig)

# %%
# Now loop through all the times in ds and make a plot for each one
# Also, check the portion of the data that is not nan

# first check the portion of the data that is not nan by making a variable that is 1 if 
# the data is not nan and 0 if it is nan
sst_not_nan = np.isfinite(ds.sst.squeeze().values).astype(int)

# I wanted to choose a region to check for NaNs, but it was 
# complicated so jsut check over whole domain
frac_not_nan = sst_not_nan.mean(axis=1).mean(axis=1)
plt.figure()
plt.plot(ds.time,frac_not_nan)
# Make time axis look nicer
plt.gca().xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))
# make time labels look nicer
plt.xticks(rotation=45)
# add minor ticks for each day
plt.gca().xaxis.set_minor_locator(matplotlib.dates.DayLocator())
plt.title('Fraction of data that is not nan')


# %%
# Now loop through all the times in ds and make a plot for each one that has frac_not_nan>0.5
# Also, when frac_not_nan>0.5, copy the data to a "good data" folder
for i in range(len(ds.time)):
    if frac_not_nan[i] > 0.5:
        ds_sub = ds.isel(time=i)
        plot_sst(ds_sub, V, xmin, xmax, ymin, ymax, savefig=True)
        # copy the data to a "good data" folder
        ds_sub.to_netcdf(output_dir + '/' + ds_sub.time.dt.strftime('%Y-%m-%d_%H%M%S').values + '.nc')

# %%
# Add frac_not_nan as a variable in ds
ds['frac_not_nan'] = xr.DataArray(frac_not_nan, dims=['time'], coords={'time':ds.time})
# add description
ds['frac_not_nan'].attrs['description'] = 'Fraction of data that is not nan'
# Now write ds to a netcdf file
ds.to_netcdf(output_dir + '/VIIRS_AVHRR_SST_merged.nc')

# %%
print("finished running plot_WFIP_SST_v0.py")
# %%
