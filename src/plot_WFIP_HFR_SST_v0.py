# Plot daily SST fith daily HFR
#
#
# @jtomfarrar, jfarrar@whoi.edu
# Started 2024-03-09

# %%
# make sure the working directory is the same as the directory as this script
import os
os.chdir('/home/jtomf/Python/WFIP3_SST/src')

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
import datetime

# %%
# %matplotlib inline
%matplotlib widget
# %matplotlib qt5
plt.rcParams['figure.figsize'] = (5,4)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 400

savefig = True # set to true to save plots as file

__figdir__ = '../img/HFR_plots/'
sst_figdir = '../img/SST_plots/'
os.system('mkdir  -p ' + __figdir__) #make directory if it doesn't exist
savefig_args = {'bbox_inches':'tight', 'pad_inches':0.2}
plotfiletype='png'

# %% SST plotting function for use with HFR
def plot_sst_v2(ax,sst, levels, xmin, xmax, ymin, ymax, savefig=False):
    day_str=sst.time.dt.strftime("%Y-%m-%d %H:%M").values[0]
    day_str2=sst.time.dt.strftime("%Y-%m-%d").values[0]
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

    ax.set_title('SST + HFR, ' + day_str2 , size = 10.)

    return ax

# %%
########################
# HFR data
HFR_dir = '../data/processed/HFR'
# Load the HFR data
ds2 = xr.open_dataset(HFR_dir + '/NES_HFR_daily.nc')
ds2
########################
# SST data
SST_dir = '../data/processed/SST/good_data'
nc_files = SST_dir + '/VIIRS_AVHRR_SST_merged.nc'
nc_files
ds_sst = xr.open_dataset(nc_files)
ds_sst

 # %%
# make a list of strings for the days in ds2.Time.values (the HFR data)
targets = ds2.Time.dt.strftime('%Y-%m-%d').values
for target_day in targets:
    try:
        sst = ds_sst.sel(time=target_day)
    except:
        print('No SST file for ' + target_day)
        continue # skip rest of loop and move to next target_day
    ###########################
    # Select the HFR data for the same day
    Umean = ds2.East_vel.sel(Time=target_day).squeeze()
    Vmean = ds2.North_vel.sel(Time=target_day).squeeze()
    spd = np.sqrt(Umean**2 + Vmean**2)
    tstr1 = ds2.Time.sel(Time=target_day).dt.strftime('%Y-%m-%d').values

    # create a mask where spd is above a threshold
    mask = spd < 0.2
    Umean = Umean.where(mask, drop=False)
    Vmean = Vmean.where(mask, drop=False)
    spd = spd.where(mask, drop=False)

    lon_min = ds2.Longitude.min().values
    lon_max = ds2.Longitude.max().values
    lat_min = ds2.Latitude.min().values
    lat_max = ds2.Latitude.max().values

    V = [2.5, 7]
    #Set color scale range (define V with input file)
    levels = np.linspace(V[0],V[1],21)
    [xmin, xmax, ymin, ymax] = [-72, -69, 40, 42.5]

    skip = 2
    tol=0.25

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    extent = [lon_min-tol, lon_max+tol,lat_min-tol, lat_max+tol]
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    plot_sst_v2(ax,sst, V, xmin, xmax, ymin, ymax, savefig)
    plt.quiver(ds2.Longitude[0:-1:skip], ds2.Latitude[0:-1:skip], Umean[0:-1:skip,0:-1:skip], Vmean[0:-1:skip,0:-1:skip], scale=5, transform=ccrs.PlateCarree())

    # Export figure
    if savefig:
        plt.savefig(__figdir__+'SST_HFR_'  + target_day + '.' +plotfiletype,**savefig_args)


# %%
# make a movie of the SST and HFR data using ffmpeg
# ffmpeg -framerate 2 -pattern_type glob -i './HFR_plots/*.png' -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -crf 23 -profile:v baseline -level 3.0 -pix_fmt yuv420p -c:a aac -ac 2 -movflags faststart HFR_SST_daily_movie.mp4


#  %%
# mask ds2 to values below .2 m/s
mask = ds2.total_err < 0.15
ds2 = ds2.where(mask, drop=False)
