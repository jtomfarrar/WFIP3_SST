# Manipulate/plot MVCO NES HFR data
#
#
# @jtomfarrar, jfarrar@whoi.edu
# Started 2024-03-04

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

savefig = False # set to true to save plots as file

__figdir__ = '../img/HFR_plots/'
os.system('mkdir  -p ' + __figdir__) #make directory if it doesn't exist
savefig_args = {'bbox_inches':'tight', 'pad_inches':0.2}
plotfiletype='png'

# %%
data_dir = '/mnt/d/tom_data/WFIP3/NES_HFR_data/NES2_NCout' #'../data/external/SST'
output_dir = '../data/processed/HFR'
os.system('mkdir  -p ' + output_dir) #make directory if it doesn't exist
# get a list of *.nc4 files in the data directory
nc_files = glob.glob(data_dir + '/*.nc')
nc_files


# %%
ds = xr.open_mfdataset(
    nc_files,
    concat_dim="Time",
    combine="nested",
)

# %%
# The dimension 'Time' is a series of integers and not a datetime64 variable.
# We need to replace the Dimension 'Time' with the coordinate 'datetime'
ds['Time'] = ds.datetime.squeeze()

# %%
# squeeze the data to remove the singleton dimension 't'
ds = ds.squeeze()

# %%
# Sort the data by time
ds = ds.sortby('datetime')

# %%
# Mask bad values (999) in ds.East_vel and ds.North_vel
ds['East_vel'] = ds.East_vel.where(ds.East_vel != 999, drop=False)
ds['North_vel'] = ds.North_vel.where(ds.North_vel != 999, drop=False)
ds = ds.where(ds.total_err < 0.25 , drop=False)

# %%
# Re-grid the data to daily averages
ds2 = ds.resample(Time='D').mean()

# %% 
# Modify attributes to note that daily averages were computed from 30-minute data
ds2.attrs['History'] = ds2.attrs['History']+ '; daily averages computed from 30-minute data by T. Farrar'
# Modify ds2.attrs['CREATION_DATE'] to the current date
# find today's date:
today = datetime.date.today()
ds2.attrs['CREATION_DATE'] = today.strftime('%Y-%m-%d')

# Save daily data to a netCDF file
ds2.to_netcdf(output_dir + '/NES_HFR_daily.nc')

# %%
# Reloading the data set clears out some wierd behavior from dask/xarray 
# that keeps everything as 4-D dask arrays despite a singleton dimension
ds2 = xr.open_dataset(output_dir + '/NES_HFR_daily.nc')
ds2
# %%

# Umean = ds2.East_vel.mean(dim='Time',skipna=True).squeeze()
# Vmean = ds2.North_vel.mean(dim='Time',skipna=True).squeeze()

# Select a particular time range
'''
ds = ds2.sel(Time=slice('2023-11-01', '2023-11-30'))
Umean = ds.East_vel.mean(dim='Time',skipna=True).squeeze()
Vmean = ds.North_vel.mean(dim='Time',skipna=True).squeeze()
'''

t = 10
# ds2 = ds2.where(ds2.total_err < 0.15 , drop=False)

Umean = ds2.East_vel.isel(Time=t).squeeze()
Vmean = ds2.North_vel.isel(Time=t).squeeze()
spd = np.sqrt(Umean**2 + Vmean**2)
tstr1 = ds2.Time[t].dt.strftime('%Y-%m-%d').values

# %%
# create a mask where spd is above a threshold
mask = spd < 0.2
Umean = Umean.where(mask, drop=False)
Vmean = Vmean.where(mask, drop=False)
spd = spd.where(mask, drop=False)

# %%
lon_min = ds2.Longitude.min().values
lon_max = ds2.Longitude.max().values
lat_min = ds2.Latitude.min().values
lat_max = ds2.Latitude.max().values

# %%
plt.figure(figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree())
extent = [lon_min, lon_max,lat_min, lat_max]
ax.set_extent(extent, crs=ccrs.PlateCarree())
V = [0, .2]
vmin, vmax = V
plt.pcolormesh(ds2.Longitude, ds2.Latitude, spd, vmin=vmin, vmax=vmax, cmap='turbo',transform=ccrs.PlateCarree())
plt.colorbar(label='Current speed (m/s)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.quiver(ds2.Longitude, ds2.Latitude, Umean, Vmean, scale=5, transform=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

coast = cartopy.feature.GSHHSFeature(scale="high")
ax.add_feature(coast, zorder=3, facecolor=[.6,.6,.6], edgecolor='black')

title_str = tstr1 
#title_str=title_str[0]
plt.title(title_str)
if savefig:
    plt.savefig(__figdir__ + 'NES_HFR_current_speed_' + title_str + '.' + plotfiletype, **savefig_args)



# %%
