# %%
# Figures for analog paper, ERA5
#
#
# Original: vikki.thompson 26/07/2023
# Last Editted 22 May 2024 (for revisions 1)


# %%
### Load neccessary libraries
import subprocess
import numpy as np
import iris
import sys
import matplotlib.pyplot as plt
sys.path.append('/usr/people/thompson/WP1')
import functions_get_data as gdata
import functions_plot_data as pdata
import scipy.stats as stats
import cartopy.crs as ccrs
import cartopy as cart
import cartopy.feature as cf
import glob
plt.ion(); plt.show()

# %%
## Variables
R1 = [70, 30, 30, -30] # analog region
R2 = [62, 42, 20, -10] # Hylke's region
R3 = [80, 20, 180, -180]
date = [2021, 'Jul', 14] # event date
##

## Analogues
region = R1
past_Y1 = 1950
past_Y2 = 1980
present_Y1 = 1993
present_Y2 = 2023

# %%
# Get data, find analogues based on streamfunction
psi_event = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2021,2023))), region)
event = gdata.pull_out_day_era(psi_event, date[0], date[1], date[2])
psi_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), region)
psi_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), region)
dates_past = gdata.analogs_datelist(psi_past, event)[:30] # find analogs
dates_present = gdata.analogs_datelist(psi_present, event)[:30]

# Composites of streamfunction for region
comp_past = gdata.composite_dates(psi_past, dates_past) 
comp_present = gdata.composite_dates(psi_present, dates_present)
#pdata.composite_plots(event, comp_past, comp_present)
#pdata.all_analogs(event, psi_past, dates_past) # plot past 30 analogs
#pdata.all_analogs(event, psi_present, dates_present) # plot present 30 analogs

# %%
yr = 2021
file = '/net/pc200023/nobackup/users/thompson/ERA5/tw/ERA5_tw_day_'+str(yr)+'*.nc' # numbers indicate months extracted
cube = iris.load(file)[0]
#iris.coord_categorisation.add_year(cube, 'time')
##\MyConstraint = iris.Constraint(year=lambda y: years[0] <= y <= years[-1])
#all_yr_list = cube.extract(MyConstraint)

# %%
import iris.util
from iris.util import equalise_attributes
from iris.util import unify_time_units

def era5_tw(years):
    all_yr_list = iris.cube.CubeList()
    for yr in years:
        files = '/net/pc200023/nobackup/users/thompson/ERA5/tw/ERA5_tw_day_'+str(yr)+'*.nc' # numbers indicate months extracted
        yr_cubes = iris.load(files)
        rm_att = equalise_attributes(yr_cubes)
        unify_time_units(yr_cubes)
        new_cubes = iris.cube.CubeList()
        for each in yr_cubes:
            new_cubes.append(iris.util.squeeze(each))
            if len(new_cubes)>1 and new_cubes[-1].metadata!=new_cubes[0].metadata:
                new_cubes[-1].metadata = new_cubes[0].metadata
        yr_cube = new_cubes.concatenate_cube()
        if len(all_yr_list)>1 and all_yr_list[-1].metadata!=all_yr_list[0].metadata:
            all_yr_list[-1].metadata = all_yr_list[0].metadata
        all_yr_list.append(yr_cube)
    return all_yr_list


# %%
psi_event_tw = gdata.extract_region(gdata.extract_JJA(era5_tw(range(2021,2023))), region)
event_tw = gdata.pull_out_day_era(psi_event_tw, date[0], date[1], date[2])
psi_past_tw = gdata.extract_region(gdata.extract_JJA(era5_tw(range(past_Y1,past_Y2))), region)
psi_present_tw= gdata.extract_region(gdata.extract_JJA(era5_tw(range(present_Y1,present_Y2))), region)

# Composites of streamfunction for region
comp_past = gdata.composite_dates(psi_past_tw, dates_past) 
comp_present = gdata.composite_dates(psi_present_tw, dates_present)

# %%

comp_era_sig = gdata.composite_dates_change_significance(psi_past_tw, dates_past, psi_present_tw, dates_present)

# %%
## COMPOSITES FIGURE
fig, axs = plt.subplots(nrows=2, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,4))
lats=event_tw.coord('latitude').points
lons=event_tw.coord('longitude').points
con_lev = np.linspace(5, 45, 20)
c = axs[0,0].contourf(lons, lats, event_tw.data, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')


c = axs[1,0].contourf(lons, lats, comp_past.data, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,0].set_title('ERA-5 past', loc='right')
c = axs[1,1].contourf(lons, lats, comp_present.data, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,1].set_title('ERA-5 present', loc='right')

#anom = (comp_present.data - comp_past.data)
#levs = np.linspace(-abs(np.min(anom.data)), abs(np.min(anom.data)), 20)
#c2 = axs[1,2].contourf(lons, lats, anom, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
#axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
#axs[1,2].set_title('ERA-5 change', loc='right')


stip_levs = [-2, 0, 2]
anom = (comp_present.data / comp_past.data)*100 -100
c2 = axs[1,2].contourf(lons, lats, anom, levels=np.linspace(-12, 12,20), cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
c3 = axs[1,2].contourf(lons, lats, comp_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,2].set_title('ERA-5 change', loc='right')

## colorbars
fig.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig.add_axes([0.27, 0.1, 0.25, 0.03])
fig.colorbar(c, cax=cbar_ax, ticks=[10, 20, 30, 40], orientation='horizontal')
cbar_ax.set_xlabel('Total column water vapour (kg/m2)')

cbar_ax2 = fig.add_axes([0.66, 0.1, 0.25, 0.03])
fig.colorbar(c2, cax=cbar_ax2, ticks=[-10, 0 ,10], orientation='horizontal')
cbar_ax2.set_xlabel('Absolute Change')

axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[1,1].set_title('(c)', loc='left')
axs[1,2].set_title('(d)', loc='left')

np.mean((comp_present.data / comp_past.data)*100 -100)


