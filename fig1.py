# %%
# Case Study: Maps of climatic variables of event
#
#
# Original: vikki.thompson 10/05/2023
# Last Editted 10 May 2023

# %%
### Load neccessary libraries
import subprocess
import numpy as np
import sys
sys.path.append('/usr/people/thompson/WP1')
import functions_get_data as gdata
import functions_plot_data as pdata
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

# %%
### Define Variables
sel_year = 2021
sel_month = 'Jul'
sel_day = 14
R1 = [60, 40, 20, -10] # region
NH = [90, 0, 180, -180]
date = [sel_year, sel_month, sel_day]

# %%
## Figure
R1 = [70, 30, 30, -30] # analog region
date = (sel_year, sel_month, sel_day)
event_psi250 = gdata.var_event_data('psi250', R1, date)
event_prec = gdata.var_event_data('tp', R1, date)
fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,4))
lats=event_prec.coord('latitude').points
lons=event_prec.coord('longitude').points
c = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 100, 11), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='max')
cbar = fig.colorbar(c, ticks=np.arange(0, 110, 10))
cbar.ax.set_ylabel('Total Precipitation (mm)', rotation=270, fontsize=12)
cbar.ax.set_yticklabels(['0', '', '20','','40','','60','','80','','100'])

lats=event_psi250.coord('latitude').points
lons=event_psi250.coord('longitude').points
c2 = axs.contour(lons, lats, event_psi250.data/1000000, levels=np.arange(-80, 0, 5), cmap = plt.cm.get_cmap('viridis'), transform=ccrs.PlateCarree(), extend='both')
axs.clabel(c2, inline=1, fontsize=12)
axs.add_feature(cf.BORDERS)
axs.add_feature(cf.COASTLINE)

axs.set_xticks([-20, 0, 20], crs=ccrs.PlateCarree())
axs.set_yticks([30, 50, 70], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
axs.xaxis.set_major_formatter(lon_formatter)
axs.yaxis.set_major_formatter(lat_formatter)


