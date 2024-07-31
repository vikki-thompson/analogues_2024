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
## FOR LENTIS
psi_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2013,2023))), R1)
psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R1)
psi_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_mydata(run='F'))), R1)
psi_lentis = gdata.regrid(psi_lentis, event)
psi_lentis_future = gdata.regrid(psi_lentis_future, event)

LEN_P_date_list = ['10472008803', '10572002823', '10472008804', '10692005726', '10692005725', '11132002725', '11132002726', '10572002822', '10542000715', '10962000807', '11122001719', '11332002806', '10342005819', '11242000719', '11102000725', '10832001830', '10732001808', '10572007803', '10982001721', '10692001811', '10752009817', '10692001810', '10842008817', '10752000831', '10862008720', '10192006816', '10872009817', '10742009816', '11122001720', '10472008802']
LEN_F_date_list = ['50292084824', '51452083727', '50122076823', '50362083827', '51522082815', '50412075721', '50222077724', '51332081802', '51362083808', '50732077831', '50632080803', '50372079821', '50142077804', '50512083823', '50752078804', '50322083805', '51192081816', '51342078822', '50722081811', '51102077814', '51692077731', '50792077726', '50262084825', '51242075728', '51202078805', '51052084815', '51482078820', '51372079714', '50932080803', '51432082822']

compD_1 = gdata.composite_dates_lentis(psi_lentis, LEN_P_date_list)
compD_2 = gdata.composite_dates_lentis(psi_lentis_future, LEN_F_date_list)

# %%
## COMPOSITES FIGURE
fig, axs = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,4))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-80, -10, 20)
c = axs[0,0].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')

c = axs[1,0].contourf(lons, lats, comp_past.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,0].set_title('ERA-5 past', loc='right')
c = axs[1,1].contourf(lons, lats, comp_present.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,1].set_title('ERA-5 present', loc='right')

#anom = comp_present.data - comp_past.data
# levs = np.linspace(-abs(np.min(anom.data)), abs(np.min(anom.data)), 20)
anom = (comp_present.data / comp_past.data) *100 # percentage change
levs = np.linspace(85, 115, 40)
c = axs[1,2].contourf(lons, lats, anom, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,2].set_title('ERA-5 change', loc='right')

c = axs[2,0].contourf(lons, lats, compD_1.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,0].set_title('LENTIS present', loc='right')
c = axs[2,1].contourf(lons, lats, compD_2.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,1].set_title('LENTIS future', loc='right')

#anom2 = compD_2.data - compD_1.data
anom2 = (compD_2.data / compD_1.data) *100 # percentage change
c2 = axs[2,2].contourf(lons, lats, anom2, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[2,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,2].set_title('LENTIS change', loc='right')

## colorbars
fig.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig.add_axes([0.27, 0.1, 0.25, 0.03])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
#cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

cbar_ax2 = fig.add_axes([0.66, 0.1, 0.25, 0.03])
fig.colorbar(c2, cax=cbar_ax2, ticks=[90, 100, 110], orientation='horizontal')
cbar_ax2.set_xlabel('Percentage Change')
cbar_ax2.set_xticklabels(['-10', '0', '10'])

axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[1,1].set_title('(c)', loc='left')
axs[1,2].set_title('(d)', loc='left')
axs[2,0].set_title('(e)', loc='left')
axs[2,1].set_title('(f)', loc='left')
axs[2,2].set_title('(g)', loc='left')



# %%
## COMPOSITES FIGURE - ABSOLUTE CHANGE NOT PERCENTAGE
fig, axs = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12,7))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-80, -10, 20)
c = axs[0,0].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')

c = axs[1,0].contourf(lons, lats, comp_past.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,0].set_title('ERA-5 past', loc='right')
c = axs[1,1].contourf(lons, lats, comp_present.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,1].set_title('ERA-5 present', loc='right')

anom = (comp_present.data - comp_past.data)/1000000
levs = np.linspace(-abs(np.min(anom.data)), abs(np.min(anom.data)), 20)
c = axs[1,2].contourf(lons, lats, anom, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,2].set_title('ERA-5 change', loc='right')

c = axs[2,0].contourf(lons, lats, compD_1.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,0].set_title('LENTIS present', loc='right')
c = axs[2,1].contourf(lons, lats, compD_2.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,1].set_title('LENTIS future', loc='right')

anom2 = (compD_2.data - compD_1.data)/1000000
c2 = axs[2,2].contourf(lons, lats, anom2, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[2,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,2].set_title('LENTIS change', loc='right')

## colorbars
fig.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig.add_axes([0.3, 0.1, 0.2, 0.03])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
#cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

cbar_ax2 = fig.add_axes([0.65, 0.1, 0.2, 0.03])
fig.colorbar(c2, cax=cbar_ax2, ticks=[-3, 0, 3], orientation='horizontal')
cbar_ax2.set_xlabel('Absolute Change  (x$10^6$ m$^2$/s)')
cbar_ax2.set_xticklabels(['-3', '0', '3'])

axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[1,1].set_title('(c)', loc='left')
axs[1,2].set_title('(d)', loc='left')
axs[2,0].set_title('(e)', loc='left')
axs[2,1].set_title('(f)', loc='left')
axs[2,2].set_title('(g)', loc='left')

from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

for each in [axs[0,0], axs[1,0], axs[1,1], axs[1,2], axs[2,0], axs[2,1], axs[2,2]]:
    each.set_xticks([-20, 0, 20], crs=ccrs.PlateCarree())
    each.set_yticks([30, 50, 70], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    each.xaxis.set_major_formatter(lon_formatter)
    each.yaxis.set_major_formatter(lat_formatter)

plt.subplots_adjust(top = 1, bottom=0.2, hspace=0.5, wspace=-.2)

# %%
## COMPOSITES FIGURE WITH SIGNIFICANCE

# Composites of streamfunction for region
comp_past = gdata.composite_dates(psi_past, dates_past) 
comp_present = gdata.composite_dates(psi_present, dates_present)
comp_era_sig = gdata.composite_dates_change_significance(psi_past, dates_past, psi_present, dates_present)

fig, axs = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12,7))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-80, -10, 20)
stip_levs = [-2, 0, 2]

c = axs[0,0].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')

c = axs[1,0].contourf(lons, lats, comp_past.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,0].set_title('ERA-5 past', loc='right')
c = axs[1,1].contourf(lons, lats, comp_present.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[1,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,1].set_title('ERA-5 present', loc='right')



compD_1 = gdata.composite_dates_lentis(psi_lentis, LEN_P_date_list)
compD_2 = gdata.composite_dates_lentis(psi_lentis_future, LEN_F_date_list)
comp_len_sig = gdata.composite_dates_change_significance_lentis(psi_lentis, LEN_P_date_list, psi_lentis_future, LEN_F_date_list)

anom = (comp_present.data - comp_past.data)/1000000 # absolute change
levs = np.linspace(-abs(np.min(anom.data)), abs(np.min(anom.data)), 20)

c = axs[1,2].contourf(lons, lats, anom, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
c3 = axs[1,2].contourf(lons, lats, comp_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,2].set_title('ERA-5 change', loc='right')

c = axs[2,0].contourf(lons, lats, compD_1.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,0].set_title('LENTIS present', loc='right')
c = axs[2,1].contourf(lons, lats, compD_2.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[2,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,1].set_title('LENTIS future', loc='right')

anom2 = (compD_2.data - compD_1.data)/1000000
c2 = axs[2,2].contourf(lons, lats, anom2, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
c3 = axs[2,2].contourf(lons, lats, comp_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,2].set_title('LENTIS change', loc='right')

## colorbars
fig.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig.add_axes([0.3, 0.1, 0.2, 0.03])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
#cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

cbar_ax2 = fig.add_axes([0.65, 0.1, 0.2, 0.03])
fig.colorbar(c2, cax=cbar_ax2, ticks=[-3, 0, 3], orientation='horizontal')
cbar_ax2.set_xlabel('Absolute Change  (x$10^6$ m$^2$/s)')
cbar_ax2.set_xticklabels(['-3', '0', '3'])

axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[1,1].set_title('(c)', loc='left')
axs[1,2].set_title('(d)', loc='left')
axs[2,0].set_title('(e)', loc='left')
axs[2,1].set_title('(f)', loc='left')
axs[2,2].set_title('(g)', loc='left')

from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

for each in [axs[0,0], axs[1,0], axs[1,1], axs[1,2], axs[2,0], axs[2,1], axs[2,2]]:
    each.set_xticks([-20, 0, 20], crs=ccrs.PlateCarree())
    each.set_yticks([30, 50, 70], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    each.xaxis.set_major_formatter(lon_formatter)
    each.yaxis.set_major_formatter(lat_formatter)

plt.subplots_adjust(top = 1, bottom=0.2, hspace=0.5, wspace=-.2)



