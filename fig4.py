# %%
# Figures for analog paper, ERA5
# Large scale streamfunction composites
#
# Original: vikki.thompson 26/07/2023
# Last Editted 23 May 2024 (revisions v1)

# %%
import subprocess
import numpy as np
import iris
import sys
import matplotlib.pyplot as plt
sys.path.append('/usr/people/thompson/WP1')
import functions_get_data as gdata
import functions_plot_data as pdata
import cartopy.crs as ccrs
import cartopy as cart
import cartopy.feature as cf
import glob
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
plt.ion(); plt.show()

# %%
## Variables
R1 = [70, 30, 30, -30] # analog region
R2 = [62, 42, 20, -10] # Hylke's region
R3 = [80, 20, 180, -180]
date = [2021, 'Jul', 14] # event date

## Analogues
region = R1
past_Y1 = 1950
past_Y2 = 1980
present_Y1 = 1993
present_Y2 = 2023

LEN_P_date_list = ['10472008803', '10572002823', '10472008804', '10692005726', '10692005725', '11132002725', '11132002726', '10572002822', '10542000715', '10962000807', '11122001719', '11332002806', '10342005819', '11242000719', '11102000725', '10832001830', '10732001808', '10572007803', '10982001721', '10692001811', '10752009817', '10692001810', '10842008817', '10752000831', '10862008720', '10192006816', '10872009817', '10742009816', '11122001720', '10472008802']
LEN_F_date_list = ['50292084824', '51452083727', '50122076823', '50362083827', '51522082815', '50412075721', '50222077724', '51332081802', '51362083808', '50732077831', '50632080803', '50372079821', '50142077804', '50512083823', '50752078804', '50322083805', '51192081816', '51342078822', '50722081811', '51102077814', '51692077731', '50792077726', '50262084825', '51242075728', '51202078805', '51052084815', '51482078820', '51372079714', '50932080803', '51432082822']

# %%
## Dates lists for composites, based on European domain (R1)
psi_event = gdata.pull_out_day_era(gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2021,2023))), R1), date[0], date[1], date[2])
psi_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), R1)
psi_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), R1)
dates_past = gdata.analogs_datelist(psi_past, psi_event)[:30] 
dates_present = gdata.analogs_datelist(psi_present, psi_event)[:30]

# Composites of streamfunction for Northern Hemisphere
## EVENT: removing the zonal mean
psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(1950,2023))), R3)
psi_zonal = psi.collapsed(('time'), iris.analysis.MEAN).collapsed(('longitude'), iris.analysis.MEAN)
event = gdata.pull_out_day_era(psi, date[0], date[1], date[2])
event = event - psi_zonal

## ERA5: composites of top 30 analogues, removing zonal mean
psi_past_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), R3)
psi_present_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), R3)
psi_present_era_zonal = psi_present_era.collapsed(('longitude', 'time'), iris.analysis.MEAN)
#psi_present_era_zonal = psi_present_era.collapsed(('time'), iris.analysis.MEAN)
psi_present_era = psi_present_era - psi_present_era_zonal
psi_past_era_zonal = psi_past_era.collapsed(('longitude', 'time'), iris.analysis.MEAN)
#psi_past_era_zonal = psi_past_era.collapsed(('time'), iris.analysis.MEAN)
psi_past_era = psi_past_era - psi_past_era_zonal
comp_present_era = gdata.composite_dates(psi_present_era, dates_present)
comp_present_era = comp_present_era - np.mean(comp_present_era.data)
comp_present_era_sig = gdata.composite_dates_significance(psi_present_era, dates_present)

comp_past_era = gdata.composite_dates(psi_past_era, dates_past)
comp_past_era = comp_past_era - np.mean(comp_past_era.data)
comp_past_era_sig = gdata.composite_dates_significance(psi_past_era, dates_past)

## LENTIS: composites of top 30 analogues, removing zonal mean
psi_present_len = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R3)
psi_future_len = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata(run='F')), R3)
psi_present_len = gdata.regrid(psi_present_len, event)
psi_future_len = gdata.regrid(psi_future_len, event)
psi_present_len_zonal = psi_present_len[0,...].collapsed(('realization', 'longitude', 'time'), iris.analysis.MEAN)
#psi_present_len_zonal = psi_present_len[0,...].collapsed(('realization', 'time'), iris.analysis.MEAN)
psi_present_len = psi_present_len - psi_present_len_zonal
psi_future_len_zonal = psi_future_len[0,...].collapsed(('realization', 'longitude', 'time'), iris.analysis.MEAN)
#psi_future_len_zonal = psi_future_len[0,...].collapsed(('realization', 'time'), iris.analysis.MEAN)
psi_future_len = psi_future_len - psi_future_len_zonal
comp_present_len = gdata.composite_dates_lentis(psi_present_len, LEN_P_date_list)
comp_future_len = gdata.composite_dates_lentis(psi_future_len, LEN_F_date_list)
comp_present_len = comp_present_len - np.mean(comp_present_len.data)
comp_future_len = comp_future_len - np.mean(comp_future_len.data)
comp_future_len_sig = gdata.composite_dates_significance_lentis(psi_future_len, LEN_F_date_list)
comp_present_len_sig = gdata.composite_dates_significance_lentis(psi_present_len, LEN_P_date_list)

# %%
# statistical testing: one sample t-test
# H0: mean is zero
import calendar
import scipy.stats as stats

def composite_dates_ttest(psi, date_list):
    '''
    Returns single composite of all dates
    Inputs required:
      psi = list of cubes, 1 per year - as used to calc D/date_list
      date_list = list of events to composite
    '''
    n = len(date_list)
    field_list = iris.cube.CubeList([])
    for each in range(n):
        year = int(date_list[each][:4])
        month = calendar.month_abbr[int(date_list[each][4:-2])]
        day = int(date_list[each][-2:])
        field_list.append(gdata.pull_out_day_era(psi, year, month, day))
    sig_field = field_list[0].data
    a, b = np.shape(field_list[0].data)
    for i in range(a):
        print(i)
        for j in range(b):
            loc_list = []
            for R in range(n):
                loc_list.append(field_list[R].data[i,j])
            t_stat, p_val = stats.ttest_1samp(loc_list, 0)
            if p_val < 0.05:
                sig_field[i,j] = 1
            else:
                sig_field[i,j] = 0
    result_cube = field_list[0]
    result_cube.data = sig_field
    return result_cube


def composite_dates_ttest_lentis(psi, date_list):
    '''
    Returns single composite of all dates
    Inputs required:
      psi = list of cubes, 1 per year - as used to calc D/date_list
      date_list = list of events to composite
    '''
    n = len(date_list)
    field_list = iris.cube.CubeList([])
    for each in range(n):
        real = int(date_list[each][:4])
        year = int(date_list[each][4:8])
        month = calendar.month_abbr[int(date_list[each][8])]
        day = int(date_list[each][-2:])
        field_list.append(gdata.pull_out_day_lentis(psi, real, year, month, day))
    sig_field = field_list[0].data
    a, b = np.shape(field_list[0].data)
    for i in range(a):
        print(i)
        for j in range(b):
            loc_list = []
            for R in range(n):
                loc_list.append(field_list[R].data[i,j])
            t_stat, p_val = stats.ttest_1samp(loc_list, 0)
            if p_val < 0.05:
                sig_field[i,j] = 1
            else:
                sig_field[i,j] = 0
    result_cube = field_list[0]
    result_cube.data = sig_field
    return result_cube


# %%
# redo significance using ttest
comp_present_era_sig = composite_dates_ttest(psi_present_era, dates_present)
comp_past_era_sig = composite_dates_ttest(psi_past_era, dates_past)

comp_future_len_sig = composite_dates_ttest_lentis(psi_future_len, LEN_F_date_list)
comp_present_len_sig = composite_dates_ttest_lentis(psi_present_len, LEN_P_date_list)

# %%
## COMPOSITES FIGURE
fig, axs = plt.subplots(nrows=3, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-30, 30, 20)
stip_levs = [-2, 0, 2]

c = axs[0,0].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[0,0].set_title('Event')
pdata.plot_box(axs[0,0], R1)

c = axs[0,1].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')

fig.subplots_adjust(bottom=-.23, hspace=-0.8)
cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-20, -10, 0, 10, 20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
#cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

c = axs[1,0].contourf(lons, lats, comp_past_era.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[1,0].contourf(lons, lats, comp_past_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[1,0].set_title('ERA-5 past')
pdata.plot_box(axs[1,0], R1)

c = axs[2,0].contourf(lons, lats, comp_present_era.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[2,0].contourf(lons, lats, comp_present_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[2,0].set_title('ERA-5 present')
pdata.plot_box(axs[2,0], R1)

c = axs[1,1].contourf(lons, lats, comp_present_len.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[1,1].contourf(lons, lats, comp_present_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,1].add_feature(cf.COASTLINE, linewidth=.5)
axs[1,1].set_title('LENTIS present')
pdata.plot_box(axs[1,1], R1)

c = axs[2,1].contourf(lons, lats, comp_future_len.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[2,1].contourf(lons, lats, comp_future_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,1].add_feature(cf.COASTLINE, linewidth=.5)
axs[2,1].set_title('LENTIS future')
pdata.plot_box(axs[2,1], R1)


axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[2,0].set_title('(c)', loc='left')
axs[1,1].set_title('(d)', loc='left')
axs[2,1].set_title('(e)', loc='left')


for each in [axs[0,0], axs[1,0], axs[1,1], axs[2,0], axs[2,1]]:
    each.set_xticks([-120, -60, 0, 60, 120], crs=ccrs.PlateCarree())
    each.set_yticks([20, 40, 60, 80], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    each.xaxis.set_major_formatter(lon_formatter)
    each.yaxis.set_major_formatter(lat_formatter)

plt.tight_layout()
plt.subplots_adjust(top = 1, bottom=0.15, hspace=0.1)


# %%

        

## Plot Climatological Mean / Anomaly

## EVENT: removing the zonal mean
psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(1950,2023))), R3)
psi_zonal = psi.collapsed(('time'), iris.analysis.MEAN)
event = gdata.pull_out_day_era(psi, date[0], date[1], date[2])
event = event - psi_zonal

## ERA5: composites of top 30 analogues, removing zonal mean
psi_past_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), R3)
psi_present_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), R3)
psi_present_era_zonal = psi_present_era.collapsed(( 'time'), iris.analysis.MEAN)
psi_present_era = psi_present_era - psi_present_era_zonal
psi_past_era_zonal = psi_past_era.collapsed(('time'), iris.analysis.MEAN)
psi_past_era = psi_past_era - psi_past_era_zonal
comp_present_era = gdata.composite_dates(psi_present_era, dates_present)
comp_present_era = comp_present_era - np.mean(comp_present_era.data)
comp_present_era_sig = gdata.composite_dates_significance(psi_present_era, dates_present)

comp_past_era = gdata.composite_dates(psi_past_era, dates_past)
comp_past_era = comp_past_era - np.mean(comp_past_era.data)
comp_past_era_sig = gdata.composite_dates_significance(psi_past_era, dates_past)

## LENTIS: composites of top 30 analogues, removing zonal mean
psi_present_len = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R3)
psi_future_len = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata(run='F')), R3)
psi_present_len = gdata.regrid(psi_present_len, event)
psi_future_len = gdata.regrid(psi_future_len, event)
psi_present_len_zonal = psi_present_len[0,...].collapsed(('realization', 'time'), iris.analysis.MEAN)
psi_present_len = psi_present_len - psi_present_len_zonal
psi_future_len_zonal = psi_future_len[0,...].collapsed(('realization', 'time'), iris.analysis.MEAN)
psi_future_len = psi_future_len - psi_future_len_zonal
comp_present_len = gdata.composite_dates_lentis(psi_present_len, LEN_P_date_list)
comp_future_len = gdata.composite_dates_lentis(psi_future_len, LEN_F_date_list)
comp_present_len = comp_present_len - np.mean(comp_present_len.data)
comp_future_len = comp_future_len - np.mean(comp_future_len.data)
comp_future_len_sig = gdata.composite_dates_significance_lentis(psi_future_len, LEN_F_date_list)
comp_present_len_sig = gdata.composite_dates_significance_lentis(psi_present_len, LEN_P_date_list)


## COMPOSITES FIGURE
fig, axs = plt.subplots(nrows=3, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-30, 30, 20)
stip_levs = [-2, 0, 2]

c = axs[0,0].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[0,0].set_title('Event')
pdata.plot_box(axs[0,0], R1)

c = axs[0,1].contourf(lons, lats, event.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')

fig.subplots_adjust(bottom=-.23, hspace=-0.8)
cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-20, -10, 0, 10, 20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
#cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

c = axs[1,0].contourf(lons, lats, comp_past_era.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[1,0].contourf(lons, lats, comp_past_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[1,0].set_title('ERA-5 past')
pdata.plot_box(axs[1,0], R1)

c = axs[2,0].contourf(lons, lats, comp_present_era.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[2,0].contourf(lons, lats, comp_present_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,0].add_feature(cf.COASTLINE, linewidth=.5)
axs[2,0].set_title('ERA-5 present')
pdata.plot_box(axs[2,0], R1)

'''
#anom = comp_present_era.data/1000000 - comp_past_era.data/1000000
anom = (comp_present_era.data / comp_past_era.data) *100 # percentage change
levs = np.linspace(50, 150, 40)
c = axs[3,0].contourf(lons, lats, anom, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[3,0].add_feature(cf.BORDERS)
axs[3,0].add_feature(cf.COASTLINE)
axs[3,0].set_title('ERA-5 change')
pdata.plot_box(axs[3,0], R1)
'''

c = axs[1,1].contourf(lons, lats, comp_present_len.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[1,1].contourf(lons, lats, comp_present_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,1].add_feature(cf.COASTLINE, linewidth=.5)
axs[1,1].set_title('LENTIS present')
pdata.plot_box(axs[1,1], R1)

c = axs[2,1].contourf(lons, lats, comp_future_len.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
c2 = axs[2,1].contourf(lons, lats, comp_future_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,1].add_feature(cf.COASTLINE, linewidth=.5)
axs[2,1].set_title('LENTIS future')
pdata.plot_box(axs[2,1], R1)

'''
#anom2 = comp_future_len.data - comp_present_len.data
anom2 = (comp_future_len.data / comp_present_len.data) *100 # percentage change
c = axs[3,1].contourf(lons, lats, anom2, levels=levs, cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
axs[3,1].add_feature(cf.BORDERS)
axs[3,1].add_feature(cf.COASTLINE)
axs[3,1].set_title('LENTIS change')
pdata.plot_box(axs[3,1], R1)

cbar_ax = fig.add_axes([0.05, 0.1, 0.01, 0.2])
fig.colorbar(c, cax=cbar_ax, ticks=[60, 100, 140])
cbar_ax.set_ylabel('Percentage change', rotation=270, labelpad=15)
'''


axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[2,0].set_title('(c)', loc='left')
axs[1,1].set_title('(d)', loc='left')
axs[2,1].set_title('(e)', loc='left')

'''
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
for each in [axs[0,0], axs[1,0], axs[1,1], axs[2,0], axs[2,1]]:
        gl = each.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.75)
        gl.xlabels_top = False; gl.xlabels_bottom = False
        gl.ylabels_left = False; gl.ylabels_right= False
        gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
        gl.ylocator = mticker.FixedLocator([20, 40, 60, 80])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

'''

        





psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(1950,2023))), R3)
psi_zonal = psi.collapsed(('time'), iris.analysis.MEAN)

fig, axs = plt.subplots(nrows=3, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
con_lev = np.linspace(-80, 30, 20)

c = axs[0,0].contourf(lons, lats, psi_zonal.data/1000000, levels=con_lev, cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=.5)
pdata.plot_box(axs[0,0], R1)


## ERA5: composites of top 30 analogues, climatological anomaly??
psi_past_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), R3)
psi_present_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), R3)
#psi_present_era_zonal = psi_present_era.collapsed(('longitude', 'time'), iris.analysis.MEAN)
psi_present_era_clim = psi_present_era.collapsed(('time'), iris.analysis.MEAN)
psi_present_era = psi_present_era - psi_present_era_clim
#psi_past_era_zonal = psi_past_era.collapsed(('longitude', 'time'), iris.analysis.MEAN)
psi_past_era_clim = psi_past_era.collapsed(('time'), iris.analysis.MEAN)
psi_past_era = psi_past_era - psi_past_era_clim
comp_present_era = gdata.composite_dates(psi_present_era, dates_present)
comp_past_era = gdata.composite_dates(psi_past_era, dates_past)



import calendar
fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
event_max = np.max(event.data)
# subplots
all_ax = axs.ravel()
for i, each in enumerate(dates_past):
        YY = np.int(each[:4])
        MM = calendar.month_abbr[int(each[4])]
        DD = int(each[-2:])
        event = gdata.pull_out_day_era(psi_past_era, YY, MM, DD)
        c = all_ax[i].contourf(lons, lats, event.data, levels=np.linspace(-event_max, event_max, 21), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
        pdata.plot_box(all_ax[i], R1) # plot box of analog region
        all_ax[i].add_feature(cf.BORDERS)
        all_ax[i].add_feature(cf.COASTLINE)
        all_ax[i].set_title(each)
        
fig.suptitle(title)
plt.tight_layout() 



