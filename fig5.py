# %%
# Figures for analog paper, ERA5
#
#
# Original: vikki.thompson 26/07/2023
# Last Editted 23 May 2024 (rev v1)

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
import calendar
import cartopy.crs as ccrs
import cartopy as cart
import cartopy.feature as cf
import scipy.stats as stats
import glob
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

# %%
# Get data, find analogues based on streamfunction

## For ERA-5
psi_event = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2021,2023))), region)
event = gdata.pull_out_day_era(psi_event, date[0], date[1], date[2])
psi_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), region)
psi_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), region)
dates_past = gdata.analogs_datelist(psi_past, event)[:30] # find analogs
dates_present = gdata.analogs_datelist(psi_present, event)[:30]

## FOR LENTIS
psi_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2013,2023))), R1)
psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R1)
psi_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_mydata(run='F'))), R1)
psi_lentis = gdata.regrid(psi_lentis, event)
psi_lentis_future = gdata.regrid(psi_lentis_future, event)
LEN_P_date_list = ['10472008803', '10572002823', '10472008804', '10692005726', '10692005725', '11132002725', '11132002726', '10572002822', '10542000715', '10962000807', '11122001719', '11332002806', '10342005819', '11242000719', '11102000725', '10832001830', '10732001808', '10572007803', '10982001721', '10692001811', '10752009817', '10692001810', '10842008817', '10752000831', '10862008720', '10192006816', '10872009817', '10742009816', '11122001720', '10472008802']
LEN_F_date_list = ['50292084824', '51452083727', '50122076823', '50362083827', '51522082815', '50412075721', '50222077724', '51332081802', '51362083808', '50732077831', '50632080803', '50372079821', '50142077804', '50512083823', '50752078804', '50322083805', '51192081816', '51342078822', '50722081811', '51102077814', '51692077731', '50792077726', '50262084825', '51242075728', '51202078805', '51052084815', '51482078820', '51372079714', '50932080803', '51432082822']

# %%
import calendar

def composite_max_dates(psi, date_list):
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
        month = calendar.month_abbr[int(date_list[each][4])]
        day = int(date_list[each][-2:])
        field_list.append(gdata.pull_out_day_era(psi, year, month, day))
    max_field = field_list[0].data
    a, b = np.shape(field_list[0].data)
    for i in range(a):
        print(i)
        for j in range(b):
            loc_list = []
            for R in range(n):
                loc_list.append(field_list[R].data[i,j])
            max_field[i,j] = np.max(loc_list)
    result_cube = field_list[0]
    result_cube.data = max_field
    return result_cube

# %%
# Composities of rainfall for streamfunction analogues - maximum gridpoint value

R2 = [62, 35, 25, -10] 
# dates_past and dates_present
pr_event = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(2021,2023))), R2)
event = gdata.pull_out_day_era(pr_event, date[0], date[1], date[2])
pr_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(past_Y1,past_Y2))), R2)
pr_present = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(present_Y1,present_Y2))), R2)
#comp_past = gdata.composite_dates(pr_past, dates_past)
#comp_present = gdata.composite_dates(pr_present, dates_present)
max_past = composite_max_dates(pr_past, dates_past)
max_present = composite_max_dates(pr_present, dates_present)

## LEN
pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), R2)
pr_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_data(run='2K'))), R2)
pr_lentis = gdata.regrid(pr_lentis, event)
pr_lentis_future = gdata.regrid(pr_lentis_future, event)

# %%
def max_field_comp(field_cube, date_list):
    n = len(date_list)
    a, b = np.shape(field_cube)[2:4]
    field = np.empty([n, a, b])
    for each in range(n):
        #print(date_list[each])
        real = int(date_list[each][:4])
        year = int(date_list[each][4:8])
        month = calendar.month_abbr[int(date_list[each][8])]
        day = int(date_list[each][-2:])
        field[each,...] = gdata.pull_out_day_lentis(field_cube, real, year, month, day).data
    return np.max(field, axis=0) * 86400.


# %%
max_present_lentis = max_field_comp(pr_lentis, LEN_P_date_list) 
max_future_lentis = max_field_comp(pr_lentis_future, LEN_F_date_list)

comp_era_sig = gdata.composite_dates_change_significance(pr_past, dates_past, pr_present, dates_present)
comp_len_sig = gdata.composite_dates_change_significance_lentis(pr_lentis, LEN_P_date_list, pr_lentis_future, LEN_F_date_list)

# %%
## COMPOSITES FIGURE
fig1, axs = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12,7))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
levs = np.linspace(0, 80, 10)
stip_levs = [-2, 0, 2]
c = axs[0,0].contourf(lons, lats, event.data, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')

c = axs[1,0].contourf(lons, lats, max_past.data, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[1,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,0].set_title('ERA-5 past', loc='right')
c2 = axs[1,1].contourf(lons, lats, max_present.data, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[1,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,1].set_title('ERA-5 present', loc='right')

anom = max_present.data - max_past.data
c = axs[1,2].contourf(lons, lats, anom, levels=np.linspace(-40, 40, 10), cmap = plt.cm.get_cmap('BrBG'), transform=ccrs.PlateCarree(), extend='both')
c3 = axs[1,2].contourf(lons, lats, comp_era_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[1,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[1,2].set_title('ERA-5 change', loc='right')

c = axs[2,0].contourf(lons, lats, max_present_lentis, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[2,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,0].set_title('LENTIS present', loc='right')
c = axs[2,1].contourf(lons, lats, max_future_lentis, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[2,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,1].set_title('LENTIS future', loc='right')

anom2 = max_future_lentis - max_present_lentis
c2 = axs[2,2].contourf(lons, lats, anom2, levels=np.linspace(-40, 40, 10), cmap = plt.cm.get_cmap('BrBG'), transform=ccrs.PlateCarree(), extend='both')
c3 = axs[2,2].contourf(lons, lats, comp_len_sig.data, levels=stip_levs, hatches=[None, '////'], colors='none', transform=ccrs.PlateCarree())
axs[2,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,2].set_title('LENTIS change', loc='right')


## colorbars
fig1.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig1.add_axes([0.23, 0.1, 0.25, 0.02])
fig1.colorbar(c, cax=cbar_ax, ticks=[0, 40, 80], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
cbar_ax.set_xticklabels(['0', '40','80'])

cbar_ax2 = fig1.add_axes([0.66, 0.1, 0.25, 0.02])
fig1.colorbar(c2, cax=cbar_ax2, ticks=[-40, 0, 40], orientation='horizontal')
cbar_ax2.set_xlabel('Rainfall Change (mm)')
cbar_ax2.set_xticklabels(['-40', '0', '40'])

axs[0,0].set_title('(a)', loc='left')
axs[1,0].set_title('(b)', loc='left')
axs[1,1].set_title('(c)', loc='left')
axs[1,2].set_title('(d)', loc='left')
axs[2,0].set_title('(e)', loc='left')
axs[2,1].set_title('(f)', loc='left')
axs[2,2].set_title('(g)', loc='left')

reg_west = [54, 48, 10, 4] # S, N, W, E
reg_east = [54, 48, 20, 14]

from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

all_ax = axs.ravel()
for i in range(9):
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])
    pdata.plot_box(all_ax[i], reg_west)
    pdata.plot_box(all_ax[i], reg_east)
    all_ax[i].set_xticks([0, 10, 20], crs=ccrs.PlateCarree())
    all_ax[i].set_yticks([45, 55], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    all_ax[i].xaxis.set_major_formatter(lon_formatter)
    all_ax[i].yaxis.set_major_formatter(lat_formatter)

    
plt.tight_layout()
plt.subplots_adjust(top=1, bottom=0.2, hspace=0.5, wspace=-0.2)


# %%
## LENTIS ANOMALY FIELD
pr_lentis_mean = pr_lentis.collapsed(('realization','time'), iris.analysis.MEAN)
pr_lentis_future_mean = pr_lentis_future.collapsed(('realization','time'), iris.analysis.MEAN)
## COMPOSITES FIGURE
fig1, axs = plt.subplots(nrows=3, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,4))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
levs = np.linspace(0, 80, 10)
c = axs[0,0].contourf(lons, lats, event.data, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[0,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[0,0].set_title('Event', loc='right')

c = axs[2,0].contourf(lons, lats, pr_lentis_mean.data, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[2,0].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,0].set_title('LENTIS present', loc='right')
c = axs[2,1].contourf(lons, lats, max_future_lentis, levels=levs, cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
axs[2,1].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,1].set_title('LENTIS future', loc='right')

anom2 = max_future_lentis - max_present_lentis
c2 = axs[2,2].contourf(lons, lats, anom2, levels=np.linspace(-40, 40, 10), cmap = plt.cm.get_cmap('BrBG'), transform=ccrs.PlateCarree(), extend='both')
axs[2,2].add_feature(cf.COASTLINE, linewidth=0.5)
axs[2,2].set_title('LENTIS change', loc='right')

## colorbars
fig1.subplots_adjust(bottom=0.15, wspace=0.1, hspace=.2)
cbar_ax = fig1.add_axes([0.27, 0.1, 0.25, 0.02])
fig1.colorbar(c, cax=cbar_ax, ticks=[0, 40, 80], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
cbar_ax.set_xticklabels(['0', '', '20','','40','','60','','80','','100'])

cbar_ax2 = fig1.add_axes([0.66, 0.1, 0.25, 0.02])
fig1.colorbar(c2, cax=cbar_ax2, ticks=[-40, 0, 40], orientation='horizontal')
cbar_ax2.set_xlabel('Rainfall Change (mm)')
cbar_ax2.set_xticklabels(['-40', '0', '40'])

reg_west = [54, 48, 10, 4] # S, N, W, E
reg_east = [54, 48, 20, 14]#[54, 48, 16, 10]

all_ax = axs.ravel()
for i in range(9):
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])
    pdata.plot_box(all_ax[i], reg_west)
    pdata.plot_box(all_ax[i], reg_east)


# %%
### For definied region, plot
### 1. Histogram of val from all grid boxes, all analogues
### 2. Histogram of sum over region, each analogue

def gp_and_sum(data, date_list):
    gp = []
    reg_sum = []
    a, b, c = np.shape(data[0])
    for i in range(30):
        print(i)
        real = int(date_list[i][:4])
        year = int(date_list[i][4:8])
        month = calendar.month_abbr[int(date_list[i][8])]
        day = int(date_list[i][-2:])
        field = gdata.pull_out_day_lentis(data, real, year, month, day)*86400
        gp.append(np.reshape(field.data, (b*c,)))
        reg_sum.append(np.sum(field.data))
    return gp, reg_sum

def gp_and_sum_era(data, date_list, len_regrid):
    gp = []
    reg_sum = []
    a, b, c = np.shape(len_regrid)
    for i in range(30):
        print(i)
        year = int(date_list[i][:4])
        month = calendar.month_abbr[int(date_list[i][4])]
        day = int(date_list[i][-2:])
        for each in data:
            each = gdata.regrid(each, len_regrid)
            if len(each.coords('year')) > 0:
                pass
            else:
                iris.coord_categorisation.add_year(each, 'time')
            if each.coord('year').points[0] == year:
                field = gdata.pull_out_day_era(each, year, month, day)
            else:
                print('skip')
        gp.append(np.reshape(field.data, (b*c,)))
        reg_sum.append(np.sum(field.data))
    return gp, reg_sum

# %%
reg_west = [54, 48, 10, 4] # S, N, W, E
reg_east = [54, 48, 20, 14]#[54, 48, 16, 10]

reg = reg_west
pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), reg)
pr_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_data(run='2K'))), reg)
pr_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(past_Y1,past_Y2))), reg)
pr_present = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(present_Y1,present_Y2))), reg) 
gp_P_west, reg_sum_P = gp_and_sum(pr_lentis, LEN_P_date_list)
gp_F_west, reg_sum_F = gp_and_sum(pr_lentis_future, LEN_F_date_list)
gp_P_era_west, reg_sum_P_era = gp_and_sum_era(pr_past, dates_past, pr_lentis[0])
gp_F_era_west, reg_sum_F_era = gp_and_sum_era(pr_present, dates_present, pr_lentis[0])

reg = reg_east
pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), reg)
pr_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_data(run='2K'))), reg)
pr_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(past_Y1,past_Y2))), reg)
pr_present = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(present_Y1,present_Y2))), reg) 
gp_P_east, reg_sum_P = gp_and_sum(pr_lentis, LEN_P_date_list)
gp_F_east, reg_sum_F = gp_and_sum(pr_lentis_future, LEN_F_date_list)
gp_P_era_east, reg_sum_P_era = gp_and_sum_era(pr_past, dates_past, pr_lentis[0])
gp_F_era_east, reg_sum_F_era = gp_and_sum_era(pr_present, dates_present, pr_lentis[0])


# %%

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(5,3))
nbins = np.linspace(0,70,50)
a, b = np.shape(gp_P_era_west)
axs[0,0].hist(np.reshape(gp_P_era_west, (a*b,)), bins = nbins, color='r', alpha = .5, label='Past')
axs[0,0].hist(np.reshape(gp_F_era_west, (a*b,)), bins = nbins, color='c', alpha = .5, label='Present')
axs[0,0].set_yscale('log')
axs[0,0].legend()
a, b = np.shape(gp_P_west)
axs[0,1].hist(np.reshape(gp_P_west, (a*b,)), bins = nbins, color='r', alpha = .5, label='Present')
axs[0,1].hist(np.reshape(gp_F_west, (a*b,)), bins = nbins, color='c', alpha = .5, label='Future')
axs[0,1].set_yscale('log')
axs[0,1].legend()
axs[0,0].set_ylabel('Western Europe')
axs[0,0].set_title('ERA-5', loc='right')
axs[0,1].set_title('KNMI-LENTIS', loc='right')

a, b = np.shape(gp_P_era_east)
axs[1,0].hist(np.reshape(gp_P_era_east, (a*b,)), bins = nbins, color='r', alpha = .5)
axs[1,0].hist(np.reshape(gp_F_era_east, (a*b,)), bins = nbins, color='c', alpha = .5)
axs[1,0].set_yscale('log')
a, b = np.shape(gp_P_east)
axs[1,1].hist(np.reshape(gp_P_east, (a*b,)), bins = nbins, color='r', alpha = .5)
axs[1,1].hist(np.reshape(gp_F_east, (a*b,)), bins = nbins, color='c', alpha = .5)
axs[1,1].set_yscale('log')
axs[1,0].set_ylabel('Eastern Europe')
axs[1,0].set_xlabel('Gridpoint daily rainfall (mm)')
axs[1,1].set_xlabel('Gridpoint daily rainfall (mm)')

axs[0,0].set_title('(h)', loc='left')
axs[1,0].set_title('(j)', loc='left')
axs[0,1].set_title('(i)', loc='left')
axs[1,1].set_title('(k)', loc='left')

all_ax = axs.ravel()
for i in range(4):
    all_ax[i].set_xlim([0, 70])
    all_ax[i].set_yticks([])

    
plt.tight_layout()
plt.subplots_adjust(top=1, bottom=0.2, hspace=0.5)

# %%


# Violin plots of gridpoint rainfall
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(6,7))
a, b = np.shape(gp_P_era)
v1 = ax[0,0].violinplot([np.reshape(gp_P_era, (a*b,)), np.reshape(gp_F_era, (a*b,))], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
ax[0,0].set_yscale('log')
colors = ['g','b']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)
    

#ax[0,0].axhline(np.mean(Q1_ana_v2), color='g')
#ax[0,0].axhline(np.mean(Q2_ana_v2), color='b')
ax[0,0].set_xticks([1, 1.6])
ax[0,0].set_xticklabels(['Past', 'Present'])
ax[0,0].set_ylabel('ERA-5')
ax[0,0].set_title(' Quality, x$10^10$ m$^2$/s')
ax[0,0].set_title('(a)', loc='left')

v2 = ax[1,0].violinplot([Q1_mod_v2, Q2_mod_v2], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v2['bodies'], colors):
    pc.set_facecolor(color)

ax[1,0].plot(1, Q1_mod_event/1e10, marker='o', color='r')
ax[1,0].plot(1.6, Q2_mod_event/1e10, marker='o', color='r')
ax[1,0].axhline(np.mean(Q1_mod_v2), color='b')
ax[1,0].axhline(np.mean(Q2_mod_v2), color='m')
ax[1,0].set_xticks([1, 1.6])
ax[1,0].set_xticklabels(['Present', 'Future'])
ax[1,0].set_ylabel('KNMI-LENTIS')
ax[1,0].set_title('(b)', loc='left')






#fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(5,5))
#axs[0].hist(reg_sum_P_era, bins = np.linspace(0,3000,10), color='r', alpha = .5)
#axs[0].hist(reg_sum_F_era, bins = np.linspace(0,3000,10), color='g', alpha = .5)
#axs[0].set_yscale('log')
#axs[1].hist(reg_sum_P, bins = np.linspace(0,1000,10), color='r', alpha = .5)
#axs[1].hist(reg_sum_F, bins = np.linspace(0,1000,10), color='g', alpha = .5)
#axs[1].set_yscale('log')
#plt.title('Western region, sum')






### All anologues rainfall fields - present LENTIS
fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
import calendar
# subplots
all_ax = axs.ravel()
date_list = LEN_P_date_list
full_field = pr_lentis
for i in range(30):
    print(i)
    real = int(date_list[i][:4])
    year = int(date_list[i][4:8])
    month = calendar.month_abbr[int(date_list[i][8])]
    day = int(date_list[i][-2:])
    field = gdata.pull_out_day_lentis(full_field, real, year, month, day)*86400
    c = all_ax[i].contourf(field.coord('longitude').points, field.coord('latitude').points, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.BORDERS)
    all_ax[i].add_feature(cf.COASTLINE)
    all_ax[i].set_title(date_list[i])



### All anologues rainfall fields - future
fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
lats=event.coord('latitude').points
lons=event.coord('longitude').points
# subplots
all_ax = axs.ravel()
date_list = LEN_F_date_list
full_field = pr_lentis_future
for i in range(30):
    print(i)
    real = int(date_list[i][:4])
    year = int(date_list[i][4:8])
    month = calendar.month_abbr[int(date_list[i][8])]
    day = int(date_list[i][-2:])
    field = gdata.pull_out_day_lentis(full_field, real, year, month, day)*86400
    c = all_ax[i].contourf(lons, lats, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.BORDERS)
    all_ax[i].add_feature(cf.COASTLINE)
    all_ax[i].set_title(date_list[i])



pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), R2)
psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R2)

date = '11132002726'

real = int(date[:4])
year = int(date[4:8])
month = calendar.month_abbr[int(date[8])]
day = int(date[-2:])
event_prec = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
event_psi250 = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)

fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats = pr_lentis.coord('latitude').points
lons = pr_lentis.coord('longitude').points
c = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
lats = psi_lentis.coord('latitude').points
lons = psi_lentis.coord('longitude').points
c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
axs.add_feature(cf.BORDERS)
axs.add_feature(cf.COASTLINE)


date = '10862008720'

real = int(date[:4])
year = int(date[4:8])
month = calendar.month_abbr[int(date[8])]
day = int(date[-2:])
event_prec = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
event_psi250 = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)

fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats = pr_lentis.coord('latitude').points
lons = pr_lentis.coord('longitude').points
c1 = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
lats = psi_lentis.coord('latitude').points
lons = psi_lentis.coord('longitude').points
c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
axs.add_feature(cf.BORDERS)
axs.add_feature(cf.COASTLINE)



date = '10472008804'

real = int(date[:4])
year = int(date[4:8])
month = calendar.month_abbr[int(date[8])]
day = int(date[-2:])
event_prec = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
event_psi250 = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)

fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats = pr_lentis.coord('latitude').points
lons = pr_lentis.coord('longitude').points
c1 = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
lats = psi_lentis.coord('latitude').points
lons = psi_lentis.coord('longitude').points
c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
axs.add_feature(cf.BORDERS)
axs.add_feature(cf.COASTLINE)




date = '10572002823'

real = int(date[:4])
year = int(date[4:8])
month = calendar.month_abbr[int(date[8])]
day = int(date[-2:])
event_prec = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
event_psi250 = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)

fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
lats = pr_lentis.coord('latitude').points
lons = pr_lentis.coord('longitude').points
c1 = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
lats = psi_lentis.coord('latitude').points
lons = psi_lentis.coord('longitude').points
c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
axs.add_feature(cf.BORDERS)
axs.add_feature(cf.COASTLINE)



## Rainfall over region (Benelux?)
R_lim = R1#[52, 48, 7, 5]
pr_max_Epast = []; pr_max_Epres = []
for i in range(30):
    print(i)
    year = int(dates_past[i][:4])
    month = calendar.month_abbr[int(dates_past[i][4])]
    day = int(dates_past[i][-2:])
    field = gdata.pull_out_day_era(pr_past, year, month, day)
    pr_field = gdata.extract_region(field, R_lim)
    pr_max_Epast.append(np.max(pr_field.data))
    year = int(dates_present[i][:4])
    month = calendar.month_abbr[int(dates_present[i][4])]
    day = int(dates_present[i][-2:])
    field = gdata.pull_out_day_era(pr_present, year, month, day)
    pr_field = gdata.extract_region(field, R_lim)
    pr_max_Epres.append(np.max(pr_field.data))

pr_max_LP = []; pr_max_LF = []
for i in range(30):
    print(i)
    real = int(LEN_P_date_list[i][:4])
    year = int(LEN_P_date_list[i][4:8])
    month = calendar.month_abbr[int(LEN_P_date_list[i][8])]
    day = int(LEN_P_date_list[i][-2:])
    field = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
    pr_field = gdata.extract_region(field, R_lim)
    pr_max_LP.append(np.max(pr_field.data))
    real = int(LEN_F_date_list[i][:4])
    year = int(LEN_F_date_list[i][4:8])
    month = calendar.month_abbr[int(LEN_F_date_list[i][8])]
    day = int(LEN_F_date_list[i][-2:])
    field = gdata.pull_out_day_lentis(pr_lentis_future, real, year, month, day)*86400
    pr_field = gdata.extract_region(field, R_lim)
    pr_max_LF.append(np.max(pr_field.data))

pr_event = gdata.extract_region(event, R_lim)
pr_max_event = np.max(pr_event.data)



fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,4))
v1 = ax[0].violinplot([pr_max_Epast, pr_max_Epres], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['g','b']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

#ax[0].plot(1, Pevent, marker='o', color='r')
#ax[0].plot(1.6, Pevent, marker='o', color='r')
#ax[0].axhline(Pevent, color='r')
ax[0].axhline(np.mean(pr_max_Epast), color='g')
ax[0].axhline(np.mean(pr_max_Epres), color='b')
ax[0].set_xticks([1, 1.6])
ax[0].set_xticklabels(['Past', 'Present'])
ax[0].set_ylabel('Maximum Rainfall, mm')
#ax[0].set_title('Persistence')
ax[0].set_title('(a)', loc='left')

v2 = ax[1].violinplot([pr_max_LP, pr_max_LF], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v2['bodies'], colors):
    pc.set_facecolor(color)

#ax[1].plot(1, Pevent, marker='o', color='r')
#ax[1].plot(1.6, Pevent, marker='o', color='r')
#ax[1].axhline(Pevent, color='r')
ax[1].axhline(np.mean(pr_max_LP), color='b')
ax[1].axhline(np.mean(pr_max_LF), color='m')
ax[1].set_xticks([1, 1.6])
ax[1].set_xticklabels(['Present', 'Future'])
ax[1].set_ylabel('Maximum Rainfall, mm')
ax[1].set_title('(b)', loc='left')



## Presentation Plot
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,4))
v1 = ax[0].violinplot([pr_max_Epast, pr_max_Epres], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['g','b']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

ax[0].axhline(Pevent, color='r')
#ax[0].axhline(np.mean(P_past), color='g')
#ax[0].axhline(np.mean(P_present), color='b')
ax[0].set_xticks([1, 1.6])
ax[0].set_ylim([0,100])
ax[0].set_yticks([0, 20, 40, 60, 80])
ax[0].set_xticklabels(['Past', 'Present'])
ax[0].set_ylabel('ERA5')

v2 = ax[1].violinplot([pr_max_LP, pr_max_LF], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v2['bodies'], colors):
    pc.set_facecolor(color)

#ax[1].axhline(np.mean(PL_present), color='b')
#ax[1].axhline(np.mean(PL_future), color='m')
ax[1].set_xticks([1, 1.6])
ax[1].set_ylim([0,100])
ax[1].set_yticks([0, 20, 40, 60, 80])
ax[1].set_xticklabels(['Present', 'Future'])
ax[1].set_ylabel('KNMI-LENTIS')



