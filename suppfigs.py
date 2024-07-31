# %%
# Figures for analog paper, ERA5
#
# All analogues fields - streamfunction & precip
#
# Original: vikki.thompson 26/07/2023
# Last Editted 13 May 2024

# %%
import warnings
warnings.filterwarnings('ignore')

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
import glob
plt.ion(); plt.show()

# %%
## Variables
region = [70, 30, 30, -30] # analog region
date = [2021, 'Jul', 14] # event date

## Analogues
past_Y1 = 1950
past_Y2 = 1980
present_Y1 = 1993
present_Y2 = 2023

# %%
# ERA-5 DATA AND ANALOGUES DATES
psi_event = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2021,2023))), region)
event = gdata.pull_out_day_era(psi_event, date[0], date[1], date[2])
lats=event.coord('latitude').points
lons=event.coord('longitude').points

psi_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), region)
psi_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), region)
pr_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(past_Y1,past_Y2))), region)
pr_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_data('tp', range(present_Y1,present_Y2))), region)
dates_past = gdata.analogs_datelist(psi_past, event)[:30] # find analogs
dates_present = gdata.analogs_datelist(psi_present, event)[:30]

# %%
# LENTIS DATA 
R1 = region
psi_era = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2013,2023))), R1)
psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R1)
psi_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_mydata(run='F'))), R1)
psi_lentis = gdata.regrid(psi_lentis, event)
psi_lentis_future = gdata.regrid(psi_lentis_future, event)

pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), R1)
pr_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_data(run='2K'))), R1)
pr_lentis = gdata.regrid(pr_lentis, event)
pr_lentis_future = gdata.regrid(pr_lentis_future, event)

# LENTIS ANALOGUE DATES
LEN_P_dates = ['10472008803', '10572002823', '10472008804', '10692005726', '10692005725', '11132002725', '11132002726', '10572002822', '10542000715', '10962000807', '11122001719', '11332002806', '10342005819', '11242000719', '11102000725', '10832001830', '10732001808', '10572007803', '10982001721', '10692001811', '10752009817', '10692001810', '10842008817', '10752000831', '10862008720', '10192006816', '10872009817', '10742009816', '11122001720', '10472008802']
LEN_F_dates = ['50292084824', '51452083727', '50122076823', '50362083827', '51522082815', '50412075721', '50222077724', '51332081802', '51362083808', '50732077831', '50632080803', '50372079821', '50142077804', '50512083823', '50752078804', '50322083805', '51192081816', '51342078822', '50722081811', '51102077814', '51692077731', '50792077726', '50262084825', '51242075728', '51202078805', '51052084815', '51482078820', '51372079714', '50932080803', '51432082822']
#LEN_F_dates = ['11402004726', '10392007726', '11152003710', '10752009722', '10922007722', '10332003807', '10902006728', '11582008719', '10662007728', '11472009704', '10162002727', '10662002712', '11642002714', '11632002706', '11642005823', '10732004724']



# %%
for i in range(len(LEN_F_dates)):
    if LEN_F_dates[i] in LEN_P_dates:
        print(Y)

# %%

### LENTIS Analogue Plots
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    real = int(LEN_F_dates[i][:4])
    year = int(LEN_F_dates[i][4:8])
    month = calendar.month_abbr[int(LEN_F_dates[i][8])]
    day = int(LEN_F_dates[i][-2:])
    print(real, year, month, day)
    field = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
    c = all_ax[i].contourf(lons, lats, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='max')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])


# %%

fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=-.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[0, 30, 60], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
plt.suptitle('LENTIS Future Analogues: Rainfall')

##
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    real = int(LEN_P_dates[i][:4])
    year = int(LEN_P_dates[i][4:8])
    month = calendar.month_abbr[int(LEN_P_dates[i][8])]
    day = int(LEN_P_dates[i][-2:])
    field = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
    c = all_ax[i].contourf(lons, lats, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='max')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])

fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=-.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[0, 30, 60], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
plt.suptitle('LENTIS Present Analogues: Rainfall')

# %%
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    real = int(LEN_P_dates[i][:4])
    year = int(LEN_P_dates[i][4:8])
    month = calendar.month_abbr[int(LEN_P_dates[i][8])]
    day = int(LEN_P_dates[i][-2:])
    field = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)
    c = all_ax[i].contourf(lons, lats, field.data/1000000, levels=np.linspace(-80, -10, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))

fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
plt.suptitle('LENTIS Present Analogues: 250 hPa Streamfunction')

# %%
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(5):
    print(i)
    real = int(LEN_F_dates[i][:4])
    year = int(LEN_F_dates[i][4:8])
    month = calendar.month_abbr[int(LEN_F_dates[i][8])]
    day = int(LEN_F_dates[i][-2:])
    field = gdata.pull_out_day_lentis(psi_lentis_future, real, year, month, day)
    c = all_ax[i].contourf(lons, lats, field.data/1000000, levels=np.linspace(-80, -10, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))
    

fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
plt.suptitle('LENTIS Future Analogues: 250 hPa Streamfunction')


# %%
### ERA-5 Analogues Plots
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    year = int(dates_past[i][:4])
    month = calendar.month_abbr[int(dates_past[i][4])]
    day = int(dates_past[i][-2:])
    analog_date = [year, month, day]
    field = gdata.var_event_data('psi250', region, analog_date)
    field = gdata.regrid(field, event)
    c = all_ax[i].contourf(lons, lats, field.data/1000000, levels=np.linspace(-80, -10, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))


fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
plt.suptitle('ERA-5 Past Analogues: 250 hPa Streamfunction')

# %%

##
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    year = int(dates_present[i][:4])
    month = calendar.month_abbr[int(dates_present[i][4])]
    day = int(dates_present[i][-2:])
    analog_date = [year, month, day]
    field = gdata.var_event_data('psi250', region, analog_date)
    field = gdata.regrid(field, event)
    c = all_ax[i].contourf(lons, lats, field.data/1000000, levels=np.linspace(-80, -10, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))


fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[-80, -60, -40, -20], orientation='horizontal')
cbar_ax.set_xlabel('250 hPa Streamfunction (x$10^6$ m$^2$/s)')
plt.suptitle('ERA-5 Present Analogues: 250 hPa Streamfunction')


# %%
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    year = int(dates_past[i][:4])
    month = calendar.month_abbr[int(dates_past[i][4])]
    day = int(dates_past[i][-2:])
    analog_date = [year, month, day]
    field = gdata.var_event_data('tp', region, analog_date)
    field = gdata.regrid(field, event)
    c = all_ax[i].contourf(lons, lats, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='max')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])


fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[0, 30, 60], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
plt.suptitle('ERA-5 Past Analogues: Rainfall')


# %%
fig, axs = plt.subplots(nrows=6, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,10))
all_ax = axs.ravel()
for i in range(30):
    print(i)
    year = int(dates_present[i][:4])
    month = calendar.month_abbr[int(dates_present[i][4])]
    day = int(dates_present[i][-2:])
    analog_date = [year, month, day]
    field = gdata.var_event_data('tp', region, analog_date)
    field = gdata.regrid(field, event)
    c = all_ax[i].contourf(lons, lats, field.data, levels=np.linspace(0, 60, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='max')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i].add_feature(cf.COASTLINE, linewidth=.5)
    all_ax[i].set_title(str(day)+' '+month+' '+str(year))
    all_ax[i].set_xlim([-5, 22])
    all_ax[i].set_ylim([43, 60])


fig.subplots_adjust(bottom=0.12, wspace=0.2, hspace=.2)
cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.01])
fig.colorbar(c, cax=cbar_ax, ticks=[0, 30, 60], orientation='horizontal')
cbar_ax.set_xlabel('Daily Rainfall (mm)')
plt.suptitle('ERA-5 Present Analogues: Rainfall')


