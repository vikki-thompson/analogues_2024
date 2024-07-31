'''
Created 27 Feb 2023
Editted 27 Feb 2023


Functions for plotting data
>conda activate butterfly 
'''

# Load neccessary libraries
import subprocess
import iris
import iris.coord_categorisation as icc
from iris.coord_categorisation import add_season_membership
import numpy as np
import matplotlib.pyplot as plt
#import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy as cart
import cartopy.feature as cf
import glob
import matplotlib.cm as mpl_cm
import scipy.io
import xarray as xr
import netCDF4 as nc
import iris.coords
from iris.util import equalise_attributes
from iris.util import unify_time_units
import cartopy.feature as cf
import matplotlib.ticker as tick
import calendar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import sys
sys.path.append('/usr/people/thompson/WP1')
import functions_get_data as gdata



def mymap(axis, cube):
    lats=cube.coord('latitude').points
    lons=cube.coord('longitude').points
    c = axis.contourf(lons, lats, cube.data, labels=np.linspace(), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axis.add_feature(cf.BORDERS)
    axis.add_feature(cf.COASTLINE)
    #cbar_ax = fig.add_axes([0.85, 0.3, 0.02, 0.4])
    #fig.colorbar(c, cax=cbar_ax, ticks=[crange[0], 0, crange[-1]])
    return

def plot_box(axs, bdry):
    axs.plot([bdry[3], bdry[2]], [bdry[1], bdry[1]],'k')
    axs.plot([bdry[3], bdry[2]], [bdry[0], bdry[0]],'k')
    axs.plot([bdry[3], bdry[3]], [bdry[1], bdry[0]],'k')
    axs.plot([bdry[2], bdry[2]], [bdry[1], bdry[0]],'k')
    return

def all_var_map(R1, date):
    # var data
    event_psi250 = gdata.var_event_data('psi250', R1, date)
    event_slp = gdata.var_event_data('msl', R1, date)/1000
    event_prec = gdata.var_event_data('tp', R1, date)
    # plot data
    plt.ion(); plt.show()
    fig, axs = plt.subplots(nrows=1, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
    lats=event_psi250.coord('latitude').points
    lons=event_psi250.coord('longitude').points
    c = axs[0].contourf(lons, lats, event_psi250.data, levels=np.linspace(-8.2e7, 1e7, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[0].set_title('psi250')
    cb = fig.colorbar(c, ax=axs[0], orientation="horizontal", pad=0.2)
    tick_locator1 = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator1
    cb.update_ticks()
    lats=event_slp.coord('latitude').points
    lons=event_slp.coord('longitude').points
    c = axs[1].contourf(lons, lats, event_slp.data, levels=np.linspace(100, 103.3, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[1].set_title('slp')
    cd = fig.colorbar(c, ax=axs[1], orientation="horizontal", pad=0.2)
    tick_locator2 = ticker.MaxNLocator(nbins=5)
    cd.locator = tick_locator2
    cd.update_ticks()
    c = axs[2].contourf(lons, lats, event_prec.data, levels=np.linspace(0, 100, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[2].set_title('prec')
    ce = fig.colorbar(c, ax=axs[2], orientation="horizontal", pad=0.2)
    tick_locator3 = ticker.MaxNLocator(nbins=5)
    ce.locator = tick_locator3
    ce.update_ticks()
    for ax in axs:
        ax.add_feature(cf.BORDERS)
        ax.add_feature(cf.COASTLINE)
    #axs[0].colorbar(c, fraction=0.046, pad=0.04)
    plt.suptitle(str(date[2])+date[1]+str(date[0]))
    return
    
def single_var_map(R1, date, var):
    # var : 'psi250', 'msl', 'tp'
    if var == 'msl':
        event_var = gdata.var_event_data(var, R1, date)/1000
    else:
        event_var = gdata.var_event_data(var, R1, date)
    # plot data
    plt.ion(); plt.show()
    fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
    lats=event_var.coord('latitude').points
    lons=event_var.coord('longitude').points
    lev_max = np.max([-np.min(event_var.data), np.max(event_var.data)])
    c = axs.contourf(lons, lats, event_var.data, levels=np.linspace(-lev_max, lev_max, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    cb = fig.colorbar(c, ax=axs, orientation="horizontal", pad=0.2)
    tick_locator1 = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator1
    cb.update_ticks()
    axs.add_feature(cf.BORDERS)
    axs.add_feature(cf.COASTLINE)
    #axs.colorbar(c, fraction=0.046, pad=0.04)
    plt.suptitle(var+', '+str(date[2])+date[1]+str(date[0]))
    return


def all_var_map_v2(R1, date):
    '''
    Plots two maps:
    SLP contours with TP fill, and shaded streamfunction.
    For given region and date
    '''
    # var data
    event_psi250 = gdata.var_event_data('psi250', R1, date)
    event_prec = gdata.var_event_data('tp', R1, date)
    # plot data
    plt.ion(); plt.show()
    fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
    lats=event_prec.coord('latitude').points
    lons=event_prec.coord('longitude').points
    c = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 100, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
    lats=event_psi250.coord('latitude').points
    lons=event_psi250.coord('longitude').points
    c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
    axs.add_feature(cf.BORDERS)
    axs.add_feature(cf.COASTLINE)
    return
    

def all_var_map_lentis_v2(R1, date):
    '''
    Plots two maps:
    SLP contours with TP fill, and shaded streamfunction.
    For given region and date
    '''
    # var data
    real = int(date[:4])
    year = int(date[4:8])
    month = calendar.month_abbr[int(date[8])]
    day = int(date[-2:])
    pr_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_data()), R1)
    event_prec = gdata.pull_out_day_lentis(pr_lentis, real, year, month, day)*86400
    psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), R1)
    event_psi250 = gdata.pull_out_day_lentis(psi_lentis, real, year, month, day)*86400 
    # plot data
    plt.ion(); plt.show()
    fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))
    lats=event_prec.coord('latitude').points
    lons=event_prec.coord('longitude').points
    c = axs.contourf(lons, lats, event_prec.data, levels=np.linspace(0, 100, 10), cmap = plt.cm.get_cmap('Blues'), transform=ccrs.PlateCarree(), extend='both')
    lats=event_psi250.coord('latitude').points
    lons=event_psi250.coord('longitude').points
    c = axs.contour(lons, lats, event_psi250.data, levels=np.linspace(-6e7, -2E7, 12), cmap = plt.cm.get_cmap('autumn'), transform=ccrs.PlateCarree(), extend='both')
    axs.add_feature(cf.BORDERS)
    axs.add_feature(cf.COASTLINE)
    return
    

## Composites
def composite_plots(event, comp1, comp2):
    # comp1 -> past (or present)
    # comp2 -> present (or future)
    fig, axs = plt.subplots(nrows=1, ncols=4, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,4))
    lats=event.coord('latitude').points
    lons=event.coord('longitude').points
    c = axs[0].contourf(lons, lats, event.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[0].add_feature(cf.BORDERS)
    axs[0].add_feature(cf.COASTLINE)
    axs[0].set_title('Event')
    fig.colorbar(c, ax=axs[0], format=tick.FormatStrFormatter('%.0f'), ticks = tick.FixedLocator([]), fraction=0.046, pad=0.04)
    c = axs[1].contourf(lons, lats, comp1.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[1].add_feature(cf.BORDERS)
    axs[1].add_feature(cf.COASTLINE)
    axs[1].set_title('Composite, Period 1')
    fig.colorbar(c, ax=axs[1], format=tick.FormatStrFormatter('%.0f'), ticks = tick.FixedLocator([]), fraction=0.046, pad=0.04)
    c = axs[2].contourf(lons, lats, comp2.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[2].add_feature(cf.BORDERS)
    axs[2].add_feature(cf.COASTLINE)
    axs[2].set_title('Composite, Period 2')
    fig.colorbar(c, ax=axs[2], format=tick.FormatStrFormatter('%.0f'), ticks = tick.FixedLocator([]), fraction=0.046, pad=0.04)
    anom = comp2.data - comp1.data
    c = axs[3].contourf(lons, lats, anom, levels=np.linspace(-abs(np.min(anom.data)), abs(np.min(anom.data)), 10), cmap = plt.cm.get_cmap('PiYG'), transform=ccrs.PlateCarree(), extend='both')
    axs[3].add_feature(cf.BORDERS)
    axs[3].add_feature(cf.COASTLINE)
    axs[3].set_title('Anomaly, Period 2 - Period 1')
    fig.colorbar(c, ax=axs[3], fraction=0.046, pad=0.04)
    plt.tight_layout()
    return


def composite(event, comp):
    fig, axs = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,4))
    lats=event.coord('latitude').points
    lons=event.coord('longitude').points
    c = axs[0].contourf(lons, lats, event.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[0].add_feature(cf.BORDERS)
    axs[0].add_feature(cf.COASTLINE)
    axs[0].set_title('Event')
    fig.colorbar(c, ax=axs[0], format=tick.FormatStrFormatter('%.0f'), ticks = tick.FixedLocator([]), fraction=0.046, pad=0.04)
    c = axs[1].contourf(lons, lats, comp.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs[1].add_feature(cf.BORDERS)
    axs[1].add_feature(cf.COASTLINE)
    axs[1].set_title('Composite of analogues')
    fig.colorbar(c, ax=axs[1], format=tick.FormatStrFormatter('%.0f'),  ticks = tick.FixedLocator([]), fraction=0.046, pad=0.04)
    plt.tight_layout()
    return



### Analogue Measures
def quality(Q1_all, Q2_all, Q1_ana, Q2_ana, Q1_event, Q2_event):
    ' Violin plot ' 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    v1 = ax.violinplot([Q1_all, Q2_all], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
    for b in v1['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further right than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
        b.set_color('b')
    v2 = ax.violinplot([Q1_ana, Q2_ana], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
    for b in v2['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further right than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
        b.set_color('orange')
    ax.plot(1, Q1_event, marker='o', color='r')
    ax.plot(1.6, Q2_event, marker='o', color='r')
    ax.set_xticks([1, 1.6])
    ax.set_xticklabels(['period 1', 'period 2'])
    ax.set_ylabel('250hPa Streamfunction, kg/ms')
    ax.set_title('Analog Quality')
    return

def persistence(P1_all, P2_all, P1_ana, P2_ana, P_event):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    v1 = ax.violinplot([P1_all, P2_all], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
    for b in v1['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further right than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
        b.set_color('b')
    v2 = ax.violinplot([P1_ana, P2_ana], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
    for b in v2['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further right than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
        b.set_color('orange')
    ax.axhline(P_event, color='r')
    ax.set_xticks([1, 1.6])
    ax.set_xticklabels(['period 1', 'period 2'])
    ax.set_ylabel('Days')
    ax.set_title('Analog Persistence')
    ax.set_ylim([0, 20])
    ax.text(1.7, 18, 'All Days', color='b')
    ax.text(1.7, 17, 'Analogues', color='orange')
    ax.text(1.7, 19, 'Event', color='red')
    return

def map_box(R1, event):
    fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,4))
    lats=event.coord('latitude').points
    lons=event.coord('longitude').points
    c = axs.contourf(lons, lats, event.data, levels=np.linspace(np.min(event.data), np.max(event.data), 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    axs.add_feature(cf.BORDERS)
    axs.add_feature(cf.COASTLINE)
    axs.set_title('Event Region: '+str(R1[1])+' to '+str(R1[0])+' N, '+str(R1[3])+' to '+str(R1[2])+' E')
    plot_box(axs, R1)
    axs.set_xlim([-60, 60])
    axs.set_ylim([10, 80])
    return


## All analogs

def plot_all_analogs(var, date_list, reg, date, title, Y1, Y2):
    '''
    Plots field 'var' for all dates in 'date_list' over region 'reg'
    'date' = date of key event (used to define plot limits)
    '''
    my_var = ['psi250', 'u250', 'v250', 'z500', 'msl']
    if var in my_var:
        var_cube = gdata.era5_mydata(var, range(2021, 2023))
    else:
        var_cube = gdata.era5_data(var, range(2021, 2023))
    var_cube = gdata.extract_region(var_cube, reg)
    fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
    event = gdata.pull_out_day_era(var_cube, date[0], date[1], date[2])
    lats=event.coord('latitude').points
    lons=event.coord('longitude').points
    event_min = np.min(event.data)
    event_max = np.max(event.data)
    if var in my_var:
        var_cube = gdata.era5_mydata(var, range(Y1, Y2))
    else:
        var_cube = gdata.era5_data(var, range(Y1, Y2))
    var_cube = gdata.extract_region(var_cube, reg)
    # subplots
    all_ax = axs.ravel()
    for i, each in enumerate(date_list):
        YY = np.int(each[:4])
        MM = calendar.month_abbr[int(each[4])]
        DD = int(each[-2:])
        event = gdata.pull_out_day_era(var_cube, YY, MM, DD)
        c = all_ax[i].contourf(lons, lats, event.data, levels=np.linspace(event_min, event_max, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
        plot_box(all_ax[i], reg) # plot box of analog region
        all_ax[i].add_feature(cf.BORDERS)
        all_ax[i].add_feature(cf.COASTLINE)
        all_ax[i].set_title(each)
    fig.suptitle(title)
    plt.tight_layout()   
    return


def all_analogs(event, psi, date_list):
    '''
    Plots field 'var' for all dates in 'date_list' over region 'reg'
    'date' = date of key event (used to define plot limits)
    '''
    fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
    lats=event.coord('latitude').points
    lons=event.coord('longitude').points
    event_min = np.min(event.data)
    event_max = np.max(event.data)
    # subplots
    all_ax = axs.ravel()
    for i, each in enumerate(date_list):
        YY = np.int(each[:4])
        MM = calendar.month_abbr[int(each[4])]
        DD = int(each[-2:])
        analog = gdata.pull_out_day_era(psi, YY, MM, DD)
        c = all_ax[i].contourf(lons, lats, analog.data, levels=np.linspace(event_min, event_max, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
        all_ax[i].add_feature(cf.BORDERS)
        all_ax[i].add_feature(cf.COASTLINE)
        all_ax[i].set_title(each)
    plt.tight_layout()   
    return 


def plot_all_analogs_lentis(cube_lentis, date_list, reg, title):
    '''
    Plots field 'var' for all dates in 'date_list' over region 'reg'
    'date' = date of key event (used to define plot limits)
    '''
    fig, axs = plt.subplots(nrows=5, ncols=6, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,10))
    lats=cube_lentis.coord('latitude').points
    lons=cube_lentis.coord('longitude').points
    RRR = np.int(date_list[0][:3])
    YY = np.int(date_list[0][3:7])
    MM = calendar.month_abbr[int(date_list[0][7])]
    DD = int(date_list[0][-2:])
    event = gdata.pull_out_day_lentis(cube_lentis, RRR, YY, MM, DD)
    event_min = np.min(event.data)
    event_max = np.max(event.data)
    # subplots
    all_ax = axs.ravel()
    for i, each in enumerate(date_list):
        RRR = np.int(each[:3])
        YY = np.int(each[3:7])
        MM = calendar.month_abbr[int(each[7])]
        DD = int(each[-2:])
        event = gdata.pull_out_day_lentis(cube_lentis, RRR, YY, MM, DD)
        c = all_ax[i].contourf(lons, lats, event.data, levels=np.linspace(event_min, event_max, 10), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
        plot_box(all_ax[i], reg) # plot box of analog region
        all_ax[i].add_feature(cf.BORDERS)
        all_ax[i].add_feature(cf.COASTLINE)
        all_ax[i].set_title(each)
    fig.suptitle(title)
    plt.tight_layout()   
    return 
