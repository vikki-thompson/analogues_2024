# %%
# Figures for analog paper, assessing persistence further for reviewers
#
# Original: vikki.thompson 26/07/2023
# Last Editted 16 May 2024

# %%
## Load neccessary libraries
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
import calendar
from scipy.stats.stats import pearsonr
plt.ion(); plt.show()

# %%
## Variables
region = [70, 30, 30, -30] # analog region
date = [2021, 'Jul', 14] # event date
coeff = 0.1

## Analogues
past_Y1 = 1950
past_Y2 = 1980
present_Y1 = 1993
present_Y2 = 2023

# %%
# Get ERA5 data: PSI, SLP, PR
psi_event = gdata.pull_out_day_era(gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(2021,2023))), region), date[0], date[1], date[2])

# Get data, find analogues based on streamfunction
psi_past = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(past_Y1,past_Y2))), region)
psi_present= gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(present_Y1,present_Y2))), region)
dates_past = gdata.analogs_datelist(psi_past, psi_event)[:30] # find analogs
dates_present = gdata.analogs_datelist(psi_present, psi_event)[:30]

# %%
def persistence(date, region, coeff):
    print(date)
    year = np.int(date[:4])
    mon = calendar.month_abbr[int(date[4])]
    day = np.int(date[5:7])
    psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(year,year+1))), region)
    psi_event = gdata.pull_out_day_era(psi, year, mon, day)
    C = gdata.correlation_coeffs(psi, psi_event)
    W = C.index(C[np.abs(C).argmin()])
    n = 0
    for i in np.arange(W-10, W+10):
        if i > len(C)-1:
            print('Too near end')
            return n
        elif C[i] < coeff: n+= 1
    return n

# %%
## ERA-5 persistence
Pevent = persistence('2021714', region, coeff)

P_past = []
for each_date in dates_past:
    P_past.append(persistence(each_date, region, coeff))

P_present = []
for each_date in dates_present:
    P_present.append(persistence(each_date, region, coeff))  

# %%
# how many days remaining below 0.1 CC? for event it is 4
# Violin plots of persistence
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
v2 = ax.violinplot([P_past, P_present], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
ax.plot(1, Pevent, marker='o', color='r')
ax.plot(1.6, Pevent, marker='o', color='r')
ax.axhline(Pevent, color='r')
ax.axhline(np.mean(P_past), color='g')
ax.axhline(np.mean(P_present), color='b')
ax.set_xticks([1, 1.6])
ax.set_xticklabels(['period 1', 'period 2'])
ax.set_ylabel('Days <'+str(coeff)+'CC')
ax.set_title('Persistence')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
nbins = np.arange(0,np.max([P_past, P_present]))
ax.hist(P_past, bins=nbins, color='g',alpha=0.5,label='Past')
ax.axvline(np.mean(P_past), color='g')
ax.hist(P_present, bins=nbins, color='b', alpha=0.5,label='Present')
ax.axvline(np.mean(P_present), color='b')
ax.axvline(Pevent, color='r',label='Event')
ax.set_xlabel('Days <'+str(coeff)+'CC')
ax.set_title('Persistence')
ax.legend()

# %%
## 14 day persistence - maps
### REFERENCE EVENT
date = '2021714'
reg = R1
year = np.int(date[:4])
mon = calendar.month_abbr[int(date[4])]
day = np.int(date[5:7])
psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(year,year+1))), reg)
psi_event = gdata.pull_out_day_era(psi, year, mon, day)
C = gdata.correlation_coeffs(psi, psi_event)
fig, axs = plt.subplots(nrows=1, ncols=9, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,3))
lats=psi_event.coord('latitude').points
lons=psi_event.coord('longitude').points
event_min = np.min(psi_event.data)
event_max = np.max(psi_event.data)
all_ax = axs.ravel()
for i in range(10,19):
    field = gdata.pull_out_day_era(psi, year, mon, i)
    c = all_ax[i-10].contourf(lons, lats, field.data, levels=np.linspace(event_min, event_max, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i-10].add_feature(cf.COASTLINE, linewidth=0.5)
    CORR = '{0:.2f}'.format(pearsonr(psi_event.data.reshape(40*61), field.data.reshape(40*61))[0])
    all_ax[i-10].set_title(str(i)+mon+str(year))
    all_ax[i-10].text(-27,25,'C: '+CORR)

date = dates_present[21] #20 is the 16day persistence event
reg = R1
year = np.int(date[:4])
mon = calendar.month_abbr[int(date[4])]
day = np.int(date[5:7])
psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(year,year+1))), reg)
psi_event = gdata.pull_out_day_era(psi, year, mon, day)
C = gdata.correlation_coeffs(psi, psi_event)
fig, axs = plt.subplots(nrows=1, ncols=9, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,3))
lats=psi_event.coord('latitude').points
lons=psi_event.coord('longitude').points
event_min = np.min(psi_event.data)
event_max = np.max(psi_event.data)
all_ax = axs.ravel()
for i in range(day-4,day+5):
    field = gdata.pull_out_day_era(psi, year, mon, i)
    c = all_ax[i-(day-4)].contourf(lons, lats, field.data, levels=np.linspace(event_min, event_max, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i-(day-4)].add_feature(cf.COASTLINE, linewidth=0.5)
    CORR = '{0:.2f}'.format(pearsonr(psi_event.data.reshape(40*61), field.data.reshape(40*61))[0])
    all_ax[i-(day-4)].set_title(str(i)+mon+str(year))
    all_ax[i-(day-4)].text(-27,25,'C: '+CORR)
   
date = dates_present[1] #20 is the 16day persistence event
reg = R1
year = np.int(date[:4])
mon = calendar.month_abbr[int(date[4])]
day = np.int(date[5:7])
psi = gdata.extract_region(gdata.extract_JJA(gdata.era5_mydata('psi250', range(year,year+1))), reg)
psi_event = gdata.pull_out_day_era(psi, year, mon, day)
C = gdata.correlation_coeffs(psi, psi_event)
fig, axs = plt.subplots(nrows=1, ncols=9, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,3))
lats=psi_event.coord('latitude').points
lons=psi_event.coord('longitude').points
event_min = np.min(psi_event.data)
event_max = np.max(psi_event.data)
all_ax = axs.ravel()
for i in range(day-4,day+5):
    field = gdata.pull_out_day_era(psi, year, mon, i)
    c = all_ax[i-(day-4)].contourf(lons, lats, field.data, levels=np.linspace(event_min, event_max, 20), cmap = plt.cm.get_cmap('RdBu_r'), transform=ccrs.PlateCarree(), extend='both')
    #pdata.plot_box(all_ax[i], R1) # plot box of analog region
    all_ax[i-(day-4)].add_feature(cf.COASTLINE, linewidth=0.5)
    CORR = '{0:.2f}'.format(pearsonr(psi_event.data.reshape(40*61), field.data.reshape(40*61))[0])
    all_ax[i-(day-4)].set_title(str(i)+mon+str(year))
    all_ax[i-(day-4)].text(-27,25,'C: '+CORR)

# %%
## LENTIS
psi_lentis = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), region)
psi_lentis_future = gdata.extract_region((gdata.extract_JJA(gdata.lentis_mydata(run='F'))), region)
psi_lentis = gdata.regrid(psi_lentis, psi_event)
psi_lentis_future = gdata.regrid(psi_lentis_future, psi_event)


# %%
Z = np.load('date_list_present_2400.npy')
'''
# swap for new data (which is presorted)
date_list = ['0072003809', '0112009729', '0142004809', '0162002728', '0192002726', '0212005726', '0262008802', '0292003718', '0292006816', '0292008811', '0302002804', '0312002725', '0322000719', '0342007814', '0362003807', '0382000812', '0412001717', '0492003716', '0502001818', '0502003730', '0522001730', '0522002713', '0542002812', '0592008828', '0612001810', '0612005725', '0622003825', '0642004723', '0652002821', '0652007803', '0692009807', '0702001813', '0712003801', '0752006721', '0762000728', '0772003729', '0772005801', '0782001830', '0792008720', '0802008819', '0812005813', '0822003808', '0832007831', '0852004726', '0862003808', '0862008713', '0892004730', '0912005728', '0932004824', '0942002806', '0942007721', '0982000729', '1002000806', '1002009720', '1012002721', '1012004815', '1022003725', '1032006820', '1042000726', '1052003724', '1092006803', '1102004818', '1122002827', '1132002721', '1142008725', '1152002802', '1162000715', '1172001814', '1202001802', '1212009817', '1242002818', '1252002729', '1282009725', '1302000725', '1322009727', '1342003728', '1372004823', '1392004720', '1412004727', '1442000831', '1442008712', '1442009817', '1452005819', '1452007813', '1452008811', '1452009817', '1482001718', '1492003827', '1522002716', '1532008829', '1542003822', '1542004717', '1542009815', '1572002811', '1572007803']
euc_dist = [383769470.0, 398308450.0, 384022370.0, 363571840.0, 398164540.0, 391334000.0, 341080220.0, 360858750.0, 334284450.0, 370511230.0, 384829100.0, 302420860.0, 313260830.0, 385307000.0, 387248860.0, 372859500.0, 391929180.0, 388205500.0, 396093500.0, 374797600.0, 384613630.0, 362002560.0, 378028100.0, 368107970.0, 330554620.0, 299685980.0, 384467400.0, 392783400.0, 369541250.0, 324847780.0, 373609920.0, 348826180.0, 347052830.0, 392992300.0, 375121300.0, 398415140.0, 386062750.0, 322342370.0, 334164600.0, 394218140.0, 392157200.0, 352917440.0, 385773470.0, 380210240.0, 390425200.0, 397124400.0, 367475780.0, 352604320.0, 393712350.0, 312966700.0, 369285730.0, 388517900.0, 386747400.0, 364849020.0, 386885280.0, 363963600.0, 397069660.0, 345820800.0, 378979800.0, 374814340.0, 353899170.0, 382420060.0, 379821570.0, 378096200.0, 397395420.0, 346168130.0, 307852960.0, 376201570.0, 392388960.0, 334605630.0, 352790530.0, 365396300.0, 382002720.0, 321811680.0, 375482340.0, 380925820.0, 398847170.0, 366802620.0, 360340450.0, 332393730.0, 397358750.0, 329283740.0, 313233440.0, 390194100.0, 396939260.0, 381506430.0, 367474660.0, 370281630.0, 352665540.0, 398021200.0, 390068260.0, 396419420.0, 342712300.0, 388975940.0, 397890560.0]
Z = [x for _,x in sorted(zip(euc_dist, date_list))]
date_list_future = ['0002075814', '0002080807', '0012077814', '0022076722', '0022084803', '0032076801', '0032077731', '0042076812', '0072078805', '0072079715', '0092080805', '0122076726', '0122077803', '0122082815', '0132076803', '0132077818', '0142081809', '0172076717', '0182083808', '0182084824', '0192077719', '0202076817', '0222079714', '0232077804', '0252081811', '0262075720', '0302080809', '0312078716', '0312084821', '0322083721', '0332083808', '0342079731', '0382076710', '0392078804', '0412076815', '0412080801', '0422080802', '0432076710', '0442077822', '0442081817', '0462081802', '0482076817', '0492078807', '0492081801', '0502079814', '0542080803', '0542084721', '0552084806', '0562084730', '0572078813', '0572082806', '0582082820', '0592083810', '0602077721', '0612076813', '0612079810', '0622083823', '0632079802', '0632080803', '0642079817', '0652075728', '0672084805', '0682075724', '0682083821', '0682084728', '0692081809', '0692082730', '0702077720', '0702079814', '0702082822', '0752076813', '0752077731', '0772077803', '0792076807', '0792080802', '0802079820', '0822076823', '0852078726', '0852083827', '0862075726', '0862084711', '0872081816', '0882076823', '0892076817', '0892083812', '0892084815', '0902081727', '0922079817', '0962080819', '0962084801', '0982084825', '0992083805', '1002078823', '1012077831', '1032078822', '1032079806', '1032082821', '1042076720', '1052077724', '1052078721']
euc_dist_future = [387434560.0, 372560640.0, 339371200.0, 383333220.0, 379542800.0, 385648420.0, 370060640.0, 384461660.0, 345797600.0, 372931140.0, 380576000.0, 376262050.0, 371955200.0, 298734820.0, 385781440.0, 353624740.0, 381472930.0, 375835900.0, 385588130.0, 275282200.0, 385672670.0, 353011600.0, 347564130.0, 323102720.0, 335224480.0, 385430000.0, 373379700.0, 378114020.0, 374401660.0, 380882240.0, 311983330.0, 377431330.0, 366168540.0, 327430000.0, 378197400.0, 387038900.0, 380519140.0, 360085470.0, 383034850.0, 371432640.0, 311750140.0, 369696450.0, 356651940.0, 357363420.0, 365603870.0, 315797340.0, 371705300.0, 369408260.0, 374296320.0, 376013150.0, 364923000.0, 384782100.0, 375115420.0, 385457380.0, 374187900.0, 371577760.0, 326246460.0, 377671700.0, 348616260.0, 354374080.0, 344209020.0, 369619170.0, 360292350.0, 354216030.0, 376531700.0, 362650460.0, 374187040.0, 365925000.0, 371743100.0, 351688830.0, 368676030.0, 340654560.0, 384840640.0, 381814240.0, 377411940.0, 386327170.0, 289593400.0, 372517440.0, 290284830.0, 362183460.0, 374905000.0, 334149570.0, 379100000.0, 366104350.0, 362510340.0, 346076960.0, 361016670.0, 374965200.0, 354852450.0, 368357730.0, 342274600.0, 333252060.0, 360338750.0, 314502880.0, 334556500.0, 379501760.0, 372532900.0, 370394600.0, 311106620.0, 386914700.0]
Z_future = [x for _,x in sorted(zip(euc_dist_future, date_list_future))]
'''

# %%
def persistence_lentis(psi, date, region, coeff):
    print(date)
    ens = np.int(date[:4])
    yr = np.int(date[4:8])
    mon = calendar.month_abbr[int(date[8])]
    day = np.int(date[9:11])
    #psi = gdata.regrid(psi, psi_event) # DOES IT MATTER IF GRID DIFFERENT?
    psi_event = gdata.pull_out_day_lentis(psi, ens, yr, mon, day)
    ens = np.int(date[1:4])
    psi = psi[ens,...]
    iris.coord_categorisation.add_year(psi, 'time')
    psi.extract(iris.Constraint(year=yr))
    C = gdata.correlation_coeffs(psi, psi_event)
    W = C.index(C[np.abs(C).argmin()])
    n = 0
    for i in np.arange(W-10, W+10):
        if i > len(C)-1:
            print('Too near end')
            return n
        elif C[i] < coeff: n+= 1
    return n

# %%
PL_present = []
psi = gdata.extract_region(gdata.extract_JJA(gdata.lentis_mydata()), region)  
count = 0 

'''
Z = np.load('/usr/people/thompson/WP1/revisions_figures/date_list_present_2400.npy')
import random
date_list2 = random.sample(Z.tolist(), 30)

for each_date in date_list2:
    print(count)
    PL_present.append(persistence_lentis(psi, each_date, region, coeff))
    print(PL_present)
    count+=1

np.save('/usr/people/thompson/WP1/revisions_figures/PL_present_2400.npy', PL_present)
'''

P1 = np.load('/usr/people/thompson/WP1/revisions_figures/PL_present_2400.npy')
P1_v2 = np.load('/usr/people/thompson/WP1/revisions_figures/PL_present_2400_v2.npy')
P1 = np.append(P1,P1_v2)[:30]


# %%

Z_future = np.load('/usr/people/thompson/WP1/revisions_figures/date_list_future_2400.npy')

import random
date_list2 = random.sample(Z_future.tolist(), 30)


PL_future = []
count = 0
psi = gdata.extract_region((gdata.extract_JJA(gdata.lentis_mydata(run='F'))), region)
for each_date in date_list2:
    print(count)
    PL_future.append(persistence_lentis(psi, each_date, region, coeff))  
    count+=1
    print(PL_future)

np.save('/usr/people/thompson/WP1/revisions_figures/PL_future_2400.npy', PL_future)
'''

P2 = np.load('/usr/people/thompson/WP1/revisions_figures/PL_future_2400.npy')
P2_v2 = np.load('/usr/people/thompson/WP1/revisions_figures/PL_future_2400_v2.npy')
#P2 = np.append(P2,P2_v2)[:30]
'''

# %%
P2 = [3, 3, 9, 11, 11, 9, 6, 6, 9, 1, 9, 6, 3, 6, 9, 8, 3, 8, 11, 6, 2, 14, 2, 2, 2, 3, 20, 6, 3, 12, 1]   

# %%

# how many days remaining below 0.1 CC? for event it is 4
# Violin plots of persistence
P_event = 4

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.5,2.5))
v1 = ax.violinplot([P_past, P_present], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

ax.plot(1, P_event, marker='o', color='r')
ax.plot(1.6, P_event, marker='o', color='r')
ax.axhline(np.mean(P1), color='b')
ax.axhline(np.mean(P2), color='m')
ax.set_xticks([1, 1.6])
ax.set_xticklabels(['Present', 'Future'])
ax.set_ylabel('Days')
ax.set_yticks([0, 5, 10, 15, 20])
ax.set_title('Persistence')
ax.set_ylim([0,17])
#ax.set_title('(a)', loc='left')


# %%
# how many days remaining below 0.1 CC? for event it is 4
# Violin plots of persistence
P_event = 4

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))

v1 = ax[0].violinplot([P_past, P_present], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['g','b']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

ax[0].plot(1, P_event, marker='o', color='r')
ax[0].plot(1.6, P_event, marker='o', color='r')
ax[0].axhline(np.mean(P1), color='g')
ax[0].axhline(np.mean(P2), color='b')
ax[0].set_xticks([1, 1.6])
ax[0].set_xticklabels(['Past', 'Present'])
ax[0].set_ylabel('Persistence (days)')
ax[0].set_yticks([0, 5, 10, 15, 20])
ax[0].set_title('ERA-5')
ax[0].set_ylim([0,17])

PL_present = [4, 8, 6, 7, 10, 5, 4, 4, 5, 5, 3, 2, 6, 7, 11, 4, 3, 5, 3, 5, 9, 3, 4, 2, 4, 7, 5, 2, 7, 6] 
PL_future = [5, 2, 6, 11, 9, 5, 6, 1, 5, 4, 3, 4, 6, 6, 11, 6, 5, 7, 5, 10, 8, 7, 3, 2, 5, 8, 8, 3, 1, 7]

v1 = ax[1].violinplot([PL_present, PL_future], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

ax[1].plot(1, P_event, marker='o', color='r')
ax[1].plot(1.6, P_event, marker='o', color='r')
ax[1].axhline(np.mean(P1), color='b')
ax[1].axhline(np.mean(P2), color='m')
ax[1].set_xticks([1, 1.6])
ax[1].set_xticklabels(['Present', 'Future'])
ax[1].set_yticks([0, 5, 10, 15, 20])
ax[1].set_title('KNMI-LENTIS 1')
ax[1].set_ylim([0,17])


v1 = ax[2].violinplot([P1, P2], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v1['bodies'], colors):
    pc.set_facecolor(color)

ax[2].plot(1, P_event, marker='o', color='r')
ax[2].plot(1.6, P_event, marker='o', color='r')
ax[2].axhline(np.mean(P1), color='b')
ax[2].axhline(np.mean(P2), color='m')
ax[2].set_xticks([1, 1.6])
ax[2].set_xticklabels(['Present', 'Future'])
ax[2].set_yticks([0, 5, 10, 15, 20])
ax[2].set_title('KNMI-LENTIS 2')
ax[2].set_ylim([0,17])

ax[0].set_title('(a)', loc='left')
ax[1].set_title('(b)', loc='left')
ax[2].set_title('(c)', loc='left')


# %%

v2 = ax[1].violinplot([PL_present, PL_future], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v2['bodies'], colors):
    pc.set_facecolor(color)

ax[1].plot(1, Pevent, marker='o', color='r')
ax[1].plot(1.6, Pevent, marker='o', color='r')
ax[1].axhline(np.mean(PL_present), color='b')
ax[1].axhline(np.mean(PL_future), color='m')
ax[1].set_xticks([1, 1.6])
ax[1].set_xticklabels(['Present', 'Future'])
ax[1].set_ylabel('Days')
ax[1].set_title('(b)', loc='left')


v2 = ax[2].violinplot([PL_present2, PL_future2], [1, 1.6], showmeans=False, showextrema=False, showmedians=False)
colors = ['b','m']
for pc, color in zip(v2['bodies'], colors):
    pc.set_facecolor(color)

ax[2].plot(1, Pevent, marker='o', color='r')
ax[2].plot(1.6, Pevent, marker='o', color='r')
ax[2].axhline(np.mean(PL_present2), color='b')
ax[2].axhline(np.mean(PL_future2), color='m')
ax[2].set_xticks([1, 1.6])
ax[2].set_xticklabels(['Present', 'Future'])
ax[2].set_ylabel('Days')
ax[2].set_title('(c)', loc='left')

plt.savefig(persistence.png)


