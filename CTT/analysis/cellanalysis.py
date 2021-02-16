## This python scripts collects functions for the analysis of tracked MCS

import numpy as np
import pandas as pd


# function to get seasonal curve of features
def get_seasonal_features(tracks):
    tracks['month']= tracks.timestr.dt.month
    seasonal=[]
    for m in np.arange(1,13):
        monthly_count = tracks[tracks.month== m].shape[0]   # for frequency 
        #meanvalue= np.nanmean(preciptracks.total_precip.values)    # for mean values of specific variable 
        seasonal.append(monthly_count)
    return seasonal



def count_tracks(tracks):
    count = 0 
    for y in np.unique(tracks.timestr.dt.year):
        subset = tracks[tracks.timestr.dt.year == y]
        count += np.unique(subset.cell.values).shape[0]
    return count 



def get_seasonal_curve(tracks):
    tracks['month']= tracks.timestr.dt.month
    td= np.timedelta64(0, 'ns')
    init_features = tracks[tracks.time_cell == td]
    seasonal=[]    
    for m in np.arange(1,13):
        monthly_count = init_features[init_features.month== m].shape[0] 
        seasonal.append(monthly_count)
    return seasonal



def get_lifetime(tracks):
    lt= []
    for y in np.arange(2000,2020):
        ytracks = tracks[tracks.timestr.dt.year== y]
        for cell in np.unique(ytracks.cell.values):
            hours= ytracks[ytracks.cell== cell].shape[0] * 0.5
            lt.append(hours)
    lt = np.array(lt)
    lt= np.histogram(lt, bins= np.arange(3,51)[::2]) 
    print('lifetime histo calculated.')
    return lt



def get_init(tracks):
    init_lats= []
    init_lons= []
    diss_lats= []
    diss_lons= []
    for y in np.arange(2000,2020):
        ytracks = tracks[tracks.timestr.dt.year== y]
        for cell in np.unique(ytracks.cell.values):
            subset= tracks[tracks.cell == cell]
            init_lats.append(subset.latitude.values[0])
            init_lons.append(subset.longitude.values[0])
            diss_lats.append(subset.latitude.values[-1])
            diss_lons.append(subset.longitude.values[-1])

    return np.array(init_lats), np.array(init_lons), np.array(diss_lats), np.array(diss_lons)



def get_area(tracks):
    a= []
    for y in np.arange(2000,2020):
        ytracks = tracks[tracks.timestr.dt.year== y]
        for cell in np.unique(ytracks.cell.values):
            area= np.nanmean(np.array(ytracks[ytracks.cell== cell].ncells.values))
            a.append(area)
    a = np.array(a)
    a = np.histogram(a, bins=np.arange(0,4050,50))
    print('area histo calculated.')
    return a


# divide into north-moving, east-moving and other 

def propagation_dir(tracks):
    pd.options.mode.chained_assignment = None 
    
    tracks['dir'] = 0 
    for c in np.unique(tracks.cell.values):
        cell= tracks[tracks.cell == c]

        west_east= np.nanmean(cell.longitude.values[-3:-1])  - np.nanmean(cell.longitude.values[-0:2])
        north_south= np.nanmean(cell.latitude.values[-3:-1])  - np.nanmean(cell.latitude.values[-0:2])
                                                                         
        if north_south > west_east:
            if np.nanmean(cell.latitude.values[0:2]) < np.nanmean(cell.latitude.values[-3:-1]):
                tracks['dir'][tracks.cell == c] =  'N'
            elif np.nanmean(cell.latitude.values[0:2]) > np.nanmean(cell.latitude.values[-3:-1]):
                tracks['dir'][tracks.cell == c] =  'S'
                
        elif north_south < west_east:
            if np.nanmean(cell.longitude.values[0:2]) < np.nanmean(cell.longitude.values[-3:-1]):
                tracks['dir'][tracks.cell == c] =  'E'
            elif np.nanmean(cell.longitude.values[0:2]) > np.nanmean(cell.longitude.values[-3:-1]):
                tracks['dir'][tracks.cell == c] =  'W'
    return tracks 





def get_v(tracks):
    v= []
    for cell in np.unique(tracks.cell.values):
        ps = np.nanmean(tracks[tracks.cell== cell].v.values)
        v.append(ps)
    v = np.array(v)
    v = np.histogram(v, bins=np.arange(6,31)[::2]) 
    print('propagation speed histo calculated.')
    return v


def get_diurnal_cycle(preciptracks):
    preciptracks['hour']= preciptracks.timestr.dt.hour
    diurnal=[]
    for h in np.arange(0,23):
        count = preciptracks[preciptracks.hour == h].shape[0]
        #meanvalue= np.nanmean(preciptracks.total_precip.values)
        diurnal.append(count)
    return diurnal



def get_diurnal_init(preciptracks):
    preciptracks['hour']= preciptracks.timestr.dt.hour
    diurnal=[]
    for cell in np.unique(preciptracks.cell.values):
        init_hour = preciptracks[preciptracks.cell == cell].hour.values[0]
        diurnal.append(init_hour)
    return diurnal



def get_diurnal_diss(preciptracks):
    preciptracks['hour']= preciptracks.timestr.dt.hour
    diurnal=[]
    for cell in np.unique(preciptracks.cell.values):
        init_hour = preciptracks[preciptracks.cell == cell].hour.values[-1]
        diurnal.append(init_hour)
    return diurnal



def get_max_values(tracks):
    tracks['hour']= tracks.timestr.dt.hour
    peak_frame = pd.DataFrame(columns = tracks.columns)
    rain_peak = []
    for cell in np.unique(tracks.cell.values):
        subset= tracks[tracks.cell == cell]
        peak = np.nanmax(subset.total_precip.values)
        hour = subset[subset.total_precip == peak].hour.values[0]
        # add row to dataframe 
        #peak_frame = pd.concat([peak_frame, rain_peak ], ignore_index=True)
        rain_peak.append(hour)
    rain_histo = np.histogram(rain_peak, bins = np.arange(0,24))
    return rain_histo[0]


def divide_data(tracks):
    tp_tracks= pd.DataFrame(columns = tracks.columns)
    surrounding_tracks= pd.DataFrame(columns = tracks.columns)
    for y in np.unique(tracks.timestr.dt.year):
        print(y)
        subset = tracks[tracks.timestr.dt.year == y]
        for i in np.unique(subset.cell.values):
            cell = subset[subset.cell == i]
            if np.sum(cell.tp_flag.values) < 100:
                surrounding_tracks= surrounding_tracks.append(cell)
            else:
                tp_tracks= tp_tracks.append(cell)
                
    return tp_tracks, surrounding_tracks

