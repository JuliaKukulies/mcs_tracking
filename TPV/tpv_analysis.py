""" 
Functions for the analysis of the TPV TRACK database and its comparison to MCS observations. 

Created by Julia Kukulies,Feb 2021. 

"""
import os 
import re 
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr 
import cftime 

def get_composites(times):
    '''
    Get a list of datetime.datetime objects from array with numpy.datetime times   
    '''
    composites = []
    for t in times:
        t= datetime.utcfromtimestamp(t.astype(int) * 1e-9)
        year = t.year
        month = t.month
        day = t.day
        hour = t.hour
        dt= datetime(year,month,day,hour) 
        composites.append(dt)
    return composites 




def get_tracks(filename, year):
    """
    Read in text files with TPV tracks and create an organized pandas dataframe from it.

    Args: 
    filename(str): filename of the TPV file to read in 
    
    Returns: 
    
    tpv: pandas.DataFrame with tracs ordered by their ID
    
    """
    cols = ['time', 'lon', 'lat', 'vorticity', 'precip',  'lon_geo', 'lat_geo',  'geopm', 'id']

    tpv =  pd.DataFrame(columns=cols)

    with open(filename, 'r') as td:
        for line in td:
            if 'TRACK_ID' in line:
                for i in line.split():
                    if len(i) <= 4:
                        trackid= int(i) 

            if len(line) > 50:
                columns = []
                for i in line.split():
                    if i != '&':
                        if str(year) in i:
                            columns.append(i)
                        else:
                            columns.append(float(i))
                columns.append(trackid)
                tpv.loc[len(tpv)]= columns

    # get right data types for dataframe 
    tpv['time']=pd.to_datetime(tpv['time'],format='%Y%m%d%H')
    
    return tpv






def check_overlap_tpv(tpv,mcs, rad):
    """
    Check the overlap between TPV (before they are moving off) and MCS tracks. 
    
    Args: 
    tpv: pandas.DataFrame with TPV tracks 
    mcs: pandas.DataFrame with MCS tracks 
    
    Returns:
    tpv overlap: TPV with MCS in vicinity 
    tpv_no_mcs: 
    
    """
    from bisect import bisect_left
    tpv_overlap = 0         
    
    for tpv_id in np.unique(tpv.id.values):
        tpv_case= tpv[tpv.id == tpv_id]
        start = tpv_case.time.values[0]
        end = tpv_case.time.values[-1] + np.timedelta64(1,'D')
        tpv_case['lon'] =pd.to_numeric(tpv_case['lon'])
        tpv_case['lat'] =pd.to_numeric(tpv_case['lat'])

        # loop through mcs dates
        for cell in np.unique(mcs.cell.values):
            # do the whole thing per year to really get individual cell IDs
            subset = mcs[mcs.cell == cell]
            timeoverlap = subset[ (subset.timestr >= start) & (subset.timestr <= end)]
            if timeoverlap.shape[0]> 0:
                for t in timeoverlap.timestr.values:
                    mcslat = timeoverlap[timeoverlap.timestr== t].latitude.values[0]
                    mcslon = timeoverlap[timeoverlap.timestr== t].longitude.values[0]
                    # get the closes timestep in TPV set!
                    i = bisect_left(tpv_case.time.values, t)
                    tpvt= min(tpv_case.time.values[max(0, i-1): i+2], key=lambda ts: abs(t - ts))
                    tpvlon= tpv_case[tpv_case.time == tpvt].lon.values[0]
                    tpvlat= tpv_case[tpv_case.time == tpvt].lat.values[0]
                    # check for overlap in space
                    if (mcslat > tpvlat - rad) & (mcslat < tpvlat + rad) & (mcslon < tpvlon + rad) & (mcslon < tpvlon + rad):
                        tpv_overlap += 1
                        break 
            
    return tpv_overlap, tpv_no_mcs


def check_overlap(tpv,mcs):
    """
    Check the overlap between TPV (before they are moving off) and MCS tracks. 
    
    Args: 
    tpv: pandas.DataFrame with TPV tracks 
    mcs: pandas.DataFrame with MCS tracks 
    
    Returns: 

      
    """
    
    tpv_overlap = np.zeros(np.unique(tpv.id.values.shape[0]))            
    tpv_no_mcs = 0
    i  = 0
    unique_cells_tpv= np.array(())
    
    for tpv_id in np.unique(tpv.id.values):
        tpv_case= tpv[tpv.id == tpv_id]
        start = tpv_case.time.values[0]
        tpv_case['lon'] =pd.to_numeric(tpv_case['lon'])
        end = tpv_case.time.values[-1] + np.timedelta64(1,'D')

            
        # loop through mcs dates
        for cell in np.unique(mcs.cell.values):
            # do the whole thing per year to really get individual cell IDs
            subset = mcs[mcs.cell == cell]
            for t in np.arange(subset.shape[0]):
                time  = subset.timestr.values[t]
                if (start <= time <= end) == True:
                    tpv_overlap[i] += 1
                    unique_cells_tpv = np.append(unique_cells_tpv, cell)  
                    break
    
        if tpv_overlap[i] == 0:
            tpv_no_mcs += 1
        i += 1

    all_mcs = np.unique(mcs.cell.values).shape[0]
    all_tpv = np.unique(tpv.id.values).shape[0] 
    mcs_overlap= np.unique(np.array(unique_cells_tpv)).shape[0]
    mcs_no_tpv= np.setxor1d(mcs.cell.values, unique_cells_tpv).shape[0]
    
    return tpv_overlap, tpv_no_mcs, mcs_no_tpv, all_mcs, all_tpv, mcs_overlap



def check_overlap_tcs(tpv,mcs):
    """
    Check the overlap between TPV (before they are moving off) and MCS tracks. 
    
    Args: 
    tpv: pandas.DataFrame with TPV tracks 
    mcs: pandas.DataFrame with MCS tracks 
    
    Returns: 
      
    """

        
    mcs_count = np.zeros(np.unique(tpv.id.values.shape[0]))
    mcs_count_off = np.zeros(np.unique(tpv.id.values.shape[0]))
    tpv_no_mcs = 0 
    i  = 0
    unique_cells_tpv = np.array(())
    unique_cells_tpv_off = np.array(())

    
    for tpv_id in np.unique(tpv.id.values):
        tpv_case= tpv[tpv.id == tpv_id]
        start = tpv_case.time.values[0]
        tpv_case['lon'] =pd.to_numeric(tpv_case['lon'])
        if tpv_case.lon.values[-1] >= 105:
            tplon =  tpv_case[tpv_case.lon < 105].lon.values[-1]
            end = tpv_case[tpv_case.lon == tplon].time.values[0]
            
        else:
            end = tpv_case.time.values[-1]

            
        # loop through mcs dates
        for cell in np.unique(mcs.cell.values):
            # do the whole thing per year to really get individual cell IDs
            subset = mcs[mcs.cell == cell]
            for t in np.arange(subset.shape[0]):
                time  = subset.timestr.values[t]
                if (start <= time <= end) == True:
                    if tpv_case.lon.values[-1] > 105:
                        mcs_count_off[i] += 1
                        unique_cells_tpv = np.append(unique_cells_tpv, cell)  
                        break
                    else:
                        mcs_count[i] += 1
                        unique_cells_tpv_off = np.append(unique_cells_tpv_off, cell)
                        break
  
        all_mcs = np.unique(mcs.cell.values).shape[0]

        if mcs_count[i] == 0 and mcs_count_off[i] == 0 :
            tpv_no_mcs += 1
        i += 1
        
    all_tpv = np.unique(tpv.id.values).shape[0]
    overlap_mcs= np.unique(np.array(unique_cells_tpv)).shape[0]
    overlap_mcs_off = np.unique(np.array(unique_cells_tpv_off)).shape[0]
    mcs_no_tpv= np.setxor1d(mcs.cell.values, np.append(unique_cells_tpv, unique_cells_tpv_off)).shape[0] 

    return mcs_count, mcs_count_off, tpv_no_mcs, mcs_no_tpv, all_mcs, all_tpv, overlap_mcs, overlap_mcs_off 


def zero_padding(time):
    if time >= 10:
        time = str(time)
    else:
        time = '0' + str(time)
    return time




def fixed_location(time,loclon,loclat, path = None):
    '''
    Create netcdf files with fixed locations given a weather systems
    center location and time. These netcdf files can then be used for composite analysis.
    
    Args:
    time(datetime.datetime): timepoint for which the composite should be created
    loclon(float): longitude value of center location 
    loclat(float): latitude value sof center location 
    path(str): path where to store the created file 

    '''
    # get string corresponding to datetime value 
    month = zero_padding(time.month)
    day = zero_padding(time.day)
    hour = zero_padding(time.hour)

    # deal with leap years 
    if time.day == 29:
        day = zero_padding(28)
    
    #################### get corresponding files ##########################################
    # GPM (at 0.1 degree)
    gpm_file= '/media/juli/Elements/gpm_v06/'+ str(time.year)+'/gpm_imerg_' +str(time.year) + month+'_monthly.nc4'

    # Tb (at 4 km)
    datestr= str(time.year) + month+ day+hour
    ncep_file= '/media/juli/Data/projects/data/satellite_data/ncep/ctt/merg_'+ datestr +'_4km-pixel.nc4'

    # check if files exist
    if os.path.isfile(gpm_file) is True and os.path.isfile(ncep_file) is True:
        try:
            precip = xr.open_dataset(gpm_file)
            tb = xr.open_dataset(ncep_file)

            ################## Extract region given TPV center location and timestep ###############
            rad = 5

            # check if location is still in downloaded observation files
            if loclon+rad > precip.lon.values.max() or loclat+rad > precip.lat.values.max() or loclon-rad < precip.lon.values.min() or loclat-rad < precip.lat.values.min():
                print('TPV center location out of dimension of the observation files.')
            else:

                # get closest timesteps
                timestep1 = precip.sel(time = cftime.DatetimeNoLeap(time.year, time.month, int(day), time.hour,0,0,0,0,0), method ='nearest')
                timestep2 = precip.sel(time = cftime.DatetimeNoLeap(time.year, time.month, int(day), time.hour,30,0,0,0,0), method ='nearest')
                # average over timesteps in same hour 
                timestep = (timestep1 + timestep2) / 2
                # select fixed location with 3 deg radius
                prec_fixed = timestep.precipitationCal.isel(lat = (precip.lat <= loclat+ rad) & (precip.lat >=loclat- rad), lon = (precip.lon <= loclon+ rad) & (precip.lon >=loclon- 3))

                tb_fixed = tb.Tb.isel(lat = (tb.lat <= loclat+ rad) & (tb.lat >=loclat- rad), lon = (tb.lon <= loclon+ rad) & (tb.lon >=loclon- rad) ) 
                # average over timesteps in same hour 
                tb_fixed = tb_fixed.mean(dim = 'time')

                # save as netcdf
                if path is None:
                    path = ''
                tb_fixed.to_netcdf(path + datestr + '_tb.nc4')
                prec_fixed.to_netcdf(path + datestr + '_precip.nc4')

                precip.close()
                tb.close()
                tb_fixed.close()
                prec_fixed.close()
                        

                
        except IOError:
            print(datestr, ' --> file corrupted')
            
            
