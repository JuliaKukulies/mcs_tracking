""" 
Functions for the analysis of the TPV TRACK database and its comparison to MCS observations. 

Created by Julia Kukulies,Feb 2021. 

"""

import re 
import numpy as np
import pandas as pd



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


