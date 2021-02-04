""" 
Functions for the analysis of the TPV TRACK database and its comparison to MCS observations. 

Created by Julia Kukulies,Feb 2021. 

"""

import re 
import numpy as np
import pandas as pd


def get_tracks(filename):
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
                    if len(i) == 4:
                        trackid= int(i) 

            if len(line) > 50:
                columns = []
                for i in line.split():
                    if i != '&':
                        if '2008' in i:
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
    Check the overlap between TPV and MCS tracks. 
    
    Args: 
    tpv: pandas.DataFrame with TPV tracks 
    mcs: pandas.DataFrame with MCS tracks 
    
    Returns: 
    mcs_count: number of MCS that occur when there is a TPV
    no_mcs_count: number of TPV that occur without any MCS 
    no_tpv_count: number of MCS that occur without any TPV 
    """
    
    mcs_count = np.zeros(np.unique(tpv.id.values).size)
    no_mcs_count = 0 
    i= 0 

    for tpv_id in np.unique(tpv.id.values):
        tpv_case= tpv[tpv.id == tpv_id]
        start = tpv_case.time.values[0]
        end = tpv_case.time.values[-1]
        # loop through mcs dates 

        for cell in np.unique(mcs.cell.values):
            # do the whole thing per year to really get individual cell IDs
            subset = mcs[mcs.cell == cell]
            for t in np.arange(subset.shape[0]):
                time  = subset.timestr.values[t]
                if (start <= time <= end) == True:
                    mcs_count[i] += 1
                    break

        if mcs_count[i] == 0:
            no_mcs_count += 1

        i += 1

    no_tpv_count = np.unique(mcs.cell.values).shape[0] - mcs_count.sum()
    return mcs_count, no_mcs_count, no_tpv_count