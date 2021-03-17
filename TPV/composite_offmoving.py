## This script creates fixed location composites of precipitation and brightness temperatures
## for Tibetan Plateau vortices that move off the TP 

import pandas as pd
from tpv_analysis import get_tracks, fixed_location, get_composites

path = 'composites/off-moving/'

for y in np.arange(2000,2016):
    filename='tpv_files_JuliaK/tpv_'+str(y)+'.real_time'
    tpv = get_tracks(filename, y)
    tpv['geopm'] = pd.to_numeric(tpv.geopm)
    for cell in np.unique(tpv.id.values):
        print(y,cell)
        subset = tpv[tpv.id == cell]
        geopmin= subset.geopm.values.min()
        # get only off-moving
        if subset.lon.values[-1]> 105:
            ### maturity ###
            time = get_composites(subset[subset.geopm == geopmin].time.values)[0]
            loclon = subset[subset.geopm == geopmin].lon.values[0]
            loclat = subset[subset.geopm == geopmin].lat.values[0]
            # create fixed location files 
            fixed_location(time, loclon,loclat, path)
            
            ### initiation ###
            time = get_composites(subset.time.values)[0]
            loclon = subset.lon.values[0]
            loclat = subset.lat.values[0]
            # create fixed location files 
            fixed_location(time, loclon,loclat, path)

            ### dissipation ###
            time = get_composites(subset.time.values)[-1]
            loclon = subset.lon.values[-1]
            loclat = subset.lat.values[-1]
            # create fixed location files 
            fixed_location(time, loclon,loclat, path)

