## This python script performs recombines detected and segmented features based on cloud top temperatures and combines those with precipitation to filter convective systems

import xarray as xr 
import warnings
warnings.simplefilter(action = "ignore", category = RuntimeWarning)

import glob 
import iris
import numpy as np
import pandas as pd
import os,sys
import datetime
from netCDF4 import Dataset
import tobac

from scipy import ndimage
from scipy.ndimage import generate_binary_structure

##########################################################################

## Import elevation file for 3000 m boundary 

dem = '/media/juli/Data/projects/data/elevation/elevation_600x350.nc'
elevations = xr.open_dataarray(dem)

# mask as coordinates 
dem_mask = elevations.where((elevations >= 3000) & (elevations.lat> 27) & (elevations.lat< 40) &(elevations.lon> 70)& (elevations.lon< 105))
dem_mask.coords['mask'] = (('lon', 'lat'), dem_mask)


#################################### Tracking parameters ########################


# Dictionary containing keyword arguments for the linking step:                                         
parameters_linking={}          
dt = 1800
dxy = 14126.0

s = generate_binary_structure(2,2)

parameters_linking['adaptive_stop']=0.2                                                                 
parameters_linking['adaptive_step']=0.95                                                                
parameters_linking['extrapolate']=0                                                                     
parameters_linking['order']=1

parameters_linking['subnetwork_size']= 1000000 # maximum size of subnetwork used for linking               
parameters_linking['memory']= 1                                                                          
parameters_linking['time_cell_min']= 6*dt                                                              
parameters_linking['method_linking']='predict'                                                          
#parameters_linking['method_detection']='threshold'                                                    
parameters_linking['v_max']= 100                                                                        
#parameters_linking['d_min']=2000                                                                      
#parameters_linking['d_min']=4*dxy # four times the grid spacing ?                                       
            
## Recombination of feature dataframes (update framenumbers)
savedir= '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs'


# get brightness temp file 
from netCDF4 import Dataset
f = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/2001/merg_200106.nc4'
ds= Dataset(f)
tbb = np.array(ds['Tb']) 

years = np.arange(2019,2020)
years = years.astype(str)

for year in years:
    # remove nan values to only save the linked features
    Track = pd.read_hdf(os.path.join(savedir,'Tracks__'+ year +'.h5'),'table')
    Track = Track[Track.cell >= 0]

    tracks= Track
    tracks_cold_core = Track.copy()        
    removed = 0
    tracks['rain_flag'] = 0 
    tracks['tp_flag'] = 0
    tracks['total_precip']= 0
    tracks['convective_precip'] = 0

    tracks['mean_temp'] = 0
    tracks_cold_core['mean_temp'] = 0
    tracks_cold_core['tp_flag'] = 0 

    pd.options.mode.chained_assignment = None

    tracks['timestr'] = pd.to_datetime(tracks.time)

    # loop through cells in detected feature frame 
    for cell in np.unique(tracks.cell.values):
        subset = tracks[tracks.cell == cell]
        #print('checking heavy rain cores for cell:', cell, subset.shape)
        precipitation_flag = 0
        tbb_flag = 0 

        # loop through timesteps of features for specific cell 
        for idx in subset.idx.values: 
            # idx is the timestep index for respective timestep or mask file
            # open corresponding precip and mask file
            year = subset.timestr[subset.idx == idx].dt.year.values[0] 
            month = subset.timestr[subset.idx == idx].dt.month.values[0]
            if len(str(month))== 1: 
                month= '0' + str(month)
            if int(month)<= 6:

                # check whether precip is in area of segmentation mask, where segmentation mask == feature number 
                maskfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs/Mask_Segmentation_'+str(year) + str(month) + '.nc'
                precipfile = '/media/juli/Elements/gpm_v06/'+str(year)+'/gpm_imerg_'+ str(year)+str(month)+'_monthly.nc4'
                tbbfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/'+str(year)+'/merg_' + str(year)+str(month)+ '.nc4'
                tbbdata= xr.open_dataarray(tbbfile)
                tbbdata = tbbdata[:,1:,1:].T

                mask = xr.open_dataarray(maskfile)
                mask= mask[:,:,:].T        
                precip = xr.open_dataarray(precipfile)
                precip = precip[:,1:,1:].T

                # check whether corresponding file exist and whether it has the same shape
                if mask.shape[2] == precip.shape[2]:

                    # get right timestep frames 
                    seg= mask[:,:, idx]
                    prec = precip[:,:, idx]
                    tbb = tbbdata[:,:,idx].T

                    # get feature ID for frame 
                    featureid= subset.feature[subset.idx== idx].values[0]

                    # get locations from all contiguous fields in segmentation mask which belong to cell which contain the featureid from tracked cells
                    # note: this is necessary, because the segmentation mask still contains different feature ids within one cloud cell 
                    labels, nr = ndimage.label(seg, structure = s)

                    if featureid not in seg:
                        np.savetxt(savedir+ 'features_'+ str(year) +str(month)+ str(cell) + '.txt', [idx, featureid])
                        continue
                    else:

                        label = np.unique(labels[ seg == featureid])[0]
                        seg_mask = seg.where(labels == label)

                        # create mask as coordinates 
                        seg_mask.coords['mask'] = (('lon', 'lat'), seg_mask)
                        # apply mask on precip data to extract precip values for feature in cell 
                        precip_values = prec.T.where(seg_mask.coords['mask'].values > 0)
                        arr= precip_values.values.flatten()
                        values = arr[~np.isnan(arr)] # values contains the amount of grid cells with precip
                        total_precip = np.nansum(values[values > 0]) * 0.5
                        tracks['total_precip'][(tracks.feature == featureid) & (tracks.idx == idx) & (tracks.cell== cell)] = total_precip 
                        rain_features = values[values >= 5].shape[0]
                        tracks['convective_precip'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)] = np.nansum(values[values >= 5])*0.5
                        tracks['rain_flag'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)]  = rain_features

                        # brightness temperatures cold core filter 
                        tbb_values = tbb.T.where(seg_mask.coords['mask'].values > 0)
                        arr= tbb_values.values.flatten()
                        values = arr[~np.isnan(arr)] # values contains the amount of grid cells with precip
                        tracks['mean_temp'][(tracks.feature == featureid) & (tracks.idx == idx) & (tracks.cell== cell)] = np.nanmean(values[values > 0])
                        tracks_cold_core['mean_temp'][(tracks_cold_core.feature == featureid) & (tracks_cold_core.idx == idx) & (tracks_cold_core.cell== cell)] = np.nanmean(values[values > 0])

                        if values[values > 0].shape[0]> 0:
                             if values[values > 0].min() <=200:
                                  tbb_flag += 1 

                        # Elevation mask  
                        elevation_values = dem_mask.where(seg_mask.coords['mask'].values > 0)
                        arr= elevation_values.values.flatten()
                        values = arr[~np.isnan(arr)]

                        mountain_features = values[values >=3000].shape[0]
                        tracks['tp_flag'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)] =  mountain_features
                        tracks_cold_core['tp_flag'][(tracks_cold_core.feature == featureid) & (tracks_cold_core.idx == idx)& (tracks_cold_core.cell== cell)] =  mountain_features

                        if rain_features >= 3: 
                            precipitation_flag += rain_features

                else:
                    np.savetxt(savedir+ 'shape_'+ str(year) +str(month)+'txt', [precip.shape, mask.shape])


            if tbb_flag == 0:
                tracks_cold_core = tracks_cold_core.drop(tracks_cold_core[tracks_cold_core.cell == cell].index)
                tracks = tracks.drop(tracks[tracks.cell == cell].index)


            if precipitation_flag  ==  0:
                # remove corresponding cell from track dataframe 
                tracks = tracks.drop(tracks[tracks.cell == cell].index)
                #print(cell, ' removed.')
                removed += 1 
            #else:
                #print('heavy rain core present in:  ', cell, rain_features)

        #print(removed, ' cells removed in total.')
        tracks_cold_core.to_hdf(os.path.join(savedir,'tracks/Tracks_'+ str(year) +'_cold_core.h5'),'table' )            
        tracks.to_hdf(os.path.join(savedir,'tracks/Tracks_'+ str(year) +'_heavy_rain_core.h5'),'table' ) 

        print('trajectory linking for year  '+ str(year) +'performed.') 

