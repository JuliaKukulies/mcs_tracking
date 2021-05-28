## This python script performs recombines detected and segmented features based on cloud top temperatures and combines those with precipitation to filter convective systems. The tracking is performed per year, which means that MCS tracks at the boundary between two years are recognized as separate MCS tracks. 

import os,sys
import datetime
import glob
import numpy as n
import iris
import numpy as np
import pandas as pd
import tobac
import xarray as xr 
from scipy import ndimage
from scipy.ndimage import generate_binary_structure
from parameters import parameters_linking, dt, dxy, savedir, data_dir, precip_area, precip_threshold, cold_feature 

##########################################################################

## Import elevation file for 3000 m boundary 
dem = '/media/juli/Data/projects/data/elevation/elevation_600x350.nc'
elevations = xr.open_dataarray(dem)
# TP mask as coordinates 
dem_mask = elevations.where((elevations >= 3000) & (elevations.lat> 27) & (elevations.lat< 40) &(elevations.lon> 70)& (elevations.lon< 105))
dem_mask.coords['mask'] = (('lon', 'lat'), dem_mask)

#########################################################################

## Recombination of feature dataframes in each months (update framenumbers so that these are consecutive)
years = np.arange(2000,2020)
years = years.astype(str)


for year in years: 
    # read in HDF5 files with saved features for the respective year 
    file_list= glob.glob(savedir  + '/Features_cells_gpm_imerg_'+year+'??.h5')  
    file_list.sort()
    print('nr. of monthly feature files:', len(file_list), 'for year', year)

    i = 0 
    frames = 0 

    for file in file_list: 
        if i == 0:
            Features = pd.read_hdf(file, 'table')
            # read in data mask with segments for tracked cells 
            date= file[len(file)-12: len(file)-6]
            ds = Dataset(savedir+ '/Mask_Segmentation_gpm_imerg_'+date+'.nc')
            mask = np.array(ds['segmentation_mask'])  
            # update total nr of frames 
            frames += np.shape(mask)[0] -1
            i = 1 
            print('file for: ',date, 'rows: ',Features.shape[0], 'frames: ', frames)

        features = pd.read_hdf(file, 'table')
        # update frame number and make sure they are sequential
        features['idx'] = features['frame']
        features['frame'] = features['frame']  + frames

        # append dataframes 
        Features = Features.append(features, ignore_index=True)      
        # read in data mask with segments for tracked cells 
        date= file[len(file)-12: len(file)-6]
        ds = Dataset(savedir+ '/Mask_Segmentation_gpm_imerg_'+date+'.nc')
        mask = np.array(ds['segmentation_mask'])  
        #update total nr of frames
        frames += np.shape(mask)[0]
        print('file for: ',date, 'rows: ',features.shape[0], 'frames: ', frames)

    ## Perform trajectory linking with trackpy 
    Track=tobac.linking_trackpy(Features, np.zeros((10,10)),dt=dt,dxy=dxy,**parameters_linking)
    # remove nan values to only save the linked features and remove non-linked features from the list  
    tracks = Track[Track.cell >= 0]
    tracks.to_hdf(os.path.join(savedir,'Tracks_'+ str(year) +'_precip_tpflag.h5'),'table' ) 
    print('trajectory linking for year  '+ str(year) +'performed.')
    
    ########################################################################################################

    ###### Filter tracks with extra criteria (here: cold feature and hevay rain area during lifetime) ######

    ########################################################################################################
    tracks_cold_core = tracks.copy()        
    removed = 0
    # get extra information on brightness temperatures and precipitation 
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

            # check whether precip is in area of segmentation mask, where segmentation mask == feature number 
            maskfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tbbtracking_revised/Mask_Segmentation_'+str(year) + str(month)+str(month) + '.nc'
            precipfile = '/media/juli/Elements/gpm_v06/'+str(year)+'/gpm_imerg_'+ str(year)+str(month)+'_monthly.nc4'
            tbbfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/'+str(year)+'/merg_' + str(year)+str(month)+ '.nc4'
            tbbdata= xr.open_dataarray(tbbfile)
            tbbdata = tbbdata[:,1:,1:].T

            mask = xr.open_dataarray(maskfile)
            mask= mask[:,1:,1:].T        
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
                    rain_features = values[values >= precip_threshold].shape[0]
                    tracks['convective_precip'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)] = np.nansum(values[values >= precip_threshold])*0.5
                    tracks['rain_flag'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)]  = rain_features

                    # brightness temperatures cold core filter 
                    tbb_values = tbb.T.where(seg_mask.coords['mask'].values > 0)
                    arr= tbb_values.values.flatten()
                    values = arr[~np.isnan(arr)] # values contains the amount of grid cells with precip
                    tracks['mean_temp'][(tracks.feature == featureid) & (tracks.idx == idx) & (tracks.cell== cell)] = np.nanmean(values[values > 0])
                    tracks_cold_core['mean_temp'][(tracks_cold_core.feature == featureid) & (tracks_cold_core.idx == idx) & (tracks_cold_core.cell== cell)] = np.nanmean(values[values > 0])

                    if values[values > 0].min() <= cold_feature:
                        tbb_flag += 1 

                    # Elevation mask  
                    elevation_values = dem_mask.where(seg_mask.coords['mask'].values > 0)
                    arr= elevation_values.values.flatten()
                    values = arr[~np.isnan(arr)]

                    mountain_features = values[values >=3000].shape[0]
                    tracks['tp_flag'][(tracks.feature == featureid) & (tracks.idx == idx)& (tracks.cell== cell)] =  mountain_features
                    tracks_cold_core['tp_flag'][(tracks_cold_core.feature == featureid) & (tracks_cold_core.idx == idx)& (tracks_cold_core.cell== cell)] =  mountain_features

                    if rain_features >= precip_area: 
                        precipitation_flag += rain_features
                        
        if tbb_flag == 0:
            tracks_cold_core = tracks_cold_core.drop(tracks_cold_core[tracks_cold_core.cell == cell].index)
            tracks = tracks.drop(tracks[tracks.cell == cell].index)

        if precipitation_flag  ==  0:
            # remove corresponding cell from track dataframe 
            tracks = tracks.drop(tracks[tracks.cell == cell].index)
            removed += 1
            
    # save track files 
    tracks_cold_core.to_hdf(os.path.join(savedir,'Tracks_'+ str(year) +'_cold_core.h5'),'table' )            
    tracks.to_hdf(os.path.join(savedir,'Tracks_'+ str(year) +'_heavy_rain_core.h5'),'table' )
    print('tracks filtered for year  '+ str(year) +'performed.') 
