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
dem_mask = elevations.where(elevations >= 3000)
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

years = np.arange(2000,2020)
years = years.astype(str)

# perform trajectory linking per year

for year in years: 
    # read in HDF5 files with saved features
    file_list= glob.glob(savedir  + '/Features_cells_'+year+'*.h5')  
    file_list.sort()
    print('nr. of monthly feature files:', len(file_list), 'for year', year)

    i = 0 
    frames = 0 

    for file in file_list: 
        if i == 0:
            Features = pd.read_hdf(file, 'table')
            # read in data mask with segments for tracked cells 
            date= file[len(file)-9: len(file)-3]
            ds = Dataset(savedir+ '/Mask_Segmentation_'+date+'.nc')
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
        date= file[len(file)-9: len(file)-3]
        ds = Dataset(savedir+ '/Mask_Segmentation_'+date+'.nc')
        mask = np.array(ds['segmentation_mask'])  
        #update total nr of frames
        frames += np.shape(mask)[0]
        print('file for: ',date, 'rows: ',features.shape[0], 'frames: ', frames)


    ## Perform trajectory linking with trackpy 
    Track=tobac.linking_trackpy(Features,tbb,dt=dt,dxy=dxy,**parameters_linking)
    # remove nan values to only save the linked features 
    Track = Track[Track.cell >= 0]
    Track.to_hdf(os.path.join(savedir,'Tracks__'+ year +'.h5'),'table')

    tracks= Track 
    print(np.unique(tracks.cell.values).shape[0], '  unique cloud cells. Features:', tracks.shape[0])


    ########################## Filter ###################################
    ## Cold core 
    # once during the lifetime a cold core of <= 200 K needs to be present 
    ######
    for i in np.unique(tracks.cell.values):
        subset = tracks[tracks.cell == i]
        if subset[subset.threshold_value <= 200].shape[0] == 0 :
            tracks.drop(tracks.loc[tracks['cell']== i].index, inplace=True)

    tracks.to_hdf(os.path.join(savedir,'Tracks_'+year+'_cold_core.h5'),'table')  
    print('cold core filtered.', tracks.shape)

    ########################################## Heavy rain core#######################################
    removed = 0
    tracks['rain_flag'] = 0 
    tracks['tp_flag'] = 0
    tracks['total_precip']= 0
    tracks['convective_precip'] = 0

    pd.options.mode.chained_assignment = None

    tracks['timestr'] = pd.to_datetime(tracks.time)

    # loop through cells in detected feature frame 
    for cell in np.unique(tracks.cell.values):
        subset = tracks[tracks.cell == cell]
        #print('checking heavy rain cores for cell:', cell, subset.shape)
        precipitation_flag = 0

        # loop through timesteps of features for specific cell 
        for idx in subset.idx.values: 
            # idx is the timestep index for respective timestep or mask file
            # open corresponding precip and mask file
            year = subset.timestr[subset.idx == idx].dt.year.values[0] 
            month = subset.timestr[subset.idx == idx].dt.month.values[0]
            if len(str(month))== 1: 
                month= '0' + str(month)

            # check whether precip is in area of segmentation mask, where segmentation mask == feature number 
            maskfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs/Mask_Segmentation_'+str(year) + str(month) + '.nc'
            precipfile = '/media/juli/Elements/gpm_v06/'+str(year)+'/gpm_imerg_'+ str(year)+str(month)+'_monthly.nc4'

            mask = xr.open_dataarray(maskfile)
            mask= mask[:,:,:].T        
            precip = xr.open_dataarray(precipfile)
            precip = precip[:,1:,1:].T

            # check whether corresponding file exist and whether it has the same shape
            if mask.shape[2] == precip.shape[2]:

                # get right timestep frames 
                seg= mask[:,:, idx]
                prec = precip[:,:, idx]

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
                    precip_values = prec.T.where(seg_mask.coords['mask'].values > 1)
                    arr= precip_values.values.flatten()
                    values = arr[~np.isnan(arr)] # values contains the amount of grid cells with precip
                    total_precip = np.nansum(values[values > 0]) * 0.5

                    tracks['total_precip'][tracks.feature == featureid] = total_precip 

                    rain_features = values[values >= 3].shape[0]
                    tracks['convective_precip'][tracks.feature == featureid] = np.nansum(values[values >= 5])*0.5
                    tracks['rain_flag'][tracks.feature == featureid]  = rain_features

                    # Elevation mask  
                    elevation_values = dem_mask.where(seg_mask.coords['mask'].values > 1)
                    arr= elevation_values.values.flatten()
                    values = arr[~np.isnan(arr)]

                    mountain_features = values[values >=3000].shape[0]
                    tracks['tp_flag'][tracks.feature == featureid] =  mountain_features

                    if rain_features >= 1: 
                        precipitation_flag += rain_features
            else:
                np.savetxt(savedir+ 'shape_'+ str(year) +str(month)+'txt', [precip.shape, mask.shape])

        if precipitation_flag  ==  0:
            # remove corresponding cell from track dataframe 
            tracks = tracks.drop(tracks[tracks.cell == cell].index)
            print(cell, ' removed.')
            removed += 1 

        #else:
            #print('heavy rain core present in:  ', cell, rain_features)

    #print(removed, ' cells removed in total.')
    tracks.to_hdf(os.path.join(savedir,'Tracks_'+ str(year) +'_heavyraincore3mmocc.h5'),'table' ) 

    print('trajectory linking for year  '+ str(year) +'performed.') 

