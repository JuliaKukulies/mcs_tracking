## This python script performs a tracking cloud-features based tobac (Heikenfeld et al., 2019) on NCEP brightness temperatures. The specific tracking here aims for mesoscale convective systems in the Tibetan Plateau region, whereby a meso-scale convective system is defined as a cloud complex with
# *
# *
# * 

# Note that the input files are monthly files of merged 30-min data aggregated in one file, so that the feature detection and segmentation step of the tracking can be performed per month, whereby the linking is is conducted in the way that it can even link systems at the boundaries of two months.  


# created by Julia Kukulies, julia.kukulies@gu.se 
#####################################################################################################################
import warnings

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

import xarray as xr 
import iris
import numpy as np
import pandas as pd
import datetime
from netCDF4 import Dataset
import tobac
import glob
import os 
import gc
############################### Parameters ###############################################################################################
# specify output directory 

data_dir = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/'
savedir = data_dir + 'Save/'
os.makedirs(savedir,exist_ok=True)

# temporal and spatial resolution
dt= 1800
dxy = 14126.0



## Feature detection
# Dictionary containing keyword options (could also be directly given to the function)
parameters_features={}
parameters_features['position_threshold']='weighted_diff' # diff between specific value and threshold for weighting when finding the center location (instead of just mean lon/lat)
parameters_features['min_distance']=0 # minimum distance between features 
parameters_features['sigma_threshold']=0.5 # for slightly smoothing (gaussian filter)
parameters_features['n_erosion_threshold']=0 # pixel erosion (for more robust results)
parameters_features['threshold']=[230, 225, 220, 215, 210, 205, 200, 195] #mm/h, step-wise threshold for feature detection 
parameters_features['n_min_threshold']= 40 # minimum nr of contiguous pixels for thresholds, 10 pixels = ca. 2000 km2, 50 pixel ca. 10 000 km2
parameters_features['target']= 'minimum'



## Segmentation
# Dictionary containing keyword arguments for segmentation step:
parameters_segmentation={}
parameters_segmentation['target'] = 'minimum'
parameters_segmentation['method']='watershed'
parameters_segmentation['threshold']=245  # mm/h mixing ratio (until which threshold the area is taken into account)

##########################################################################                                                                                     

## Import elevation file for 3000 m boundary                                                                                                                   

dem = '/media/juli/Data/projects/data/elevation/elevation_600x350.nc'
elevations = xr.open_dataarray(dem)
elev = elevations.data.T 

# mask as coordinates                                                                                                                                          
dem_mask = elevations.where(elevations >= 3000)
dem_mask.coords['mask'] = (('lon', 'lat'), dem_mask)



############################################################## Tracking : Feature detection and Segmentation ###################################################################################################

import glob 
# list with all files by month
file_list= glob.glob(data_dir + '????/merg_??????.nc4')  
print('files in dataset:  ', len(file_list))
file_list.sort()


for f in file_list[0::]:
    i = f[len(data_dir)+10:-4]
    month = f[len(data_dir)+14:-4]

    if month in ['01','02','12']:
        parameters_segmentation['threshold'] = 235
        print('winter month: segmentation threshold switched to 235k.')

    print('start process for file.....', i, month )
    ## DATA PREPARATION
    precip = xr.open_dataarray(f)
    precip = precip[:,1:,1:]
    
    # elevation mask
    #for t in np.arange(precip.shape[2]):
        #prec =precip[:,:,t]
        #tbb_values = prec.where(dem_mask.coords['mask'].values > 1)
        #precip.data[:,:,t] = tbb_values
        
    # fancy indexing 
    precip.data[:, elev <= 3000]  = np.nan

    #to iris cube 
    Precip = precip.to_iris()
    print(Precip.shape)
    
    # FEATURE DETECTION
    print('starting feature detection based on multiple thresholds')
    Features=tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
    print('feature detection done')
    Features.to_hdf(os.path.join(savedir,'smallscale40/Features_' + str(i) + '.h5'),'table')
    print('features saved', Features.shape)
    
    # SEGMENTATION 
    print('Starting segmentation based on surface precipitation')
    Mask,Features_Precip=tobac.segmentation_2D(Features,Precip,dxy,**parameters_segmentation)
    print('segmentation based on surface precipitation performed, start saving results to fs')
    iris.save([Mask],os.path.join(savedir,'smallscale40/Mask_Segmentation_' + str(i) + '.nc'),zlib=True,complevel=4)                
    Features_Precip.to_hdf(os.path.join(savedir,'smallscale40/Features_cells_' + str(i) + '.h5'),'table')
    print('segmentation surface precipitation performed and saved')

    del Precip
    del precip
    del Features
    del Mask
    del Features_Precip
    gc.collect()
