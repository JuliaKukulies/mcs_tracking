## This python script performs a tracking cloud-features based tobac (Heikenfeld et al., 2019) on NCEP brightness temperatures. The specific tracking here aims for mesoscale convective systems in the Tibetan Plateau region, whereby a meso-scale convective system is defined as a cloud complex with
# *
# *
# * 

# Note that the input files are monthly files of merged 30-min data aggregated in one file, so that the feature detection and segmentation step of the tracking can be performed per month, whereby the linking is is conducted in the way that it can even link systems at the boundaries of two months.  


# created by Julia Kukulies, julia.kukulies@gu.se 
#####################################################################################################################
import warnings

warnings.simplefilter(action = "ignore", category = RuntimeWarning)

import dask.array as da 
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

# temporal and spatial resolution (in seconds and meter)
dt= 1800
dxy = 14126.0


## Feature detection
# Dictionary containing keyword options (could also be directly given to the function)
parameters_features={}
parameters_features['position_threshold']='weighted_diff' # diff between specific value and threshold for weighting when finding the center location (instead of just mean lon/lat)
parameters_features['min_distance']=0 # minimum distance between features 
parameters_features['sigma_threshold']=0.5 # for slightly smoothing (gaussian filter)
parameters_features['n_erosion_threshold']=0 # pixel erosion (for more robust results)
parameters_features['threshold']=[220, 218, 216, 214, 212, 210, 205, 200,195, 190] #mm/h, step-wise threshold for feature detection 
parameters_features['n_min_threshold']=250 # minimum nr of contiguous pixels for thresholds, 10 pixels = ca. 2000 km2, 50 pixel ca. 10 000 km2
parameters_features['target']= 'minimum'


## Segmentation
# Dictionary containing keyword arguments for segmentation step:
parameters_segmentation={}
parameters_segmentation['target'] = 'minimum'
parameters_segmentation['method']='watershed'
parameters_segmentation['threshold']=221  # mm/h mixing ratio (until which threshold the area is taken into account)


## Tracking 
# Dictionary containing keyword arguments for the linking step:
parameters_linking={}
parameters_linking['adaptive_stop']=0.2
parameters_linking['adaptive_step']=0.95
parameters_linking['extrapolate']=0
parameters_linking['order']=1
parameters_linking['subnetwork_size']= 1000 # maximum size of subnetwork used for linking 
parameters_linking['memory']=0
parameters_linking['time_cell_min']= 6*dt 
parameters_linking['method_linking']='predict'
#parameters_linking['method_detection']='threshold'
parameters_linking['v_max']= 10
#parameters_linking['d_min']=2000
parameters_linking['d_min']=4*dxy # four times the grid spacing ?

############################################################## Tracking : Feature detection and Segmentation ###################################################################################################
import glob 
# list with all files by month
file_list= glob.glob(data_dir + '????/merg_??????.nc4')  
print('files in dataset:  ', len(file_list))
file_list.sort()


for f in file_list[99::]:
    i = f[len(data_dir)+10:-4]
    month = f[len(data_dir)+14:-4]
    
    # if threshold should change with season 
    #if month in ['01','02','12']:
        #parameters_segmentation['threshold'] = 235
        #print('winter month: segmentation threshold switched to 235k.')

    print('start process for file.....', i, month )
    ## DATA PREPARATION
    Precip=iris.load_cube(f, 'brightness_temperature')
    # set values to NaN
    #Precip.data[Precip.data > 300] = np.nan
    #Precip.data[Precip.data < 0 ] = np.nan

    for i in range(Precip.shape[0]):
        part = Precip.data[i,:,:]
        part[part < 0] = np.nan
        Precip.data[i] = part 

    #data = Precip.data
    #data = da.where(data < 0, data, np.nan)
    #Precip = Precip.copy(data=data)

    # FEATURE DETECTION
    print('starting feature detection based on multiple thresholds')
    Features=tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
    print('feature detection done')
    Features.to_hdf(os.path.join(savedir,'tbbtracking_revised/Features_' + str(i) + '.h5'),'table')
    print('features saved', Features.shape)
    
    # SEGMENTATION 
    print('Starting segmentation based on surface precipitation')
    Mask,Features_Precip=tobac.segmentation_2D(Features,Precip,dxy,**parameters_segmentation)
    print('segmentation based on surface precipitation performed, start saving results to fs')
    iris.save([Mask],os.path.join(savedir,'tbbtracking_revised/Mask_Segmentation_' + str(i) + '.nc'),zlib=True,complevel=4)                
    Features_Precip.to_hdf(os.path.join(savedir,'tbbtracking_revised/Features_cells_' + str(i) + '.h5'),'table')
    print('segmentation surface precipitation performed and saved')

    del Precip
    del Features
    del Mask
    del Features_Precip
    gc.collect()
