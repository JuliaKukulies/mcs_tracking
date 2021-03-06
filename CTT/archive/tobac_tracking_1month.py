## This python script performs a tracking cloud-features based tobac (Heikenfeld et al., 2019) on NCEP brightness temperatures. The specific tracking here aims for mesoscale convective systems in the Tibetan Plateau region, whereby a meso-scale convective system is defined as a cloud complex with
# *
# *
# * 

# Note that the input files are monthly files of merged 30-min data aggregated in one file, so that the feature detection and segmentation step of the tracking can be performed per month, whereby the linking is is conducted in the way that it can even link systems at the boundaries of two months.  


# created by Julia Kukulies, julia.kukulies@gu.se 
#####################################################################################################################

import iris
import numpy as np
import pandas as pd
import datetime
from netCDF4 import Dataset
import tobac
import glob
import os 

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
parameters_features['threshold']=[230, 225, 220, 215, 210, 205, 200] #mm/h, step-wise threshold for feature detection 
parameters_features['n_min_threshold']=50 # minimum nr of contiguous pixels for thresholds, 10 pixels = ca. 2000 km2, 50 pixel ca. 10 000 km2
parameters_features['target']= 'minimum'



## Segmentation
# Dictionary containing keyword arguments for segmentation step:
parameters_segmentation={}
parameters_segmentation['target'] = 'minimum'
parameters_segmentation['method']='watershed'
parameters_segmentation['threshold']=245  # mm/h mixing ratio (until which threshold the area is taken into account)


## Tracking 
# Dictionary containing keyword arguments for the linking step:
parameters_linking={}
parameters_linking['adaptive_stop']=0.2
parameters_linking['adaptive_step']=0.95
parameters_linking['extrapolate']=0
parameters_linking['order']=1
parameters_linking['subnetwork_size']= 1000 # maximum size of subnetwork used for linking 
parameters_linking['memory']=0
parameters_linking['time_cell_min']= 12*dt 
parameters_linking['method_linking']='predict'
#parameters_linking['method_detection']='threshold'
parameters_linking['v_max']= 10
#parameters_linking['d_min']=2000
parameters_linking['d_min']=4*dxy # four times the grid spacing ?

############################################################## Tracking : Feature detection and Segmentation ###################################################################################################

# list with all files by month
#file_list= glob.glob(data_dir + '????/merg_??????.nc4')  
#print('files in dataset:  ', len(file_list))
#file_list.sort()

# test feature detection and segmentation with one monthly file 
file = data_dir + '2010/merg_201008.nc4' 
i = file[len(data_dir)+10:-4]
print('start process for file.....', i )
## DATA PREPARATION
Precip=iris.load_cube(file, 'brightness_temperature')
# FEATURE DETECTION
print('starting feature detection based on multiple thresholds')
Features=tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
print('feature detection done')

Features.to_hdf(os.path.join(savedir,'Features_' + str(i) + '.h5'),'table')
print('features saved')
    
# SEGMENTATION 
print('Starting segmentation based on surface precipitation')
Mask,Features_Precip=tobac.segmentation_2D(Features,Precip,dxy,**parameters_segmentation)
print('segmentation based on surface precipitation performed, start saving results to files')
iris.save([Mask],os.path.join(savedir,'Mask_Segmentation_precip' + str(i) + '_240_test.nc'),zlib=True,complevel=4)                
Features_Precip.to_hdf(os.path.join(savedir,'Features_Precip' + str(i) + '.h5'),'table')
print('segmentation surface precipitation performed and saved')


print(i)
arr = np.array(Mask.data)
print(arr[arr > 1000 ])



############################################################################## Linking Features ##############################################################################################

## Recombindation of feature dataframes (update framenumbers)

# # read in HDF5 files with saved features
# file_list = glob.glob(savedir  + '/Features_Precip??????.h5')  
# file_list.sort()
# print('nr. of monthly feature files:', len(file_list))


# i = 0 
# frames = 0 

# for file in file_list: 
#     if i == 0:
#         Features = pd.read_hdf(file, 'table')
#         # read in data mask with segments for tracked cells 
#         date= file[len(file)-9: len(file)-3]
#         ds = Dataset(savedir+ '/Mask_Segmentation_precip'+date+'.nc')
#         mask = np.array(ds['segmentation_mask'])  
#         # update total nr of frames 
#         frames += np.shape(mask)[0] -1
#         i = 1 
#         print('file for: ',date, 'rows: ',features.shape[0], 'frames: ', frames)

#     features = pd.read_hdf(file, 'table')
#     # update frame number and make sure they are sequential
#     features['frame'] = features['frame']  + frames
#     # append dataframes 
#     Features = Features.append(features, ignore_index=True)      
#     # read in data mask with segments for tracked cells 
#     date= file[len(file)-9: len(file)-3]
#     ds = Dataset(savedir+ '/Mask_Segmentation_precip'+date+'.nc')
#     mask = np.array(ds['segmentation_mask'])  
#     #update total nr of frames
#     frames += np.shape(mask)[0]
#     print('file for: ',date, 'rows: ',features.shape[0], 'frames: ', frames)

    
# ## Perform trjactory linking with trackpy 
# Track=tobac.linking_trackpy(Features,Precip,dt=dt,dxy=dxy,**parameters_linking)
# Track.to_hdf(os.path.join(savedir,'Tracks.h5'),'table')
