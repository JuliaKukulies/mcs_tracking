## This script takes chunks of a larger dataset (e.g. monthly netCDF files with hourly timesteps) and performs feature detection and segmentation with the cloud tracking package tobac. The output data of the features detection can be read as a pandas dataframe and when recombining all the chunks into one dataframe, trajectory linking can be performed.   

# Import necessary packages  
import iris
import numpy as np
import pandas as pd
import os,sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import datetime
import urllib,zipfile,shutil


# Import tobac itself
import tobac
from tobac import utils




#Set up directory to save output and plots:
savedir='/media/juli/Elements/gpm_v06/Save/'

## parameter

# Dictionary containing keyword options (could also be directly given to the function)
parameters_features={}
parameters_features['position_threshold']='weighted_diff' # diff between specific value and threshold for weighting when finding the center location (instead of just mean lon/lat)
parameters_features['min_num']=3 #? 
parameters_features['min_distance']=0 # minimum distance between features 

parameters_features['sigma_threshold']=0.5 # for slightly smoothing (gaussian filter)
parameters_features['n_erosion_threshold']=0 # pixel erosion (for more robust results)

parameters_features['threshold']=[5, 10  ] #mm/h, step-wise threshold for feature detection 
parameters_features['n_min_threshold']=10 # minimum nr of contiguous pixels for thresholds, 10 pixels = ca. 2000 km2




# Dictionary containing keyword arguments for segmentation step:
parameters_segmentation={}
parameters_segmentation['method']='watershed'
parameters_segmentation['threshold']=1  # mm/h mixing ratio (until which threshold the area is taken into account)


# get list with all files by month
import glob
file_list= glob.glob('/media/juli/Elements/gpm_v06/????/gpm_imerg*monthly*.nc4')  
print('files in dataset:  ', len(file_list))
file_list.sort()


dxy= 14126.0
dt = 1800 

# loop through files for features detection and creation of segmentation mask 
for file in file_list[34::]:
    i = file[44:50]
    print('start process for file.....', file)
    
    # read in data 
    Precip=iris.load_cube(file,'precipitationCal')
    
    # feature detection 
    print('starting feature detection based on multiple thresholds')
    Features=tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
    print('feature detection done')
    Features.to_hdf(os.path.join(savedir,'2000_2019/Features' + str(i) + '.h5'),'table')
    print('features saved')
    
    # segmentation 
    print('Starting segmentation based on surface precipitation')
    Mask,Features_Precip=tobac.segmentation_2D(Features,Precip,dxy,**parameters_segmentation)
    print('segmentation based on surface precipitation performed, start saving results to files')
    iris.save([Mask],os.path.join(savedir,'2000_2019/Mask_Segmentation_precip' + str(i) + '.nc'),zlib=True,complevel=4)                
    Features_Precip.to_hdf(os.path.join(savedir,'2000_2019/Features_Precip' + str(i) + '.h5'),'table')
    print('segmentation surface precipitation performed and saved')
