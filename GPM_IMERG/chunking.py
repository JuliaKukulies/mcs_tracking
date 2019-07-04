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


# get list with all files by month
import glob
file_list= glob.glob('/media/juli/Data/third_pole/tobac/examples/climate-processes-tobac_example_data-b3e69ee/data/gpm/*tobac_input.nc4')


# loop through files for features detection and creation of segmentation mask 
for file in file_list:
    i = file[107:113]
    print('start process for file.....', file)
    
    # read in data 
    Precip=iris.load_cube(file,'precipitationCal')
    # define temporal and spatial resolution
    dxy,dt=tobac.get_spacings(Precip)
    
    # feature detection 
    print('starting feature detection based on multiple thresholds')
    Features=tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
    print('feature detection done')
    Features.to_hdf(os.path.join(savedir,'Features' + str(i) + '.h5'),'table')
    print('features saved')
    
    # segmentation 
    print('Starting segmentation based on surface precipitation')
    Mask,Features_Precip=tobac.segmentation_2D(Features,Precip,dxy,**parameters_segmentation)
    print('segmentation based on surface precipitation performed, start saving results to files')
    iris.save([Mask],os.path.join(savedir,'Mask_Segmentation_precip' + str(i) + '.nc'),zlib=True,complevel=4)                
    Features_Precip.to_hdf(os.path.join(savedir,'Features_Precip' + str(i) + '.h5'),'table')
    print('segmentation surface precipitation performed and saved')
