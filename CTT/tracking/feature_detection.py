## This python script performs a tracking cloud-features based tobac (Heikenfeld et al., 2019) on NCEP brightness temperatures. The specific tracking here aims for mesoscale convective systems in the Tibetan Plateau region, whereby a meso-scale convective system is defined as a cloud complex with
# * contiguous area of 50 000 km2 < 221 K 
# * persistence of this area for more than three hours 
# * existence of an area < 200 K at least once during MCS lifetime 
# * existence of at least 10 % of the cloud area with rain rates >= 5mm/hr at least once during MCS lifetime 

# Note that the input files are monthly files of merged 30-min data aggregated in one file, so that the feature detection and segmentation step of the tracking can be performed per month, whereby the linking is is conducted in the way that it can even link systems at the boundaries of two months.  

# created by Julia Kukulies, julia.kukulies@gu.se 
#####################################################################################################################

import tobac
import glob
import os 
import gc
import iris
import numpy as np
import pandas as pd
import datetime
from netCDF4 import Dataset
from parameters import parameters_features, parameters_segmentation, dt, dxy, savedir, data_dir

# get list with all files by month
file_list= glob.glob(data_dir + 'merg_??????.nc4')  
print('files in dataset:  ', len(file_list))
file_list.sort()

for f in file_list:
    year= f[len(data_dir)+10:-4]
    month = f[len(data_dir)+14:-4]
    print('start process for file.....', year, month )

    ## GET DATA (Proxy for convection could be cloud top temperature, precipitation or other parameters)
    Proxy=iris.load_cube(f, 'brightness_temperature')
    # set negative values to NaN
    for i in range(Proxy.shape[0]):
        part = Proxy.data[i,:,:]
        part[part < 0] = np.nan
        Precip.data[i] = part 

    # FEATURE DETECTION
    print('starting feature detection based on multiple thresholds')
    Features=tobac.feature_detection_multithreshold(Proxy,dxy,**parameters_features)
    print('feature detection done')
    Features.to_hdf(os.path.join(savedir,'Features_' + str(year) +str(month) + '.h5'),'table')
    print('features saved', Features.shape)
    
    # SEGMENTATION 
    print('Starting segmentation based on surface precipitation')
    Mask,Features_Proxy=tobac.segmentation_2D(Features,Proxy,dxy,**parameters_segmentation)
    print('segmentation based on surface precipitation performed, start saving results to fs')
    iris.save([Mask],os.path.join(savedir,'Mask_Segmentation_' + str(year) + str(month) + '.nc'),zlib=True,complevel=4)                
    Features_Proxy.to_hdf(os.path.join(savedir,'Features_cells_' + str(year) +  str(month) + '.h5'),'table')
    print('segmentation surface precipitation performed and saved')

    del Proxy
    del Features
    del Mask
    del Features_Proxy
    gc.collect()
