
# Import libraries
import iris
import numpy as np
import pandas as pd
import os,sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import datetime
import urllib,zipfile,shutil
from netCDF4 import Dataset
import glob
# Import tobac itself
import tobac
from tobac import utils

## This python script calculates the area of al features which are associated with a tracked cell in tobac

## function to extract only area with tracked features



def update_mask(mask, feat_select):
    import scipy.ndimage as ndi

    for frame, xy in enumerate(mask):
        segments, n = ndi.label(xy)
        for i in np.unique(segments[segments > 0]):
            keep = 0 
            for f in xy[segments == i]:
                # if one of the features is in tracking dataframe, keep marked area 
                if f in feat_select.feature.values:
                    keep += 1
                else:
                    keep += 0
                    
            if keep == 0 :
                xy[segments == i] = 0
                mask[frame] = xy 
                
    return mask

data_dir = '/media/juli/Data/third_pole/mcs_tracking/CNRR/Save'

# import tracks for GPM 
file = data_dir + '/Tracks_CNRR_2006_2016_updatedframes.h5'
Tracks = pd.read_hdf(file, 'table')
Tracks['timestr']=pd.to_datetime(Tracks['timestr'],format='%Y-%m-%d %H:%M:%S')

## calculate area

mask_chunks= glob.glob(data_dir + '/2006_2016/Mask_Segmentation_precip??????.nc')  
mask_chunks.sort()


features = Tracks
for file in mask_chunks:
    date= file[len(file)-9: len(file)-3]
    yearmonth= file[len(file)-9: len(file)-5] + '-' + file[len(file)-5: len(file)-3]
    print('calculating area for...', yearmonth)
    # read in  mask as iris cube 
    mask=iris.load_cube(file,'segmentation_mask')
    mask.attributes = None


    # read in corresponding precip cells 
    precip = iris.load_cube(data_dir + '/Precip_cells'+ date +'.nc')
    precip_arr= precip.data
    print('mask precip= ', np.shape(mask.data[mask.data > 0])[0], np.shape(precip_arr[precip_arr > 0 ])[0])

    # select feature for specific month  
    #feat_select= Tracks.loc[Tracks['timestr'].dt.strftime('%Y-%m') == yearmonth]
    #mask_arr = update_mask(mask_arr, feat_select)
    mask.data = precip_arr
    print('total precip cells segment mask:', np.shape(mask.data[mask.data > 0])[0])

    features = tobac.calculate_area(features, mask, method_area= 'latlon')
    print('mean area for track features was...', np.nanmean(features.area.data))
    

## save new dataframe
features.to_hdf(data_dir + '/gpm_tracks_area.h5','table')
