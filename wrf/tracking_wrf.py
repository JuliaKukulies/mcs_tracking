## This script tracks contiguous precipitation cells in the WRF 9km simulation for the Tibetan Plateau region,

from pathlib import Path
import numpy as np
import tobac 
import pandas as pd 
import iris 
import xarray as xr 
import matplotlib.pyplot as plt
# ignore warnings
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)


# import parameters for feature detection and segmentation
from parameters import parameters_features, parameters_segmentation, dxy, dt

# function to create netcdf file with segmentation mask 
def create_netcdf(Mask, precipfile, output):
    ''' 
    Args:  

    Mask: iris cube with masked features
    precipfile(str or path): file that contains hourly precipitation 
    output(str): name of output file  
    
    '''
    xprecip= xr.open_dataset(precipfile)
    data_vars= dict(segmentation_mask=(["time", "south_north", "west_east"], Mask.data),  lat= (["south_north", "west_east"], xprecip.lat.values), lon = (["south_north", "west_east"], xprecip.lon.values)  )
    coords = dict( time= xprecip.time, south_north = xprecip.south_north.values, west_east=xprecip.west_east.values)
    data= xr.Dataset(data_vars= data_vars, coords = coords)
    data.to_netcdf(output)


    
####################### Feature detection and Segmentation##################################
# set directories for monthly chunks of hourly data and output 
savedir = Path('wrf/')
datadir = Path('/media/WRFSVars/Hourly/Precip/')


# loop through monthly files 
years= np.arange(2000,2017)
months = np.arange(1,13)

for year in years:
    for mon in months:
        if mon <10:
            mon= '0' + str(mon)
        # open file with hourly precip
        fname= 'WRFOut_TP9km_HourlyP_'+ str(year)+'_' + str(mon)+ '.nc'
        print('processing file....', fname)
        f= datadir/ fname
        Precip=iris.load_cube(str(f), 'Prep')
        
        # feature detection
        Features= tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
        outname= 'Features_' + str(year)+ str(mon) +'.h5'
        Features.to_hdf(savedir/ outname ,'table')
        print('feature detection done.')
        
        # segmentation
        Precip=iris.load_cube(str(f), 'Prep')
        Mask, Features_cells = tobac.segmentation_2D(Features, Precip, dxy, **parameters_segmentation)
        print('segmentation done.')
        mask_out = 'Mask_segmentation_' + str(year)+ str(mon) +'.nc'
        create_netcdf(Mask, datadir/ fname, savedir / mask_out)
        #iris.save([Mask], savedir/ mask_out ,zlib=True,complevel=4)
        features_out = 'Features_cells_' + str(year)+ str(mon) +'.h5'
        Features_cells.to_hdf(savedir/ features_out ,'table')
        print('files saved.')
