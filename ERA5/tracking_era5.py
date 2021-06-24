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

dxy= 30000
parameters_features['n_min_threshold'] = 10

# function to create netcdf file with segmentation mask 
def create_netcdf(Mask, xprecip, output):
    ''' 
    Args:  

    Mask: iris cube with masked features
    xprecip(str or path): extracted xarray that contains hourly precipitation 
    output(str): name of output file  
    
    '''
    data_vars= dict(segmentation_mask=(["time", "latitude", "longitude"], Mask.data)  )
    coords = dict( time= xprecip.time, latitude = xprecip.latitude.values, longitude = xprecip.longitude.values)
    data= xr.Dataset(data_vars= data_vars, coords = coords)
    data.to_netcdf(output)

    
####################### Feature detection and Segmentation##################################
# set directories for monthly chunks of hourly data and output 
savedir = Path('era5/')
datadir = Path('/media/CDS/ERA5/Hourly/NCs/')


# loop through monthly files 
years= np.arange(2000,2017)
months = np.arange(1,13)



for year in years:
    for mon in months:
        if mon <10:
            mon= '0' + str(mon)
        # open file with hourly precip
        fname= 'ERA5_'+ str(year)+'-' + str(mon)+ '_PT_GlHrly.nc'
        print('processing file....', fname)
        f= datadir/ fname

        # extract domain and fix unit (from m to mm per hour) 
        precip= xr.open_dataset(f)
        lc= precip.coords["longitude"]
        la= precip.coords["latitude"]
        precip_tp= precip.tp.loc[dict(longitude=lc[(lc > 50) & (lc < 135)], latitude=la[(la > 10) & (la < 50)])] * 1000
        Precip= xr.DataArray.to_iris(precip_tp)

        # feature detection
        Features= tobac.feature_detection_multithreshold(Precip,dxy,**parameters_features)
        outname= 'Features_' + str(year)+ str(mon) +'.h5'
        Features.to_hdf(savedir/ outname ,'table')
        print('feature detection done.')
        
        # segmentation
        Mask, Features_cells = tobac.segmentation_2D(Features, Precip, dxy, **parameters_segmentation)
        print('segmentation done.')
        mask_out = 'Mask_segmentation_' + str(year)+ str(mon) +'.nc'
        create_netcdf(Mask, precip_tp, savedir/mask_out)
        #iris.save([Mask], savedir/ mask_out ,zlib=True,complevel=4)
        features_out = 'Features_cells_' + str(year)+ str(mon) +'.h5'
        Features_cells.to_hdf(savedir/ features_out ,'table')
        print('files saved.')
