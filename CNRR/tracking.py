# Import a range of python libraries used in this notebook:
import glob 
import iris
import numpy as np
import pandas as pd
import os,sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import urllib,zipfile
import datetime
from netCDF4 import Dataset
import tobac



# import precip data
path='/media/juli/Data/third_pole/CNRR/data/'
file= 'cnrr_cnrr_199701_TP.nc4.nc4'
file = 'substring.nc4'
ds = Dataset(path + file)
# read in lats and lons information to add to the iris cubes 
precip = np.array(ds['prcp'])
lats = np.array(ds['LAT'])
lons = np.array(ds['LON'])

Precip = iris.load_cube(path+file, "prcp")

import iris.coords
## ad 2d lats and lons with aux coords 
aux_lats = iris.coords.AuxCoord(lats, standard_name='latitude', long_name='latitude')
aux_lons = iris.coords.AuxCoord(lons, standard_name='longitude', long_name='longitude')
Precip.add_aux_coord(aux_lats, [1,2])
Precip.add_aux_coord(aux_lons, [1,2])



dxy,dt=tobac.get_spacings(Precip, grid_spacing = 18000.0) # grid spacing needs to be given as input! 

# Dictionary containing keyword arguments for the linking step:
parameters_linking={}

parameters_linking['adaptive_stop']=0.2
parameters_linking['adaptive_step']=0.95
parameters_linking['extrapolate']=0
parameters_linking['order']=1
parameters_linking['subnetwork_size']= 100 # maximum size of subnetwork used for linking 
parameters_linking['memory']=0
#parameters_linking['time_cell_min']=5*60
parameters_linking['time_cell_min']= 3*dt 
parameters_linking['method_linking']='predict'
#parameters_linking['method_detection']='threshold'
parameters_linking['v_max']= 10
#parameters_linking['d_min']=2000
parameters_linking['d_min']=2*dxy # four times the grid spacing  (! seems to be important for GPM data)

savedir='Save'
os.makedirs(savedir,exist_ok=True)
plot_dir="Plot"
os.makedirs(plot_dir,exist_ok=True)


file_list= glob.glob(savedir  + '/2006_2016/Features_CNRR_??????.h5')
file_list.sort()

## recombination of features
i = 0 
for file in file_list: 
    if i == 0:
        print(file)
        features_p = pd.read_hdf(file, 'table')
        Features = features_p
        i +=1
        end_frame =  np.max(features_p['frame'])
    else:
        print(file, end_frame)
        features = pd.read_hdf(file, 'table')
        # update frame number and make sure they are sequential! 
        features['frame'] = features['frame']  + end_frame
        # append dataframes 
        Features = Features.append(features, ignore_index=True)
        # update last number in frame 
        end_frame = np.max(features['frame'])
        i +=1 
        print(Features.shape)




## tracking

Track=tobac.linking_trackpy(Features,Precip,dt=dt,dxy=dxy,**parameters_linking)
# remove NaN tracks from feature space!!
Track = Track.loc[Track.cell > 0]
Track.to_hdf(os.path.join(savedir,'Tracks_CNRR_2006_2016_sorted.h5'),'table')

