import pandas as pd
import numpy as np
import xarray as xr





## elevation mask 
import xarray as xr 
dem = '/media/juli/Data/projects/data/elevation/elevation_600x350.nc'
elevations = xr.open_dataarray(dem)
# mask as coordinates 
dem_mask = elevations.where(elevations >= 3000)
dem_mask.coords['mask'] = (('lon', 'lat'), dem_mask)


from scipy import ndimage
from scipy.ndimage import generate_binary_structure
s= generate_binary_structure(2,2)


def plateau_mask(tracks):
    removed = 0
    tracks['tp_flag'] = 0
    tracks['total_precip'] = 0 # total precip of segmented feature 
    tracks['convective_precip'] = 0 # rain from convective area 
    tracks['rain_flag'] = 0 # pixels >= 5mm/hr
    tracks['max_precip'] = 0 # pixel with max rain rate 

    # loop through cells 
    for cell in np.unique(tracks.cell.values):
        subset = tracks[tracks.cell == cell]
        tp_flag = 0 
        # loop through timesteps of features for specific cell 
        for idx in subset.idx.values: 
            # idx is the timestep index for respective timestep or mask file 

            # open corresponding precip and mask file 
            year = subset.time.values[0].year 
            month = subset.time.values[0].month
            if len(str(month))== 1: 
                month= '0' + str(month)

            # check whether segmented feature is in area above 3000 m 
            maskfile = '/media/juli/Elements/gpm_v06/Save/2000_2019/Mask_Segmentation_precip'+str(year) + str(month) + '.nc'
            mask = xr.open_dataarray(maskfile)
            mask= mask[:,1:,1:]

            # get right timestep frames 
            seg= mask[idx,:, :]

            # get write features from segmentation mask 
            labels, nr = ndimage.label(seg, structure = s)

            # get feature ID for frame                                                                                                                     
            featureid= subset.feature[subset.idx== idx].values[0]
            
            if featureid in seg:
                label = np.unique(labels[ seg == featureid])[0]
                seg_mask = seg.where(labels == label)
                # create mask as coordinates                                                                                                               
                seg_mask.coords['mask'] = (('lon', 'lat'), seg_mask)

                # Elevation mask                                                                                                                          
                elevation_values = dem_mask.where(seg_mask.coords['mask'].values > 1)
                arr= elevation_values.values.flatten()
                values = arr[~np.isnan(arr)]     

                mountain_features = values[values >=3000].shape[0]
                tracks['tp_flag'][tracks.feature == featureid] =  mountain_features

                if mountain_features == 0 : 
                    tracks = tracks.drop(tracks[tracks.cell == cell].index)
                else:
                    tracks['tp_flag'][tracks.feature == featureid] =  mountain_features

                precipfile = '/media/juli/Elements/gpm_v06/'+str(year)+'/gpm_imerg_'+ str(year)+str(month)+'_monthly.nc4'
                precip = xr.open_dataarray(precipfile)
                precip = precip[:,1:,1:].T
                if np.shape(precip)[2] == np.shape(mask)[0]:
                    prec= precip[:, :, idx]
                    # apply mask on precip data to extract precip values for feature in cell                                                                  
                    precip_values = prec.T.where(seg_mask.coords['mask'].values > 1)
                    arr= precip_values.values.flatten()
                    values = arr[~np.isnan(arr)] # values contains the amount of grid cells with precip
                    if np.shape(values[values >=5])[0] != 0:
                        total_precip = np.nansum(values[values > 0]) * 0.5
                        tracks['total_precip'][tracks.feature == featureid] = total_precip
                        rain_features = values[values >= 5].shape[0]
                        tracks['rain_flag'][tracks.feature == featureid] = rain_features
                        tracks['convective_precip'][tracks.feature == featureid] = np.nansum(values[values >= 5]) * 0.5
                        tracks['max_precip'][tracks.feature == featureid] = np.nanmax(values[values >=5])    

    return tracks 

   
years = np.arange(2000,2019)
for y in years:
    # read in precip tracks 
    f = '/media/juli/Elements/gpm_v06/Save/2000_2019/Tracks_precipitation_GPM_'+ str(y) + '.h5'
    ptracks= pd.read_hdf(f, 'table')
    # remove nan values to only save the linked features                                                                                                      
    ptracks = ptracks[ptracks.cell >= 0]
    ptracks.timestr = pd.to_datetime(ptracks.timestr)
    tracks = plateau_mask(ptracks)
    tracks.to_hdf('/media/juli/Elements/gpm_v06/Save/2000_2019/Tracks_GPM_'+ str(y)+'_TPflag.h5','table')
    print('new tracks saved for year ', y)
