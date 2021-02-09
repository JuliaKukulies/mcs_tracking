### Import packages 

from netCDF4 import Dataset
import numpy as np
import xarray as xr 
import pandas as pd 
import glob
import os
import warnings
from scipy import ndimage
from scipy.ndimage import generate_binary_structure

# for image labeling 
s = generate_binary_structure(2,2)
# ignore warnings 
warnings.filterwarnings('ignore')


## Precip files 
precip_files = glob.glob('/media/juli/Elements/gpm_v06/200[1-9]/gpm_imerg_??????_monthly.nc4')
for i in glob.glob('/media/juli/Elements/gpm_v06/201[0-9]/gpm_imerg_??????_monthly.nc4'):
    precip_files.append(i)    
precip_files.sort()


## loop through monthly precip files
for f in precip_files[0::]:
    congrid = np.zeros((600,350))
    year = int(f[44:48])
    month = int(f[48:50])

    print(year, month)
    if len(str(month)) ==1:
        month = '0'+ str(month)

    if int(month) in [6,7,8]:
    
        # open precip file
        prec = xr.open_dataarray(f)
        
        # open corresponding mask file
        fi = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs/Mask_Segmentation_'+str(year) + str(month) +'.nc'
        # open file with tracks

        
        trackfile = '/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs/tracks/Tracks_'+ str(year)+'_heavyraincore3mm.h5'
        tracks = pd.read_hdf(trackfile, 'table')
        # get for specific month 
        subtracks = tracks[tracks.time.dt.month == int(month)]

        if os.path.isfile(fi) == True:
            mcsdata = xr.open_dataarray(fi)

            if np.shape(prec)[0] == np.shape(mcsdata)[0]:
                p = prec[:,:-1,:-1]
                pr = np.array(p.data)
                # extreme precip
                pr[pr< 2.16] = np.nan

                # feature mask to get the fields from tracked MCS 

                # loop through timesteps
                for t in np.arange(0,np.shape(prec)[0]):
                    i = 0 

                    # get precip and mask for specific timestep 
                    prect= p[t,:,:]
                    # remove very low values  
                    prect.data[prect.data < 2.16] = np.nan

                    mcst = mcsdata[t, :, :].T
                    features = mcst.data
                    # get new labels so that connected features have only one ID 
                    mcslabels, nr = ndimage.label(features, structure = s)

                    # get labels for connected precipitation features as well
                    prcplabels, nr = ndimage.label(prect.data, structure = s)
                    prcp_features = prect.copy()
                    prcp_features.data = prcplabels

                    # set mask values to zero, where no feature in table, otherwise get whole segmented object for feature 
                    trackfeatures = subtracks[subtracks.idx == t].feature.values
                    # loop through unique mask cloud objects
                    for l in np.unique(mcslabels[mcslabels > 0]):
                        # get features in mask 
                        tbbfeatures = np.unique(features[mcslabels == l])
                        # compare features in mask to tracked features and change l area to zero if feature is not tracked 
                        if any(x in tbbfeatures for x in trackfeatures) == False:
                            mcslabels[mcslabels == l] = 0


                    if np.shape(mcslabels[mcslabels>0])[0] > 0:
                        # create MCS mask 
                        mcsmask = mcst.where(mcslabels > 0)
                        mcsmask.coords['mcsmask'] = (('lon', 'lat'), mcsmask)

                        ## add mcs-associated precip of month to empty grid
                        feature_labels = prcplabels[mcslabels > 0 ]
                        #feature_labels  = np.unique(prcp_features.where(mcsmask.coords['mcsmask'].values > 0 ).values

                        for pf in feature_labels[feature_labels> 0]:
                            prect.data[prcplabels != pf] = 0

                        # add all contiguous precip features until 0.1 mm/hr 
                        stacked = np.dstack((congrid, prect.data))
                        congrid= np.nansum(stacked, axis = 2 )

                # calculate total sum of extreme precip in month (along time axis)
                extreme_precip = np.nansum(pr, axis = 0 )

            #################### write and save netcdf file##########################################

            # Creating dimensions
            data = Dataset('/media/juli/Data/projects/data/satellite_data/ncep/ctt/Save/tcs/mcs_contr_precip'+str(year)+ str(month) + '_extreme.nc4','w', format = 'NETCDF4_CLASSIC')

            lat =data.createDimension('lat',np.shape(p)[2])
            lon =data.createDimension('lon',np.shape(p)[1])    

            # Creating variables
            lats = data.createVariable('latitude',np.float64, ('lat',))
            lons = data.createVariable('longitude',np.float64, ('lon',))
            precip = data.createVariable('precip',np.float64, ('lon','lat'))
            mcs = data.createVariable('mcs',np.float64, ('lon','lat'))

            # Creating attributes

            lats.units = 'degrees_north'
            lons.units = 'degrees_east'

            # Write data to variable

            latitude=prect.lat.values
            longitude=prect.lon.values
            mcs[:]= congrid
            precip[:] = extreme_precip
            lats[:]=latitude
            lons[:]=longitude

            # Close file
            data.close()
            # close precip dataset
            prec.close()
            print('netcdf file created.', str(year), str(month))



