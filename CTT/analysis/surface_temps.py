## This function downloads ERA5 files 

def download_era5_surfacetemp(year, month): 
    import cdsapi                                                                                           
    # Open a new Client instance                                                                            
    c = cdsapi.Client()         
    # output filename 
    file= 'surfacetemps/era5_'+'_'+ year + month+'.nc'
    # Send request (download data)                                                                     
    c.retrieve('reanalysis-era5-single-levels', {                                                           
            "product_type":   "reanalysis",                                                                 
            "format":         "netcdf",                                                                     
            "area":           "45.00/70.00/25.00/105.00",                                                   
            "variable":       '2t',                                                                                              
            "year":          [year],                                                      
            "month":         [ month],                
            "day":             ["01" , "02","03","04","05","06","07","08","09","10","11",                      
                           "12","13","14","15","16","17","18","19","20","21","22",                      
                           "23","24","25","26","27","28","29","30","31" ],                               
            "time": ['00:00','01:00','02:00',                                                                    
            '03:00','04:00','05:00',                                                                    
            '06:00','07:00','08:00',                                                                    
            '09:00','10:00','11:00',                                                                    
            '12:00','13:00','14:00',                                                                    
            '15:00','16:00','17:00',                                                                    
            '18:00','19:00','20:00',                                                                    
            '21:00','22:00','23:00'     
            ]
                                                                                                   
        }, file)                                                                

    print(file, 'downloaded and saved.')   
    return file





## import  elevation data 

import xarray as xr 
demfile = '/media/juli/Data/projects/master_thesis/Master_thesis/data/DEM_TP/dem_ERA5_format.nc'
dem= xr.open_dataarray(demfile)




#############################################################


# extract information from per month and year 
import numpy as np
import os 

years = np.arange(2000,2016).astype(str)
months = np.arange(1,13).astype(str)

# loop through month in year 
for year in years:
    for month in months:
        print('start getting ERA5 surface temperatures for', year, month)

        # download file 
        file = download_era5_surfacetemp(year, month)
        # open file 
        temps= xr.open_dataarray(file)
        # extract elevations > 3000 m 
        temps.data[:, dem.data < 3000] = 999

        # save bin counts for file 
        counts, bins = np.histogram(temps, bins = np.arange(190,320,5))
        np.savetxt('surfacetemps/counts_'+ str(year)+ str(month)+'.txt', counts)

        # remove file
        os.remove(file)
        print('info extracted from ', file)






