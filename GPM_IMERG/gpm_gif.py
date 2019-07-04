import numpy as np
import os
import glob 
import random 

import matplotlib.pyplot as plt
import cartopy
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors

import scipy
from scipy import ndimage
from scipy.stats import skew
from scipy.ndimage import label, generate_binary_structure
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from collections import Counter 

import pandas as pd 


import xarray as xr
import xesmf as xe


def plot_gpm(lons,lats, prec, date, time ):
    plt.figure(figsize=(20, 10))

    cmap = plt.cm.get_cmap('plasma_r')
    bounds= np.array([ 0.1, 0.5, 1 , 2, 3, 4, 5, 6, 7])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors= 256)


    m = Basemap(projection='cyl', llcrnrlat=26.95,urcrnrlat=44.95, llcrnrlon=70.05, urcrnrlon=105.05,  resolution = 'c')
    lon, lat =np.meshgrid(lons[:,0], lats[0,:])
    xi,yi = m(lon,lat)
    cs = m.pcolormesh(xi,yi, prec.T, cmap=cmap, norm = norm, vmin= 0.1, vmax = 7 )
    cmap.set_under(color='lightyellow')

    xlabels=[70, 80, 90, 100]
    ylabels= [ 27, 30, 35, 40]

    plt.xticks([70, 80,90, 100], xlabels, fontsize=25)
    plt.yticks([27,30, 35, 40],ylabels, fontsize=25)
    plt.xlabel('Lon $^\circ$N')
    plt.ylabel('Lat $^\circ$N')
    plt.title(str(date)+ ' ' + time + 'UTC')

    # Plot TP boundary polyline from shapefile 
    shapefile='/media/juli/Data/master_thesis/Master_thesis/data/DBATP/DBATP'
    TP_bound=m.readshapefile(shapefile, 'boundary', color='black', linewidth=2.5)


    cbar = plt.colorbar(extend= 'max')
    cbar.set_label(' Rain rate (mm/hr)')
    cbar.set_ticks(bounds)
    labels = ['0.1', '0.5', '1', '2', '3', '4', '5','6',  '7']
    cbar.set_ticklabels(labels)

    plt.rcParams.update({'font.size': 25})

    plt.savefig(working_dir + 'tracks/plots/gpm_'+ str(date) + str(time) +  '.png')
    plt.close()



def read_in_netcdf(file, filename):
    date= filename[21:29]
    time= filename[31:35]
    dataset = Dataset(file)
    
    time_slot= np.array(dataset["precipitationCal"])
    prec= np.array(dataset["precipitationCal"])
    lon= np.array(dataset["lon"])
    lat= np.array(dataset["lat"])
    # fill lat and lon values over entire grid 
    lons= np.repeat(np.expand_dims(lon, axis= 1), np.shape(lat)[0], axis= 1 )
    lats= np.repeat(np.expand_dims(lat, axis= 0), np.shape(lon)[0], axis= 0)
    dataset.close()
    
    return time_slot, prec, lons, lats, date, time



def create_dic(working_dir):

    files={}
    keys=[]
    values=[]
    for day in np.arange(20,21,1):
        if day < 10:
            keys.append('0'+str(day))
        else:
            keys.append(str(day))       

    for k in keys:
        values.append(glob.glob(working_dir + '3B-HHR.MS.MRG.3IMERG.201708' +  str(k) + '*.nc4'))


    return files


###########################################################################################################################################



working_dir= '/media/juli/Elements/GPM_IMERG_F_v05/GPM_finalrun/'

files= create_dic(working_dir)


##################################################################################################################

for file in files:
    print(file)
    filename = file[50::]
    time_slot, prec,  lons, lats, date, time = read_in_netcdf(file, filename)
    gpm_plot(lons, lats, prec, date, time )


########################################################################################################################


print(files)






