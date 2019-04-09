### This program identifies and tracks mesoscale convective systems based on surface precipitation. The dataset used is GPM IMERG, which contains hourly rain rates in mm/hr for 30 min timesteps,
# at a 0.1 x 0.1 degree resolution. Mesoscale convective systems are identified based on a rain rate and area threshold value and tracked on threshold values for overlapping area and persistence through
# contiguous timesteps. The output data is stored in 2D monthly netCDF4 files, containing the mean lon and lat values and precipitation statistics for each system, whereby each system has a unique ID value.




###################################### Functions #################################################################################################

import numpy as np
import os
import glob 
import random 

import matplotlib.pyplot as plt


import scipy
from scipy import ndimage
from scipy.ndimage import label, generate_binary_structure
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from collections import Counter 

import pandas as pd




# This function creates a dictionary containing all hourly files within one month, keys are month for year and values are the corresponding files 


def create_dic(working_dir):
    files={}
    keys=[]
    values=[]
    
    for year in np.arange(2014,2019,1):
        for month in np.arange(1,13,1):
            if month < 10:
                keys.append(str(year)+'0'+str(month))
            else:
                keys.append(str(year)+str(month))       

    for k in keys:
        values.append(glob.glob(working_dir + '3B-HHR.MS.MRG.3IMERG.' +  str(k) + '*.nc4'))

    ## populate dictionary with keys and values 
    files= dict(zip(keys, values))
    
    return files 



# This function reads in netcdf files and returns numpy arrays with precipitation for the specific time step and respective lon and lat grids
# returns also strings with date and time for respective timestep 

def read_in_netcdf(file, filename):
    date= filename[21:29]
    time= filename[31:35]
    dataset = Dataset(file)
    #print(dataset.variables.keys())
    
    time_slot= np.array(dataset["precipitationCal"])
    lon= np.array(dataset["lon"])
    lat= np.array(dataset["lat"])
    # fill lat and lon values over entire grid 
    lons= np.repeat(np.expand_dims(lon, axis= 1), np.shape(lat)[0], axis= 1 )
    lats= np.repeat(np.expand_dims(lat, axis= 0), np.shape(lon)[0], axis= 0)
    dataset.close()

    return time_slot, lons, lats, date, time
    


## This function assigns new labels to the identified MCS if labels already exist in label set (to avoid double naming)

def assign_labels(mcs_labels,all_mcs_labels):
    unique_labels = np.unique(mcs_labels[mcs_labels > 0 ])
    for label in unique_labels:  # loop through unique labels 
        if label in all_mcs_labels:
            new_label = range(1,5000)  # generate new identification nr if value already exists in set 
            for new in new_label:
                if new not in all_mcs_labels:
                    break
                while new in all_mcs_labels:
                    new += 1 
            mcs_labels[mcs_labels== label] = new  # assign unique values to assigned label which already exist         
            
    return mcs_labels, all_mcs_labels



## This function updates the label set

def update_label_set(mcs_labels, all_mcs_labels):
    unique_labels = np.unique(mcs_labels[mcs_labels > 0 ])
    for val in unique_labels:
        all_mcs_labels.add(val)
        
    return all_mcs_labels



## This function identifies MCS in one time slot based on a threshold rain rate and a threshold value for contiguous area
# returns following arrays for each time slot:
#  mcs : absolute rain rates 
#  mcs_labels: numbers assigned for each MCS 
# number_of_mcs: scalar containing the total number of detected mcs 

def mcs_identification(time_slot, threshold_prec, threshold_area, s):
    
    prec_loc= np.where(time_slot > threshold_prec)
    ind_row= prec_loc[0]
    ind_col= prec_loc[1]

    im= time_slot
    im[ im < 7 ]=0
    potential_mcs, number_mcs = ndimage.label(im, structure = s) # array with nr. labels of contigous pixels above threshold and nr. of total identified MCS 
    
    x= potential_mcs[potential_mcs > 0 ]
    unique, counts = np.unique(x, return_counts=True)
    labels= np.asarray((unique, counts))[0]
    selection= np.asarray((unique, counts))[1] # np.array which contains all the assigned labels for areas which fulfill intensity threshold
    large = labels[selection > threshold_area]

    # create mask for pixel areas which fulfill area threshold 
    mask= np.isin(potential_mcs,large )
    # set all pixels which do not fulfill criteria to 0 in label matrix and rain rate matrix 
    potential_mcs[mask == False ]= 0
    time_slot[mask== False ]= 0 
    mcs = time_slot
    mcs_labels= potential_mcs


    # updated number of detected MCS 
    number_mcs = np.shape(np.unique(potential_mcs))[0]-1
    
    return mcs, mcs_labels, number_mcs


## This function compares identified MCS in the next step with MCS identified in the previous timestep
# and tracks the movement based on an overlap criterium 
# returns mcs_labels_next: updated array containing unique MCS labels (with same label for those pixel groups which have been identified belonging to the same system )


def update_labels(mcs_labels, mcs_labels_next, threshold_overlap): 
    for val in np.unique(mcs_labels_next[mcs_labels_next > 0 ]):
        overlap = 0 
        loc= np.where(mcs_labels_next == val) 
        for i, x in enumerate(loc[0]):
            y= loc[1][i]
            if mcs_labels[x,y] > 0:
                overlap += 1 
                old_val = mcs_labels[x,y]
        if overlap > threshold_overlap: 
            mcs_labels_next[mcs_labels_next == val] = old_val # assign old MCS label if group contains overlap    
    return mcs_labels_next




## This function calculates and stores the lon and lat values of the system centers defined as the mean lon/lat of all
# pixels which belong to one identified system 
# returns dictionary with the number tags for each identied MCS as key values and corresponding lat and lon values 
    
def store_statistics(mcs_labels, date, time, system_stats, lats, lons, mcs):
    for i in np.unique(mcs_labels[mcs_labels > 0]):
        area= 0 
        loc= np.where(mcs_labels == i)
        for idx,x in enumerate(loc[0]):
            y = loc[1][idx]             
            R = 6371
            lat1= np.deg2rad(lats[x,y]-0.05)
            lat2 = np.deg2rad(lats[x,y]+0.05)
            lon1= np.deg2rad(lons[x,y]-0.05)
            lon2 = np.deg2rad(lons[x,y] + 0.05)
            A= (np.sin(lat2) - np.sin(lat1)) * (lon2 - lon1) * R**2
            area+= A 
        skew = scipy.stats.skew(mcs[mcs_labels ==i], axis=0, bias=True) 
        data = [str(i), str(date), str(time), np.mean(lats[mcs_labels == i ]) , np.mean(lons[mcs_labels == i ]), np.nanmean(mcs[ mcs_labels == i ]), np.nanmax(mcs[ mcs_labels == i ]), np.nanmin(mcs[ mcs_labels == i ]), area, skew] 
        system_stats.loc[len(system_stats)] = data
        system_stats.ID = system_stats.ID.astype(int)   
        system_stats.date = system_stats.date.astype(int)
        system_stats.time = system_stats.time.astype(int)
        
    return system_stats


# This function removes all identified MCS which only persist for less than the defined time threshold value 

def timestep_con(all_mcs_labels, system_stats):
    for l in all_mcs_labels:
        mcs = system_stats.loc[system_stats['ID'] == l]  
        if mcs.shape[0] < threshold_timesteps:
            system_stats = system_stats[system_stats.ID != l]
    # sort values 
    system_stats.sort_values(['ID', 'date', 'time'], ascending=True )       
    return system_stats 




# This function removes all identified MCS which do not contain at least one step where a size larger than 10 000 km2 is reached 

def size_con(system_stats):
    for l in all_mcs_labels:
        mcs = system_stats.loc[system_stats['ID'] == l]  
        if mcs.loc[mcs['area'] > threshold_max_area].shape[0] == 0:
            system_stats = system_stats[system_stats.ID != l]        
    # sort values 
    system_stats= system_stats.sort_values(['ID', 'date', 'time'], ascending=True )       
    return system_stats



## This function saves the creates pandas dataframe with all detected MCS tracks within one month to a netcdf file 

def create_netcdf(system_stats, output_path):
    data_as_xr= system_stats.to_xarray()
    data_as_xr.to_netcdf(output_path, mode = 'w', format='netCDF4', unlimited_dims=['ID']) 
    return data_as_xr 




#### This function plots all detected MCS inzoomad 


def plot_mcs(lons, lats, mcs, mcs_labels, date, time):
    #if val not in np.array(all_mcs_labels):
     #   max_id = mcs[mcs == np.max(mcs)]
      #  val = mcs_labels[mcs == max_id][0]   # get ID from MCS with strongest precipitation in this time slot 

        for val in np.unique(mcs_labels[mcs_labels > 0]):

            loc= np.where(mcs_labels == val)

            fig = plt.figure(figsize=(20, 10))

            if np.min(loc[0]) > 9:
                r1= np.min(loc[0])-10
            else:
                r1= 0 

            if np.min(loc[1]) > 9:
                c1= np.min(loc[1])-10
            else:
                c1 = 0
                
            if np.max(loc[0])+10 > np.shape(lons)[0]:
                r2= np.max(loc[0])
            else:
                r2= np.max(loc[0])+10

            if np.max(loc[1])+10 > np.shape(lons)[1]:
                c2= np.max(loc[1])
            else:
                c2= np.max(loc[1])+10                       
                
            cmap = plt.cm.get_cmap('viridis', lut= 10)

            plt.xlabel('Lon $^\circ$N')
            plt.ylabel('Lat $^\circ$N')

            m = Basemap(projection='cyl', llcrnrlat=lats[0,c1],urcrnrlat=lats[0, c2-1], llcrnrlon= lons[r1,0], urcrnrlon=lons[r2-1,0],  resolution = 'c')
            lon, lat =np.meshgrid(lons[r1:r2,0], lats[0,c1:c2])
            xi,yi = m(lon,lat)
            cs = plt.contourf(xi,yi, np.fliplr(mcs[r1:r2, c1:c2].T), cmap=cmap)
            cmap.set_under(color='lightyellow')

            cbar = plt.colorbar(extend= 'max')
            cbar.set_label(' Rain rate (mm/hr)')
            plt.rcParams.update({'font.size': 25})

            plt.savefig(working_dir + 'tracks/plots/mcs' + str(date) + str(time) + '_' + str(val) + '.png')
            plt.close()


#### This function plots the GPM derived rain rates of a specific 30 min time slot 


def plot_gpm(lons,lats, prec, date, time ):
    plt.figure(figsize=(20, 10))

    cmap = plt.cm.get_cmap('plasma')
    bounds= np.array([0, 0.5, 1 , 2, 3, 5, 7, 10])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors= 256)


    m = Basemap(projection='cyl', llcrnrlat=26.95,urcrnrlat=44.95, llcrnrlon=70.05, urcrnrlon=105.05,  resolution = 'c')
    lon, lat =np.meshgrid(lons[:,0], lats[0,:])
    xi,yi = m(lon,lat)
    cs = m.pcolormesh(xi,yi, np.fliplr(prec.T), cmap=cmap, norm = norm, vmin= 0.01, vmax = 10 )
    cmap.set_under(color='lightyellow')


    xlabels=[70, 80, 90, 100]
    ylabels= [ 27, 30, 35, 40]

    plt.xticks([70, 80,90, 100], xlabels, fontsize=25)
    plt.yticks([27,30, 35, 40],ylabels, fontsize=25)
    plt.xlabel('Lon $^\circ$N')
    plt.ylabel('Lat $^\circ$N')

    # Plot TP boundary polyline from shapefile 
    shapefile='/media/juli/Data/master_thesis/Master_thesis/data/DBATP/DBATP'
    TP_bound=m.readshapefile(shapefile, 'boundary', color='black', linewidth=2.5)


    cbar = plt.colorbar(extend= 'max')
    cbar.set_label(' Rain rate (mm/hr)')
    cbar.set_ticks(bounds)
    labels = ['0', '0.5', '1', '2', '3', '5', '7', '10']
    cbar.set_ticklabels(labels)

    plt.rcParams.update({'font.size': 25})

    plt.savefig(working_dir + 'tracks/plots/gpm_'+ str(date) + str(time) +  '.png')
    plt.close()

