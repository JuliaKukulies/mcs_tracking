### This program identifies and tracks mesoscale convective systems based on surface precipitation. The dataset used is GPM IMERG, which contains hourly rain rates in mm/hr for 30 min timesteps,
# at a 0.1 x 0.1 degree resolution. Mesoscale convective systems are identified based on a rain rate and area threshold value and tracked on threshold values for overlapping area and persistence through
# contiguous timesteps. The output data is stored in 2D monthly netCDF4 files, containing the mean lon and lat values and precipitation statistics for each system, whereby each system has a unique ID value.



######################## Import functions ##########################################################################################


from mcs_functions import *


############################################################ Main program ###########################################################

## Initialization of global variables

i = 0

# thresholds 
threshold_prec= 7 # rain rate mm/hr 
threshold_area= 20 # number of contigous pixels
threshold_timesteps= 6 # 6 consecutive timesteps = 3hr hr number of contiguous timesteps for which the MCS is identified 
threshold_overlap = 1 # number of pixels for overlap 
s = generate_binary_structure(2,2) # structure element which defined what type of connections are allowed in cluster finding,here: diagonal connections 
threshold_max_area = 10000


# create empty pandas dataframe 
stats= ['ID', 'date', 'time', 'lon', 'lat', 'PREC_mean','PREC_tot' ,'PREC_max', 'PREC_min', 'PREC_area', 'PREC_skew', 'PF_mean', 'PF_area', 'PF_tot']
system_stats = pd.DataFrame(columns=stats)

# set working directory 
working_dir= '/media/juli/Elements/GPM_IMERG_F_v05/GPM_finalrun/' 


# set with all mcs labels 
all_mcs_labels = set({})

# create dictionary with all files to be processes
files = create_dic(working_dir) 


################################################################ processing files ###########################################################################################


## loop through all hourly files within one month

def __processing_files__(threshold_prec, threshold_area, threshold_timesteps,threshold_overlap, stats, system_stats, working_dir, all_mcs_labels, files, i):
    for month in files.keys():

        if len(files[month]) == 0:
            print('no files for month  ', month)

        else:  
            #i = 0 
            for file in files[month]:
                if i == 0: 
                    filename = file[50::]
                    # read in first netcdf file
                    time_slot, prec,  lons, lats, date, time = read_in_netcdf(file, filename)
                    time_slot= extract_high_elevations(time_slot)
                    # plot over time slot 
                    #plot_gpm(lons, lats, prec, date, time)
                    # identify MCS in first netcdf file 
                    mcs, mcs_labels, number_mcs= mcs_identification(time_slot,threshold_prec,threshold_area, s)
                    # add labels to label set and control that assigned labels are not already in MCS label set 
                    mcs_labels, all_mcs_labels = assign_labels(mcs_labels, all_mcs_labels )
                    # update label set 
                    all_mcs_labels= update_label_set(mcs_labels, all_mcs_labels)
                    # update MCS statistics 
                    system_stats = store_statistics(mcs_labels, date, time, system_stats, lats, lons, mcs, prec, s)
                    # save plots of detected MCS
                    #plot_mcs(lons, lats, mcs, mcs_labels, date, time)
                    #print('first file:', filename, system_stats.shape)

                while i < np.shape(files[month])[0]-1: 
                    # read in next timestep 
                    file_next= files[month][i+1]
                    filename_next = file_next[50::]
                    print('reading in.....', filename_next)
                    time_slot_next, prec_next, lons, lats, date, time = read_in_netcdf(file_next, filename_next)
                    time_slot_next= extract_high_elevations(time_slot_next)
                    # plot over time slot 
                    #plot_gpm(lons, lats, prec_next, date, time)
                    # identify MCS in next timestep 
                    mcs_next, mcs_labels_next, number_mcs_next = mcs_identification(time_slot_next,threshold_prec,threshold_area, s)
                    if np.shape(mcs_labels_next[mcs_labels_next > 0])[0] > 0: # if MCS are present in new timestep
                        # add labels to label set and control that assigned labels are not already in MCS label set 
                        mcs_labels_next, all_mcs_labels = assign_labels(mcs_labels_next, all_mcs_labels )
                        if np.shape(mcs_labels[mcs_labels > 0])[0] > 0: # if MCS were present in previous timestep 
                            # compare MCS to previous timestep and track systems with overlap (through assigning same label)
                            mcs_labels_next = update_labels(mcs_labels, mcs_labels_next, threshold_overlap)
                            # update label set 
                        all_mcs_labels= update_label_set(mcs_labels_next, all_mcs_labels)
                        # update MCS statistics 
                        system_stats = store_statistics(mcs_labels_next, date, time, system_stats, lats, lons, mcs_next, prec_next, s)
                        # save plots of detected MCS
                        #plot_mcs(lons, lats, mcs_next, mcs_labels_next, date, time)
                    else:
                        print('file does not contain any MCS ')

                    # control that next file is openenend and open labels become the MCS labels from previous timestep to compare with 
                    i += 1 
                    mcs_labels = mcs_labels_next 

            system_stats = system_stats.set_index('ID', drop = False )
            system_stats = system_stats.sort_values(['ID', 'date', 'time'] )
            system_stats = timestep_con(all_mcs_labels, system_stats, threshold_timesteps)
            #system_stats = size_con(all_mcs_labels, system_stats)
            output_path= working_dir + '/tracks/tracking_7mm6h20elev'+   date[0:6]+ '_' + 'tracks.nc'
            create_netcdf(system_stats, output_path)    





# process files with elevation extraction
__processing_files__(threshold_prec, threshold_area, threshold_timesteps,threshold_overlap, stats, system_stats, working_dir, all_mcs_labels, files, i)


 # process files without elevation extraction 

