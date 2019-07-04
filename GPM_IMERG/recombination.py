## This script reads in h5 tables with detected features and recombines them into one dataframe which then can be used as input for the trajectory linking in tobac.



# read in data
file = savedir + '/Features_Precip201410.h5'
features_p = pd.read_hdf(file, 'table')
np.max(features_p['frame'])



# recombination of dataframes with update of frame number 
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


