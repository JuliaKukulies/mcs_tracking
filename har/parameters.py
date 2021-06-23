# This script contains the various parameters that can be changed for the feature detection, segmentation and trajectory linking.

# temporal and spatial resolution (in seconds and meter)
dt= 3600 # hourly
dxy = 9000 # 9km



### Feature detection ###
parameters_features={}
parameters_features['position_threshold']='weighted_diff' # diff between specific value and threshold for weighting when finding the center location (instead of just mean lon/lat)
parameters_features['min_distance']=0 # minimum distance between features 
parameters_features['sigma_threshold']=0.5 # for slightly smoothing (gaussian filter)
parameters_features['n_erosion_threshold']=0 # pixel erosion (for more robust results)
parameters_features['threshold']=[3,4,5,6,7,8,9,10] #step-wise threshold for feature detection 
parameters_features['n_min_threshold']= 10 # minimum nr of contiguous pixels for thresholds
parameters_features['target']= 'maximum'

### Segmentation ###
parameters_segmentation={}
parameters_segmentation['target'] = 'maximum'
parameters_segmentation['method']='watershed'
parameters_segmentation['threshold']= 1  # threshold value until which the area is taken into account (detected feature is diluted with neighboring cells below this threshold)

### Trajectory linking ### 
parameters_linking={}
parameters_linking['adaptive_stop']=0.2
parameters_linking['adaptive_step']=0.95
parameters_linking['extrapolate']=0
parameters_linking['order']=1
parameters_linking['subnetwork_size']= 1000 # maximum size of subnetwork used for linking 
parameters_linking['memory']=0
parameters_linking['time_cell_min']= 3*dt 
parameters_linking['method_linking']='predict'
#parameters_linking['method_detection']='threshold'
parameters_linking['v_max']= 100 # maximum propagation speed allowed to connect features
parameters_linking['d_min']=4*dxy

