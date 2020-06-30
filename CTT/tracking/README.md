## This directory contains tracking files for precipitation - brightness temperature - combined tracking os meso-scale convective systems in the Tibetan plateau region.



# Data sets

- GPM v06, 2000 - 2019, 0.1 deg and 30 min reoslution 
- NCEP 2000 -2019, 4 km and 30 min resolution 

# Requirements:

- tobac (Heikenfeld et al., 2019)


# Tracking procedure

The tracking has been performed in three steps.

1. Identification of minima in smoothed brightness temperature field
2. Trajectory linking of detected features
3. Filtering out tracks, which contain a cold and heavy rain core during their lifetime. For this an overlay with the GPM precipitation data is necessary.  