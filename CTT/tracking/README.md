
# MCS tracking for the Tibetan Plateau region 

This directory contains scripts for tracking Meso-scale convective systems in the Tibetan Plateau region. 


## Data sets

The tracking has been performed using collocations of infrared brightness temperatures (NCEP/CPC) and satellite precipitation estimates from GPM IMERG. 

- GPM v06, 2000 - 2019, 0.1 deg and 30 min resolution: https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_06/summary?keywords=GPM%20IMERG%20v06 
- NCEP/CPC 2000 -2019, 4 km and 30 min resolution: https://disc.gsfc.nasa.gov/datasets/GPM_MERGIR_1/summary

## Requirements:

- The python package tobac: https://github.com/climate-processes/tobac


## MCS dataset 

The resulting MCS dataset for 2000 - 2019 can be downloaded from 

https://zenodo.org/record/4767152#.YLD_cXUzbmE 


## Tracking procedure

The tracking has been performed in three steps.

1. Identification of minima in smoothed brightness temperature field (Feature detection)
2. Trajectory linking of detected features based on their average propagation speed and location 
3. Filtering out tracks, which that do not contain a cold core or heavy rain area during their lifetime.  


The detailed method and analysis of MCS tracks are described in this paper: https://www.essoar.org/doi/10.1002/essoar.10504239.1 


## Data Pre-processing 

The feature detection was performed on monthly aggregated files of half-hourly data. The trajectory linking was performeed by linking the monthly features in one year together. 



*Contact: julia.kukulies@gu.se*
