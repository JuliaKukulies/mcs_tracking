#! /bin/bash


## This script takes the lat and lon information from a different NC file and adds it to a dataset

for file in *mon*nc4
do
    ncks -A -v LAT surface_20150112.nc ${file}
    ncks -A -v LON surface_20150112.nc ${file}
done

