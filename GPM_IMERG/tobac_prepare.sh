#! /bin/bash

# This script sets a time axis and calendar for modified GPM files, so that these can be used for e.g. cloud tracking in tobac 
## It also creates merged files for each individual year and month containing 30 min time steps (to allow for chunking in tobac)

for y in {2017..2018}
do
    for m in {01..12}
    do
	for d in {01..31}
	do
	    cdo setcalendar,365days GPM_IMERG_${y}${m}${d}_cat.nc4 gpm_imerg_ttaxis_${y}${m}${d}.nc4
	    cdo settaxis,${y}-${m}-${d},00:00,30min gpm_imerg_ttaxis_${y}${m}${d}.nc4 gpm_imerg_timeaxis_${y}${m}${d}.nc4
	    rm gpm_imerg_ttaxis_${y}${m}${d}.nc4
	done
	cdo cat gpm_imerg_timeaxis_${y}${m}??.nc4 gpm_imerg_${y}${m}_tobac_input.nc4
    done
done



