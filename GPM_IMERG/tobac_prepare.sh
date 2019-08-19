#! /bin/bash

# This script aggregates high-resolution GPM files in netcdf4 format and sets a time axis and calendar, so that these can be used for e.g. trackpy/ cloud tracking with tobac. For tracking applications, the files are aggregated into monthly files for each individual year containing 30 min time steps. 

for y in {2005..2019}
do
    for m in {01..12}
    do
	for d in {01..31}
	do
	    cdo cat ${y}/3B-HHR.MS.MRG.3IMERG.${y}${m}${d}*.V06B.HDF5.nc4 ${y}/gpm_imerg_${y}${m}${d}_cat.nc4 
	    cdo setcalendar,365days ${y}/gpm_imerg_${y}${m}${d}_cat.nc4 ${y}/gpm_imerg_ttaxis_${y}${m}${d}.nc4
	    cdo settaxis,${y}-${m}-${d},00:00,30min ${y}/gpm_imerg_ttaxis_${y}${m}${d}.nc4 ${y}/gpm_imerg_timeaxis_${y}${m}${d}.nc4
	    rm ${y}/gpm_imerg_ttaxis_${y}${m}${d}.nc4
	done
	cdo cat ${y}/gpm_imerg_timeaxis_${y}${m}??.nc4 ${y}/gpm_imerg_${y}${m}_tobac_input.nc4
	rm ${y}/gpm_imerg_timeaxis_${y}${m}??.nc4
	echo inputfile for ${y}${m} created.
    done
done



