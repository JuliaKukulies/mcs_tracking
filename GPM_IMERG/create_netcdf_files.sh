#! /bin/bash

## This bash script creates .ctl description files for CMORPH binary data files (CTR) and converts the binary data files to netCDF format  




for y in {2014..2014};
do
    for m in {01..12};
    do
	for d in {01..31};
	do
	    for h in {00..23};
	    do
		out=${y}${m}${d}${h}.ctl
		echo "DSET CMORPH_V1.0_CRT_8km-30min_"${y}${m}${d}${h} >> ${out}
		echo "OPTIONS template little_endian" >> ${out}
		echo "UNDEF  -999.0" >> ${out}
		echo "TITLE  Precipitation estimates" >> ${out}
		echo "XDEF 4948 LINEAR   0.036378335 0.072756669" >> ${out}
		echo "YDEF 1649 LINEAR -59.963614    0.072771377" >> ${out}
		echo "ZDEF   01 LEVELS 1" >> ${out}
		echo "TDEF 999999 LINEAR  00z01jan1998 30mn" >> ${out}
		echo "VARS 1" >> ${out}
		echo "cmorph   1  99  cmorph [ mm/hr ]" >> ${out}
		echo "ENDVARS" >> ${out}

   cdo -f nc import_binary -i ${out} cmorph_${y}${m}${d}${h}.nc 

done


