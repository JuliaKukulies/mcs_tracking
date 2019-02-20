#! /bin/bash



indices=($(seq 0 82))

stations=(51886 51931 52602 52633 52645 52657 52707 52713 52737 52754 52765 52818 52825 52836 52856 52866 52868 52908
	  52943 52955 52974 55228 55248 55279 55294 55299 55437 55472 55493 55569 55578 55585 55591 55655 55664 55680
	  55690 55696 55773 56004 56018 56021 56029 56033 56034 56038 56043 56046 56065 56067 56074 56079 56080 56106
	  56116 56125 56137 56144 56146 56151 56152 56167 56172 56173 56178 56182 56202 56223 56227 56247 56251 56257
	  56312 56331 56357 56374 56434 56444 56459 56462 56533 56543 56548)

lons=(90.85 81.65 93.33 98.42 99.6 100.25 93.68 95.35 97.38 100.13 101.62 94.92 96.43 98.1 100.62 101.75 101.37 93.08 99.98 100.73 102.03 80.08 84.42 90.02 91.1 92.07 81.25 88.63 91.1 87.6 88.88 90.17 91.13 85.97 87.08 89.6 91.95 92.47 89.08 92.43 95.28 95.8 96.97 98.22 97.13 98.1 100.23 99.65 101.6 101.48 102.08 102.97 102.9 93.78 95.6 96.47 97.17 98.58 100 100.75 100.33 101.12 102.23 102.55 102.35 103.6 93.28 95.83 95.77 99.1 100.32 100.27 94.33 97.83 100.3 101.97 97.47 98.92 101.27 101.5 98.67 99.75 99.28)

lats=(38.25 36.85 38.75 38.82 38.43 38.18 36.8 37.85 37.37 37.33 37.38 36.42 36.43 36.3 36.27 36.73 36.02 35.22 35.58 35.58 35.55 32.5 32.15 31.38 32.35 31.48 30.28 30.95 30.48 29.08 29.25 29.43 29.67 28.18 28.63 28.92 27.98 28.42 27.73 34.22 32.88 34.12 33 34.92 33.8 32.98 34.48 33.75 34.73 33.43 34 33.58 35 31.88 31.42 32.2 31.15 31.8 31.62 32.93 32.28 30.98 31.9 32.8 31 32.67 30.67 30.75 29.87 30 30.93 30 29.67 29.67 29.05 30.05 28.65 28.48 27.93 29 27.75 27.85 27.17)



for i in ${indices[*]}
do
    #echo ${lats[$i]} ${lons[$i]} ${stations[$i]}.txt

    ncks --no_nm_prn -H -C -v precipitationCal -d lat,${lats[$i]} -d lon,${lons[$i]} GPM_IMERG_2014_JJA.nc4 > ${stations[$i]}.txt
done 



