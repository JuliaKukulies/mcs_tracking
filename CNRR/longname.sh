# This script changes the long_name attribute of the precipitaiton variable in all files of the CNRR datasets, in order to be able to read in the data as an iris cube with tobac

for file in *mon*nc4
do
    ncatted -a long_name,prcp,o,c,"prcp" ${file} cnrr_${file}.nc4
    rm ${file}
done

