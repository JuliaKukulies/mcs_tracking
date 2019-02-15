#! /bin/bash


for year in {2006..2006}
do
    for mon in {01..04}
    do
        for day in {01..31}
        do
            for lat in {27..45}
            do
                for lon in {70..105}
                do
                    echo "searing for...." $year $lat $lon $mon $day
                    if ./CTread $year $mon $day $lat $lon 0 0 | grep -q 'DATABASE'
                    then
                        echo "CS families found!"
                        ./CTread $year $mon $day $lat $lon 0 0  >> CS_TP_${year}_${mon}.txt
                    fi
                done
            done
        done
    done
done


















