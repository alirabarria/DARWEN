#!/bin/bash
# multipac.sh
max_processes=6
input_dir=$1
count=0

for dir in $(ls $input_dir)
do
    ./pac $input_dir/$dir &
    ((count++))

    if [ $count -ge $max_processes ]; then
        wait
        count=0
    fi
done

wait
