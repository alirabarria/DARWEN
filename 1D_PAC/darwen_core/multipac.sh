#!/bin/bash
# multipac.sh
input_dir=$1
max_processes=$2
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
