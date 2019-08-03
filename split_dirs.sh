#!/bin/bash

if [ $# -eq 0 ]; then
    echo "usage: split_dirs.sh dir_size dir_prefix type [f or l]"
    exit 1
fi

dir_size=${1}
dir_name=${2}
n=$((`find . -maxdepth 1 -type ${3} | wc -l`/$dir_size+1))
for i in `seq 1 $n`;
do
    mkdir -p "$dir_name$i";
    find . -maxdepth 1 -type ${3} | head -n $dir_size | xargs -i mv "{}" "$dir_name$i"
done
