#!/bin/bash

if [ $# -lt 2 ] || [ $# -gt 4 ];then
    echo "$0 <base RegulonDB directory> <base RegulonDBObjects directory> [id2gene map] [Python script]"
    exit 2
elif [ $# -eq 3 ];then
    map=$3
elif [ $# -eq 4 ];then
    map=$3
    script=$4
fi

base=$1
output=$2
: ${map:="feature2gene"} # conditional default value
: ${script="~/CodeBase/Development/pyorganism/scripts/store_id2gene.py"}

for ver in $(ls -d -- ${base}/[0-9].[0-9]/)
do
    number=`basename ${ver}`
    python ${script} ${base}/${number}/${map}.pkl ${output}/unified_maps.h5
done

