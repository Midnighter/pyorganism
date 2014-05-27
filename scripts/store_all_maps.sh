#!/bin/bash

if [ $# -lt 1 ] || [ $# -gt 3 ];then
    echo "$0 <base RegulonDBObjects directory> [id2gene map] [Python script]"
    exit 2
elif [ $# -eq 2 ];then
    map=$2
elif [ $# -eq 3 ];then
    map=$2
    script=$3
fi

output=$1
: ${map:="feature2gene"} # conditional default value
: ${script:="`dirname $0`/store_id2gene.py"}

for ver in $(ls -d -- ${output}/[0-9].[0-9]/)
do
    number=`basename ${ver}`
    python ${script} ${output}/${number}/${map}.pkl ${output}/unified_maps.h5
done

