#!/bin/bash


if [ $# -lt 1 ] || [ $# -gt 4 ];then
    echo "$0 <base RegulonDBObjects directory> [id2gene map] [base version] [Python script] "
    exit 2
elif [ $# -eq 2 ];then
    map=$2
elif [ $# -eq 3 ];then
    map=$2
    update=$3
elif [ $# -eq 4 ];then
    map=$2
    update=$3
    script=$4
fi

output=$1
: ${map:="feature2gene"} # conditional default value
: ${update:="8.5"} # conditional default value
: ${script="~/CodeBase/Development/pyorganism/scripts/update_maps.py"}

for ver in $(ls -d -- ${output}/[0-9].[0-9]/)
do
    number=`basename ${ver}`
    python ${script} ${output}/${number}/${map}.pkl ${output}/unified_maps.h5 ${update}
done

