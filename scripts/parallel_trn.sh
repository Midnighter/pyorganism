#!/bin/bash


if [ $# -lt 1 ] || [ $# -gt 2 ];then
    echo "$0 <base RegulonDBObjects directory> [Python script]"
    exit 2
fi

output=$1
script=$2
: ${script="~/CodeBase/Development/pyorganism/scripts/construct_trn.py"}

ls -d -- ${output}/*/ | parallel --group -j+0 python ${script} ${output}/{/}

