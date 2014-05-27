#!/bin/bash


if [ $# -lt 2 ] || [ $# -gt 3 ];then
    echo "$0 <base RegulonDB directory> <base RegulonDBObjects directory> [Python script]"
    exit 2
fi

base=$1
output=$2
script=$3
: ${script="~/CodeBase/Development/pyorganism/scripts/compile_regulondb_content.py"}

ls -d -- ${base}/*/ | parallel --group -j+0 python ${script} ${base}/{/} ${output}/{/}

