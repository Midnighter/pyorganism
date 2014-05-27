#!/bin/bash

if [ $# -lt 2 ];then
    echo "$0 <base RegulonDB directory> <output directory> [Python script]"
    exit 2
fi

base=$1
output=$2
script=$3
location=`dirname $0`
: ${script="${location}/store_network_statistics.py"}

for ver in $(ls -d -- ${base}/[0-9].[0-9]/)
do
    number=`basename ${ver}`
    echo ${number}
    python ${script} ${base}/${number} ${output}
done

