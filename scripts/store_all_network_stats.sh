#!/bin/bash

if [ $# -lt 2 ] || [ $# -gt 3 ];then
    echo "$0 <RegulonDBObjects directory> <output directory> [Python script]"
    exit 2
fi

base=$1
output=$2
script=$3
: ${script:="`dirname $0`/store_network_statistics.py"}

for ver in $(ls -d -- ${base}/[0-9].[0-9]/)
do
    number=`basename ${ver}`
    python ${script} ${base}/${number} ${output}
done

