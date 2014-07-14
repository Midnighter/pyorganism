#!/bin/bash

if [ $# -gt 3 ];then
    echo "$0 [base RegulonDB directory] [base RegulonDBObjects directory] [Python script]"
    exit 2
fi

base=`dirname $0`
indir=$1
outdir=$2
: ${indir:="${PWD}/RegulonDB"}
: ${outdir:="${PWD}/RegulonDBObjects"}

${base}/parallel_compilation.sh ${indir} ${outdir}
${base}/parallel_name2gene.sh ${outdir}
${base}/parallel_feature2gene.sh ${outdir}
${base}/parallel_trn.sh ${outdir}
${base}/parallel_gpn.sh ${outdir}
${base}/store_all_maps.sh ${outdir} name2gene
${base}/store_all_maps.sh ${outdir} blattner2gene
${base}/store_all_maps.sh ${outdir} feature2gene
${base}/update_all_maps.sh ${outdir} name2gene 5.2
${base}/update_all_maps.sh ${outdir} blattner2gene 5.2
${base}/update_all_maps.sh ${outdir} feature2gene 8.5

