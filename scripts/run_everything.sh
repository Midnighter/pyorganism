#!/bin/bash

base="~/CodeBase/Development/pyorganism/scripts"
indir="RegulonDB"
outdir="RegulonDBObjects"

${base}/parallel_basics.sh ${indir} ${outdir}
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

