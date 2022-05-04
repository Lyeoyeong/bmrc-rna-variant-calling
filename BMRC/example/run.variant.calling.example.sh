#!/bin/bash

SAMPLES=$1 # e.g. samples.ids

for SAMPLE in `cat ${SAMPLES}`

do
	DATE=`date` && echo "[$DATE] .. ${SAMPLE}"
	java -Dconfig.file=/users/sansom/rrl815/work/01_AS_P170676/demultiplexing/rna_variant_calling/rna-germline-variant-calling/conf/cromwell.conf -jar /users/sansom/rrl815/devel/cromwell/cromwell-51.jar run -i /gpfs3/well/sansom/users/rrl815/work/01_AS_P170676/demultiplexing/rna_variant_calling/INPUTs/${SAMPLE}.json /gpfs3/well/sansom/users/rrl815/work/01_AS_P170676/demultiplexing/rna_variant_calling/rna-germline-variant-calling/pipeline/main.se.wdl > rna-germline-variant-calling.${SAMPLE}.log
done

DATE=`date` && echo "[$DATE] Done"
