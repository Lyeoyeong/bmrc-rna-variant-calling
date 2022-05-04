#!/bin/bash

SAMPLES=$1 # e.g. samples.ids

for SAMPLE in `cat ${SAMPLES}`

do
	DATE=`date` && echo "[$DATE] .. ${SAMPLE}"
	java -Dconfig.file=./conf/cromwell.conf -jar /users/sansom/rrl815/devel/cromwell/cromwell-51.jar run -i ./INPUTs/${SAMPLE}.json ./pipeline/main.single_end.wdl > rna-germline-variant-calling.${SAMPLE}.log
        DATE=`date` && echo "[$DATE] .. ${SAMPLE} Finished"
done

DATE=`date` && echo "[$DATE] Done"
