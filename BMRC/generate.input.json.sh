#!/bin/bash

set -e -o pipefail

SAMPLES=$1 # e.g. samples.ids
INPUTS_FOLDER=$2 # Path to where the JSon input files for each sample will be stored

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
JSON_TEMPLATE="${SCRIPT_DIR}/conf/input.template.json"

for SAMPLE in `cat ${SAMPLES}`
do
	INPUT_FILE="${INPUTS_FOLDER}/${SAMPLE}.json"
	sed "s/REPLACEME/${SAMPLE}/" ${JSON_TEMPLATE} > ${INPUT_FILE}
done


