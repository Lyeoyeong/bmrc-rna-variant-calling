#!/bin/bash

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

# Load software modules
module purge
module load Cromshell/0.4.2

# Parameters
SAMPLES=$1			# e.g. samples.ids
INPUTS_FOLDER=$2	# Path to where the JSon input files for each sample will be stored

# Some paths
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
WDL_PIPELINE="${SCRIPT_DIR}/pipeline/main.wdl"
JSON_TEMPLATE="${SCRIPT_DIR}/conf/input.template.json"
CROMWELL_OPTIONS="${SCRIPT_DIR}/conf/cromwell.options.json"

DATE=`date` && echo "[$DATE] Pipeline paths"
DATE=`date` && echo "[$DATE] .. WDL_PIPELINE     : ${WDL_PIPELINE}"
DATE=`date` && echo "[$DATE] .. JSON_TEMPLATE    : ${JSON_TEMPLATE}"
DATE=`date` && echo "[$DATE] .. CROMWELL_OPTIONS : ${CROMWELL_OPTIONS}"
echo

DATE=`date` && echo "[$DATE] Script arguments"
DATE=`date` && echo "[$DATE] .. SAMPLES       : ${SAMPLES}"
DATE=`date` && echo "[$DATE] .. INPUTS_FOLDER : ${INPUTS_FOLDER}"
echo

DATE=`date` && echo "[$DATE] Submitting Jobs to Cromwell Server"

mkdir -p $PWD/INPUTs
for SAMPLE in `cat ${SAMPLES}`
do
	INPUT_FILE="${INPUTS_FOLDER}/${SAMPLE}.json"
	sed "s/REPLACEME/${SAMPLE}/" ${JSON_TEMPLATE} > ${INPUT_FILE}
	cromshell submit ${WDL_PIPELINE} ${INPUT_FILE} ${CROMWELL_OPTIONS}
	sleep 5
done

echo
DATE=`date` && echo "[$DATE] Done"
exit 0
