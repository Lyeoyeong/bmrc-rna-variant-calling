#!/bin/bash

#$ -cwd
#$ -N joint_genotyping -j y
#$ -pe shmem 2

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can
# read/write what it's created by the script
umask 002

# Load software modules
module purge
module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

# Job Arguments
SAMPLES_FILE=$1

# Task Arguments
GROUP=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $SAMPLES_FILE)       # Prefix for merged VCF file
REF_FILE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$2}" $SAMPLES_FILE)    # Path to the reference FastA file
DBSNP_FILE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$3}" $SAMPLES_FILE)  # Path to the reference dbSNP VCF file
OUT_PATH=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$4}" $SAMPLES_FILE)    # Path where to store the merged VCF files
VCF_PATH=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$5}" $SAMPLES_FILE)    # Path to a file listing all VCF files
SAMPLES=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$6}" $SAMPLES_FILE)     # Comma separated Sample IDs

echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "SGE job ID       : "$JOB_ID
echo "SGE task ID      : "$SGE_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "Samples File     : ${SAMPLES_FILE}"
echo
echo "********************************************************"
echo "* Task Parameters"
echo "********************************************************"
echo "Group ID         : ${GROUP}"
echo "Reference file   : ${REF_FILE}"
echo "Output path      : ${OUT_PATH}"
echo "VCF files list   : ${VCF_PATH}"
echo "Samples          : ${SAMPLES}"
echo

mkdir -p ${OUT_PATH} && cd ${OUT_PATH}

DATE=`date` && echo "[$DATE] Combining gVCF files"

VCF_LIST="${OUT_PATH}/${GROUP}.vcf.list"
DATE=`date` && echo "[$DATE] .. creating gVCF file list: ${VCF_LIST}"
rm -f ${VCF_LIST}

# In order to avoid having a super-long command by providing many VCF files with
# their path, I create a temporary directory where I link the VCF files and run 
# the analysis. In the end, I copy the output from here to the final folder.
TEMP_DIR=`mktemp -d -p ${OUT_PATH}`
if [ ! -d "${TEMP_DIR}" ]; then
        echo "Unable to create temporary directory"
        exit 1
fi
cd ${TEMP_DIR}

VCF_OPT=""
for SAMPLE in `echo "${SAMPLES}" | tr ',' "\n"`
do
    VCF_FILE=`grep -w "${SAMPLE}" ${VCF_PATH}`
    ln -s ${VCF_FILE}
    ln -s ${VCF_FILE}.tbi
    echo "${VCF_FILE}" >> ${VCF_LIST}
    VCF_FILE=`basename ${VCF_FILE}`
    VCF_OPT="${VCF_OPT} -V ${VCF_FILE}"
done

DATE=`date` && echo "[$DATE] .. combining files"
gatk CombineGVCFs -O combined.g.vcf.gz -R ${REF_FILE} ${VCF_OPT}

DATE=`date` && echo "[$DATE] .. genotyping files"
gatk GenotypeGVCFs -O genotyped.vcf.gz -R ${REF_FILE} -D ${DBSNP_FILE} -V combined.g.vcf.gz

DATE=`date` && echo "[$DATE] .. filtering files"
# ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
# than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
EXCESS_HET=54.69
gatk VariantFiltration -O ${OUT_PATH}/${GROUP}.vcf.gz -V genotyped.vcf.gz --filter-name ExcessHet --filter-expression "ExcessHet > ${EXCESS_HET}"
      

# Delete temporary directory
#cd ${OUT_PATH}
#rm -Rf ${TEMP_DIR}

echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
