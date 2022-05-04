#!/bin/bash

cd /well/singlecell/P170676/RNAseq-genotyping/variant_calling/rna-germline-variant-calling
ls $PWD/*/*.variant_filtered.vcf.gz > vcf_files.list
/apps/well/bcftools/1.4.1/bin/bcftools merge -l vcf_files.list -O z -o variant_filtered.merged.vcf.gz
/apps/well/htslib/1.6/bin/tabix variant_filtered.merged.vcf.gz
