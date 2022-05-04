## Copyright (c) 2020, Santiago Revale
## 
## The current file has been built on top of the following one:
## https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/1.0.0/gatk4-rna-best-practices.wdl
## 
## Below is the license from the original file.
## 
## ==============================================================================
## Copyright (c) 2019, Bhanu Gandham
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## 1. Redistributions of source code must retain the above copyright notice, this
##    list of conditions and the following disclaimer.
## 
## 2. Redistributions in binary form must reproduce the above copyright notice,
##    this list of conditions and the following disclaimer in the documentation
##    and/or other materials provided with the distribution.
## 
## 3. Neither the name of the copyright holder nor the names of its
##    contributors may be used to endorse or promote products derived from
##    this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## ==============================================================================
##
##
## Workflow for processing RNA data for germline short variant discovery with GATK (v4) and related tools
##
## Requirements/expectations :
## - Sample Name
## - Specification file in TSV format with the following fields having one sequenced sample per row:
##     Read1
##     Read2
##     ProjectNumber
##     SampleName
##     LibraryName
##     ReadGroupName
##     SequencingCenter
##     Platform
##     PlatformModel
##     PlatformUnit
##     RunDate
##
##   e.g. (skip header row)
##   Read1	Read2	ProjectNumber	SampleName	LibraryName	ReadGroupName	SequencingCenter	Platform	PlatformModel	PlatformUnit	RunDate
##   /path/to/FASTQ/HV3M3BBXX/SampleXXX_L004_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_123456_70655002	WHG-OGC	Illumina	HiSeq4000	HV3M3BBXX.4.70655002	2018-04-24
##   /path/to/FASTQ/HTJJ2BBXX/SampleXXX_L005_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_987654_70655002	WHG-OGC	Illumina	HiSeq4000	HTJJ2BBXX.5.70655002	2018-04-26
##   /path/to/FASTQ/HTJJ2BBXX/SampleXXX_L006_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_987655_70655002	WHG-OGC	Illumina	HiSeq4000	HTJJ2BBXX.6.70655002	2018-04-26
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##

workflow RNAseq {

	String sampleName
	File sequencingDetails
	Array[Array[String]] sequencingDetailsTable = read_tsv(sequencingDetails)
	String outputFolder
	String finalOutputPath = outputFolder + '/' + sampleName

	File refFasta
	File refFastaIndex
	File refDict

	String gatk_path
	String star_path

	Array[File] knownVcfs
	Array[File] knownVcfsIndices

	File dbSnpVcf
	File dbSnpVcfIndex

	Int? minConfidenceForVariantCalling

	File? wgsCallingIntervalList
	File? exomeCallingIntervalList

	## Inputs for STAR
	Int? readLength
	File? starReferencePath
	File annotationsGTF
  
	## Optional user optimizations
	Int? haplotypeScatterCount
	Int scatterCount = select_first([haplotypeScatterCount, 6])

	if (!defined(exomeCallingIntervalList)) {
		if (!defined(wgsCallingIntervalList)) {
			call gtfToCallingIntervals {
				input:
					gtf = annotationsGTF,
					ref_dict = refDict,
					gatk_path = gatk_path
			}
		}
	}
	File intervalList = select_first([exomeCallingIntervalList, wgsCallingIntervalList, gtfToCallingIntervals.interval_list, ""])

	if (!defined(starReferencePath)) {
		call StarGenerateReferences { 
			input:
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				annotations_gtf = annotationsGTF,
				read_length = readLength,
				star_path = star_path
		}
	}
	File starReferences = select_first([starReferencePath, StarGenerateReferences.star_genome_ref, ""])

	# Read sequencing details file and loop through every run FastQ file/s
	scatter (detailsRow in sequencingDetailsTable) {
		File fastq_read1 = detailsRow[0]
		String project_number = detailsRow[1]
		String sample_name = detailsRow[2]
		String library_name = detailsRow[3]
		String read_group_name = detailsRow[4]
		String sequencing_center = detailsRow[5]
		String platform = detailsRow[6]
		String platform_model = detailsRow[7]
		String platform_unit = detailsRow[8]
		String run_date = detailsRow[9]

		call FastqToUbam {
			input:
				fastq_read1 = fastq_read1,
				sample_name = sample_name,
				library_name = library_name, 
				read_group_name = read_group_name,
				sequencing_center = sequencing_center,
				platform = platform,
				platform_model = platform_model,
				platform_unit = platform_unit,
				run_date = run_date,
				gatk_path = gatk_path
		}
	}

	call MergeUnalignedBams {
		input:
			unmapped_bams = FastqToUbam.output_bam,
			sample_name = sampleName,
			gatk_path = gatk_path
	}

	call StarAlign { 
		input: 
			star_genome_ref = starReferences,
			fastq1 = FastqToUbam.output_fastq_read1,
			base_name = sampleName + ".star",
			star_path = star_path
	}

	call MergeBamAlignment {
		input: 
			unaligned_bam = MergeUnalignedBams.output_bam,
			star_bam = StarAlign.output_bam,
			base_name = sampleName + ".merged",
			ref_fasta = refFasta,
			ref_dict = refDict,
			gatk_path = gatk_path
	}

	call MarkDuplicates {
		input:
			input_bam = MergeBamAlignment.output_bam,
			base_name = sampleName + ".dedupped",
			gatk_path = gatk_path
	}

    call SplitNCigarReads {
        input:
            input_bam = MarkDuplicates.output_bam,
            input_bam_index = MarkDuplicates.output_bam_index,
            base_name = sampleName + ".split",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            interval_list = intervalList,
            gatk_path = gatk_path
    }

	call BaseRecalibrator {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
			gatk_path = gatk_path
	}

	call ApplyBQSR {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			gatk_path = gatk_path
	}

	call ScatterIntervalList {
		input:
			interval_list = intervalList,
			scatter_count = scatterCount,
			gatk_path = gatk_path
	}

	scatter (interval in ScatterIntervalList.out) {
		call HaplotypeCaller {
			input:
				input_bam = ApplyBQSR.output_bam,
				input_bam_index = ApplyBQSR.output_bam_index,
				base_name = sampleName + ".hc",
				interval_list = interval,
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				ref_dict = refDict,
				dbSNP_vcf = dbSnpVcf,
				dbSNP_vcf_index = dbSnpVcfIndex,
				stand_call_conf = minConfidenceForVariantCalling,
				gatk_path = gatk_path
		}

		File HaplotypeCallerOutputVcf = HaplotypeCaller.output_vcf
		File HaplotypeCallerOutputVcfIndex = HaplotypeCaller.output_vcf_index
	}

	call MergeVCFs {
		input:
			input_vcfs = HaplotypeCallerOutputVcf,
			input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
			output_vcf_name = sampleName + ".g.vcf.gz",
			gatk_path = gatk_path
	}

	call VariantFiltration {
		input:
			input_vcf = MergeVCFs.output_vcf,
			input_vcf_index = MergeVCFs.output_vcf_index,
			base_name = sampleName + ".variant_filtered.vcf.gz",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			gatk_path = gatk_path
	}

	call copyFinalOutput {
		input:
			outputs = [ApplyBQSR.output_bam, ApplyBQSR.output_bam_index, MergeVCFs.output_vcf, MergeVCFs.output_vcf_index, VariantFiltration.output_vcf, VariantFiltration.output_vcf_index],
			outputs_dir = finalOutputPath
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_vcf = MergeVCFs.output_vcf
		File merged_vcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_vcf = VariantFiltration.output_vcf
		File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
	}
}


task gtfToCallingIntervals {
	File gtf
	File ref_dict

	String output_name = basename(gtf, ".gtf") + ".exons.interval_list"

	String gatk_path

	command <<<
		module purge
		module load "${gatk_path}"
		module load R/3.6.0-foss-2018b

		Rscript --no-save -<<'RCODE'
			gtf = read.table("${gtf}", sep="\t")
			gtf = subset(gtf, V3 == "exon")
			write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
		RCODE

		awk '{print $1 "\t" ($2 - 1) "\t" $3}' exome.bed > exome.fixed.bed

		gatk \
			BedToIntervalList \
			-I=exome.fixed.bed \
			-O=${output_name} \
			-SD=${ref_dict}
	>>>

	output {
		File interval_list = "${output_name}"
	}

	runtime {
		cpu: 3
	}
}

# Generate uBAM from FASTQ
task FastqToUbam {

	File fastq_read1

	String sample_name
	String library_name
	String read_group_name
	String sequencing_center
	String platform
	String platform_model
	String platform_unit
	String run_date

	String gatk_path
	String java_opt

	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}" \
			FastqToSam \
			--FASTQ ${fastq_read1} \
			--OUTPUT ${read_group_name}.unaligned.bam \
			--SAMPLE_NAME ${sample_name} \
			--LIBRARY_NAME ${library_name} \
			--READ_GROUP_NAME ${read_group_name} \
			--SEQUENCING_CENTER ${sequencing_center} \
			--PLATFORM ${platform} \
			--PLATFORM_MODEL ${platform_model} \
			--PLATFORM_UNIT ${platform_unit} \
			--RUN_DATE ${run_date}
	>>>

	output {
		File output_fastq_read1 = fastq_read1
		File output_bam = "${read_group_name}.unaligned.bam"
	}

	runtime {
		cpu: 4
	}
}

# Merge unaligned bams for the same sample into one ubam
task MergeUnalignedBams {

	Array[File] unmapped_bams
	String sample_name

	String gatk_path
	String java_opt

	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}" \
			MergeSamFiles \
			--SORT_ORDER queryname \
			--INPUT ${sep=' --INPUT ' unmapped_bams} \
			--OUTPUT ${sample_name}.unaligned.bam
	>>>

	output {
		File output_bam = "${sample_name}.unaligned.bam"
	}

	runtime {
		cpu: 4
	}
}

task StarGenerateReferences {

	File ref_fasta
	File ref_fasta_index
	File annotations_gtf
	Int? read_length  ## Should this be an input, or should this always be determined by reading the first line of a fastq input

	Int? num_threads
	Int threads = select_first([num_threads, 8])

	String star_path

	command <<<
		set -e

		module purge
		module load "${star_path}"

		mkdir STARindex
		STAR \
			--runMode genomeGenerate \
			--genomeDir STARindex \
			--genomeFastaFiles ${ref_fasta} \
			--sjdbGTFfile ${annotations_gtf} \
			${"--sjdbOverhang "+(read_length-1)} \
			--runThreadN ${threads}

		ls STAR2_5
	>>>

	output {
		Array[File] star_logs = glob("*.out")
		File star_genome_ref = "STARindex"
	}

	runtime {
		cpu: threads
	}
}

task StarAlign {

	File star_genome_ref
	Array[File] fastq1
	String base_name

	Int? num_threads
	Int threads = select_first([num_threads, 8])
	Int? star_mem_max_gb
	Int star_mem = select_first([star_mem_max_gb, 45])
	#Is there an appropriate default for this?
	Int? star_limitOutSJcollapsed

	String star_path

	command <<<
		set -e

		module purge
		module load "${star_path}"

		STAR \
			--genomeDir ${star_genome_ref} \
			--runThreadN ${threads} \
			--readFilesIn ${sep=',' fastq1} --readFilesCommand "gunzip -c" \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--limitBAMsortRAM ${star_mem+"000000000"} \
			--limitOutSJcollapsed ${default=1000000 star_limitOutSJcollapsed} \
			--outFileNamePrefix ${base_name}.
	>>>

	output {
		File output_bam = "${base_name}.Aligned.sortedByCoord.out.bam"
		File output_log_final = "${base_name}.Log.final.out"
		File output_log = "${base_name}.Log.out"
		File output_log_progress = "${base_name}.Log.progress.out"
		File output_SJ = "${base_name}.SJ.out.tab"
	}

	runtime {
		cpu: threads
	}
}

task MergeBamAlignment {

	File ref_fasta
	File ref_dict

	File unaligned_bam
	File star_bam
	String base_name

	String gatk_path

	command <<<
		module purge
		module load "${gatk_path}"

		gatk \
			MergeBamAlignment \
			--REFERENCE_SEQUENCE ${ref_fasta} \
			--UNMAPPED_BAM ${unaligned_bam} \
			--ALIGNED_BAM ${star_bam} \
			--OUTPUT ${base_name}.bam \
			--INCLUDE_SECONDARY_ALIGNMENTS false \
			--PAIRED_RUN False \
			--VALIDATION_STRINGENCY SILENT
	>>>

	output {
		File output_bam="${base_name}.bam"
	}

	runtime {
		cpu: 8
	}
}

task MarkDuplicates {

	File input_bam
	String base_name

	String gatk_path

	command <<<
		module purge
		module load "${gatk_path}"

		gatk \
			MarkDuplicates \
			--INPUT ${input_bam} \
			--OUTPUT ${base_name}.bam  \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY SILENT \
			--METRICS_FILE ${base_name}.metrics
	>>>

	output {
		File output_bam = "${base_name}.bam"
		File output_bam_index = "${base_name}.bai"
		File metrics_file = "${base_name}.metrics"
	}

	runtime {
		cpu: 8
	}
}

task SplitNCigarReads {

	File input_bam
	File input_bam_index
	String base_name
	File interval_list

	File ref_fasta
	File ref_fasta_index
	File ref_dict

	String gatk_path

	command <<<
		module purge
		module load "${gatk_path}"

		gatk \
			SplitNCigarReads \
			-R ${ref_fasta} \
			-I ${input_bam} \
			-O ${base_name}.bam 
	>>>

	output {
		File output_bam = "${base_name}.bam"
		File output_bam_index = "${base_name}.bai"
	}

	runtime {
		cpu: 8
	}
}

task BaseRecalibrator {

	File input_bam
	File input_bam_index
	String recal_output_file

	File dbSNP_vcf
	File dbSNP_vcf_index
	Array[File] known_indels_sites_VCFs
	Array[File] known_indels_sites_indices

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	String gatk_path
	String java_opt

	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}" \
			BaseRecalibrator \
			-R ${ref_fasta} \
			-I ${input_bam} \
			--use-original-qualities \
			-O ${recal_output_file} \
			-known-sites ${dbSNP_vcf} \
			-known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
	>>>

	output {
		File recalibration_report = recal_output_file
	}

	runtime {
		cpu: 4
	}
}

task ApplyBQSR {

	File input_bam
	File input_bam_index
	String base_name
	File recalibration_report

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	String gatk_path
	String java_opt

	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}" \
			ApplyBQSR \
			--add-output-sam-program-record \
			-R ${ref_fasta} \
			-I ${input_bam} \
			--use-original-qualities \
			-O ${base_name}.bam \
			--bqsr-recal-file ${recalibration_report}
	>>>

	output {
		File output_bam = "${base_name}.bam"
		File output_bam_index = "${base_name}.bai"
	}

	runtime {
		cpu: 4
	}
}

task ScatterIntervalList {

	File interval_list
	Int scatter_count

	String gatk_path
	String java_opt

	command <<<
		set -e

		module purge
		module load "${gatk_path}"

		mkdir out
		gatk --java-options "${java_opt}" \
			IntervalListTools \
			--SCATTER_COUNT=${scatter_count} \
			--SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
			--UNIQUE=true \
			--SORT=true \
			--INPUT=${interval_list} \
			--OUTPUT=out

		module load python/3.5.2-gcc5.4.0

		python3 <<CODE
		import glob, os
		# Works around a JES limitation where multiples files with the same name overwrite each other when globbed
		intervals = sorted(glob.glob("out/*/*.interval_list"))
		for i, interval in enumerate(intervals):
			(directory, filename) = os.path.split(interval)
			newName = os.path.join(directory, str(i + 1) + filename)
			os.rename(interval, newName)
		print(len(intervals))
		f = open("interval_count.txt", "w+")
		f.write(str(len(intervals)))
		f.close()
		CODE
	>>>

	output {
		Array[File] out = glob("out/*/*.interval_list")
		Int interval_count = read_int("interval_count.txt")
	}

	runtime {
		cpu: 2
	}
}

task HaplotypeCaller {

	File input_bam
	File input_bam_index
	String base_name

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	Int? stand_call_conf

	String gatk_path
	String java_opt

	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}" \
			HaplotypeCaller \
			-R ${ref_fasta} \
			-I ${input_bam} \
			-L ${interval_list} \
			-O ${base_name}.vcf.gz \
			-ERC GVCF \
			-dont-use-soft-clipped-bases \
			--standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
			--dbsnp ${dbSNP_vcf}
	>>>

	output {
		File output_vcf = "${base_name}.vcf.gz"
		File output_vcf_index = "${base_name}.vcf.gz.tbi"
	}

	runtime {
		cpu: 4
	}
}

task MergeVCFs {

	Array[File] input_vcfs
	Array[File] input_vcfs_indexes
	String output_vcf_name

	String gatk_path
	String java_opt

	# Using MergeVcfs instead of GatherVcfs so we can create indices
	# See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
	command <<<
		module purge
		module load "${gatk_path}"

		gatk --java-options "${java_opt}"  \
			MergeVcfs \
			--INPUT ${sep=' --INPUT ' input_vcfs} \
			--OUTPUT ${output_vcf_name}
	>>>

	output {
		File output_vcf = output_vcf_name
		File output_vcf_index = "${output_vcf_name}.tbi"
	}

	runtime {
		cpu: 2
	}
}

task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_name

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	String gatk_path

	command <<<
		module purge
		module load "${gatk_path}"

		gatk \
			VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "HardFiltered" \
			--filter-expression "QD < 2.0 || FS > 30.0" \
			-O ${base_name}
	>>>

	output {
		File output_vcf = "${base_name}"
		File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
		cpu: 4
	}
}

task copyFinalOutput {

	Array[String] outputs
	String outputs_dir

	command {
		mkdir -p ${outputs_dir}
		cp ${sep=' ' outputs} ${outputs_dir}
	}
}
