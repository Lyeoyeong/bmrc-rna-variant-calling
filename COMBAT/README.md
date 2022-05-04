
# RNA Germline Variant Calling and Joint Genotyping

Workflow for processing RNA data for germline short variant discovery with GATK (v4) and related tools. The workflow will produce a gVCF (genomic VCF) file per individual. Once gVCF files have been created for multiple individuals, joint genotyping could be performed as well on a group of selected individuals.


## Installation 

**1.** In order to be able to run the pipeline, the following software should already be installed on your cluster:

- Java v11.0.2
- [cromwell v51](https://github.com/broadinstitute/cromwell/releases/tag/51)
- [cromshell v0.4.2](https://github.com/broadinstitute/cromshell/releases/tag/0.4.2)
- [GATK v4.1.7.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.7.0)
- [STAR v2.7.3a](https://github.com/alexdobin/STAR/releases/tag/2.7.3a)

**2.** Clone the repository:

```bash
git clone https://github.com/COMBATOxford/rna-germline-variant-calling
cd rna-germline-variant-calling
```

**3.** Modify the `cromwell-server` configuration file to work in your cluster by replacing the placeholders `YOUR_QUEUE_NAME` and `YOUR_PROJECT_NAME` with the appropriate values:

```bash
sed 's/YOUR_QUEUE_NAME/queue_name/;s/YOUR_PROJECT_NAME/project_name/' conf/cromwell.conf
```

**4.** Set up `cromshell` for the first time by running `cromshell status`. You will be prompted to enter the URL of the `cromwell` server (e.g. `http://myserver:8000`). This will be written in the file `~/.cromshell/cromwell_server.config`. Alternatively, you could export the `CROMWELL_URL` environment variable like `export CROMWELL_URL=http://myserver:8000`.

**5.** Update the following in the `conf/input.template.json` file:

- `RNAseq.sequencingDetails`: replace the prefix path `/path/to/inputs` for the proper path to where your `TSVs` folder will be located.
- `RNAseq.outputFolder`: replace the prefix path `/path/to/outputs` for the proper path to where you want the `variant-calling` output folder to be located.
- `RNAseq.refFasta`
- `RNAseq.refFastaIndex`
- `RNAseq.refDict`
- `RNAseq.dbSnpVcf`
- `RNAseq.dbSnpVcfIndex`
- `RNAseq.knownVcfs`
- `RNAseq.knownVcfsIndices`
- `RNAseq.annotationsGTF`
- `RNAseq.wgsCallingIntervalList`
- `RNAseq.starReferencePath`: check each of the above mentioned paths to see that they are pointing to the appropriate reference file or folder.
- `RNAseq.readLength`: replace the placeholder `read_length` for the appropriate read length of the sequencing reads.

Feel free to additionally check the other parameters as well.

**6.** Update the `conf/cromwell.options.json` file by replacing the prefix path `/path/to` for the proper path to where you would like the pipeline `LOGs` folder to be located.

**7.** If you want to build and use the same reference used for the over-arching publication, check the `Preparing the references` section in the [COMBATOxford/tenx-preprocessing](https://github.com/COMBATOxford/tenx-preprocessing) repository.

**8.** Done! Now you are ready to run your pipelines.

_NOTE1: the procedure had only been tested with the above mentioned software versions._

_NOTE2: the scripts are invoking the tools using the `module` directive because this is how they have been installed in our cluster._


## Preparing configuration files

### Samples submission file

In order to make a bulk submission (it will work for one sample as well), create a plain text file (I like to call it `samples.ids`) and write all your `sample IDs` one per line, e.g.

```
COMBAT0019
COMBAT0045
COMBAT0125
COMBAT0126
```

### Sample configuration TSV file

In order for the pipeline to find where the FastQ files are for each sample as well as to know which values to provide to the BAM readgroup tags, a `tsv` (tab-separated value) file has to be created for each sample. Each line on this file will correspond to a pair of FastQ files. For example, if a sample has been sequenced 3 times in order to get enough coverage, than its `tsv` file should contain 3 lines. The requirements for this file are:

- Sample Name: the `Sample Name` or `Sample ID` should be used to name the file, e.g. `COMBAT0019.tsv`
- Specification file in TSV format with the following fields having one sequenced sample per row:
    - Read1
    - Read2
    - Project Number
    - Sample Name
    - Library Name
    - Read Group Name
    - Sequencing Center
    - Platform
    - Platform Model
    - Platform Unit
    - Run Date

```
e.g. (skip header row)
Read1	Read2	ProjectNumber	SampleName	LibraryName	ReadGroupName	SequencingCenter	Platform	PlatformModel	PlatformUnit	RunDate
/path/to/FASTQ/HV3M3BBXX/SampleXXX_L004_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_123456_70655002	WHG-OGC	Illumina	HiSeq4000	HV3M3BBXX.4.70655002	2018-04-24
/path/to/FASTQ/HTJJ2BBXX/SampleXXX_L005_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_987654_70655002	WHG-OGC	Illumina	HiSeq4000	HTJJ2BBXX.5.70655002	2018-04-26
/path/to/FASTQ/HTJJ2BBXX/SampleXXX_L006_R1.fastq.gz		PROJECT	SampleXXX	LIB5985A9	WTCHG_987655_70655002	WHG-OGC	Illumina	HiSeq4000	HTJJ2BBXX.6.70655002	2018-04-26
```

### Joint Genotyping submission file

In order to make a bulk submission (it will work for one sample as well), create a plain text file (I like to call it `samples.joint_genotyping.list`) and write one submission per line, using the format below:

- `Group prefix`: prefix to use for the resulting joint genotyping VCF file.
- `Reference genome`: path to the genome reference FastA file.
- `dbSNP`: path to the reference dbSNP VCF file.
- `Output path`: path where to store the merged VCF files.
- `VCF files list`: path to a file listing all available gVCF files (one per line, using absolute paths).
- `Samples list`: comma separated list of the `Sample IDs` to perform the joint genotyping on.

The `Sample IDs` will be used to get the gVCF files to use by grepping them from the gVCF files list, so be sure that the gVCF files are named using the `Sample ID`.

```
e.g. here are three groups to perform joint genotyping on:
gPlexA	/path/to/references/GRCh38/Sequence/genome.fa	/path/to/references/GRCh38/Annotation/Variation/Homo_sapiens_assembly38.dbsnp138.vcf.gz	/path/to/myproject/joint_genotyping	/path/to/myproject/samples.genotypes.gvcf.list	COMBAT0007,COMBAT0027,COMBAT0030,COMBAT0049,COMBAT0058,COMBAT0063,COMBAT0071,COMBAT0080,COMBAT0088,COMBAT0093,COMBAT0114,COMBAT0115,COMBAT0125,COMBAT0127
gPlexB	/path/to/references/GRCh38/Sequence/genome.fa	/path/to/references/GRCh38/Annotation/Variation/Homo_sapiens_assembly38.dbsnp138.vcf.gz	/path/to/myproject/joint_genotyping	/path/to/myproject/samples.genotypes.gvcf.list	COMBAT0004,COMBAT0025,COMBAT0028,COMBAT0054,COMBAT0059,COMBAT0070,COMBAT0073,COMBAT0074,COMBAT0077,COMBAT0104,COMBAT0112,COMBAT0118,COMBAT0126,COMBAT0128
gPlexC	/path/to/references/GRCh38/Sequence/genome.fa	/path/to/references/GRCh38/Annotation/Variation/Homo_sapiens_assembly38.dbsnp138.vcf.gz	/path/to/myproject/joint_genotyping	/path/to/myproject/samples.genotypes.gvcf.list	COMBAT0010,COMBAT0031,COMBAT0033,COMBAT0048,COMBAT0053,COMBAT0067,COMBAT0068,COMBAT0076,COMBAT0079,COMBAT0085,COMBAT0088,COMBAT0090,COMBAT0121,COMBAT0129
```


## Usage/Examples

**1.** Open two new independent shell terminals.

**2.** In the first terminal, start `Cromwell` in server mode:

```bash
export CROMJAR="${EBROOTCROMWELL}/cromwell-${EBVERSIONCROMWELL}.jar"
export CROMCONF="/path/to/rna-germline-variant-calling/conf/cromwell.conf"

module purge
module load Cromwell/51-Java-11.0.2
java -Dconfig.file=${CROMCONF} -jar ${CROMJAR} server
```

**3.** Use the second terminal to submit your Variant Calling jobs for each of your samples:

```bash
/path/to/rna-germline-variant-calling/submit_jobs.variant_calling.sh /path/to/myproject/samples.ids /path/to/myproject/INPUTs
```

**4.** Once the Variant Calling jobs are done, you could proceed to create a Joint Genotyping VCF file for a group of individuals by running:

```bash
qsub \
    -t 1-$(cat /path/to/myproject/samples.joint_genotyping.list | wc -l) \
    -P YOUR_PROJECT_NAME \
    -q YOUR_QUEUE_NAME \
    /path/to/rna-germline-variant-calling/submit_jobs.joint_genotyping.sh /path/to/myproject/samples.joint_genotyping.list
```


## Authors

- [@santiagorevale](https://www.github.com/santiagorevale)

  
## License

Copyright © 2020 [Santiago Revale](https://www.github.com/santiagorevale).
This project is [MIT](LICENSE) licensed.

  
## References

The RNAseq Germline Variant Calling workflow is based on v1.0.0 of `gatk4-rnaseq-germline-snps-indels` workflow:

- Gandham B. gatk-workflows/gatk4-rnaseq-germline-snps-indels. 2019. URL: https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/tree/1.0.0

Here are the references to the software being used by the pipeline:

- Voss K, Van der Auwera G and Gentry J. Full-stack genomics pipelining with GATK4 + WDL + Cromwell [version 1; not peer reviewed]. F1000Research 2017, 6(ISCB Comm J):1381 (slides). doi: [10.7490/f1000research.1114634.1](https://doi.org/10.7490/f1000research.1114634.1)
- Smith J, Bergelson L, et al. broadinstitute/cromshell. 2020. URL: https://github.com/broadinstitute/cromshell
- McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res, 20:1297-303. doi: [10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110).
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635). Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.
