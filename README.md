# SBG VCF COMPARISON TOOL

SBG vcf comparison tool is written using the algorithm of [vcfeval](https://github.com/RealTimeGenomics/rtg-tools "vcfeval github page") using [this](http://biorxiv.org/content/early/2015/08/02/023754) paper. Please see the command line parameters below.

### Quick Install:
```
git clone https://gitlab.sbgdinc.com/Variants/Vcf-Comparison.git
make
```
**Note:** If you are unable to compile program or getting runtime error, it is mostly because the project is still under heavy development. Please shoot me an email(berke.toptas@sbgdinc.com) for any problems, suggestions etc.

#### Parameter format:
```
./sbgVcfComp -called <called_vcf_file> -baseline <baseline_vcf> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

#### Sample executions:

```
./sbgVcfComp -called called.vcf -base base.vcf -ref reference.fa -outDir SampleResultDir -filter none
```
   
```
./sbgVcfComp -called called2.vcf -base base2.vcf -ref reference.fa -outDir SampleResultDir -filter PASS -snp_only -sampleBase sample0 -sampleCalled sample01
```


## Command line parameters:


### -help
A **single** parameter mode for sbg tool which prints all parameter options to the console. (./sbgtool -help)


### -base <Baseline_vcf_path>

A **required** parameter for sbg tool which specifies the baseline vcf file path. It supports both bcf and vcf file formats.


### -called <Query_vcf_path>

A **required** parameter for sbg tool which specifies the query vcf file path. It supports both bcf and vcf file formats.


### -ref <Reference_fasta_path>

A **required** parameter for sbg tool which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.


### -outDir <Output_Directory_path>

A **required** parameter for sbg tool which specifies the output directory for program/error logs and ga4gh output vcf file. In current version, in order not to damage multiplatform capability, **directory should be created by user**.


### -sampleBase <Baseline_vcf_sample_name>

An **optional** parameter for sbg tool which is used to select sample name from baseline vcf file. By default, first sample is selected.


### -sampleCalled <Called_vcf_sample_name>

An **optional** parameter for sbg tool which is used to select sample name from called vcf file. By default, first sample is selected.

### -filter <Filter_name>

An **optional** parameter for sbg tool which is used to filter variants with given filter name. Filter name should be same in baseline and called variants. By default, **PASS filtering** is applied to the variants. In order to disable filtering, **'-filter none'** should be used.


### -snp_only

An **optional** parameter for sbg tool which is used to eliminate all of INDELs and SVs from input vcf files.


### -indel_only

An **optional** parameter for sbg tool which is used to eliminate all of SNPs and SVs from input vcf files.


### -ref_overlap [Almost Finished]

An **optional** parameter for sbg tool which tolerates reference overlaps to include more variants to the truth set.

### -platform-mode

An **optional** parameter for sbg tool which maximizes the thread count (to 25) assuming that there is no memory and processor limit. The default execution uses up to 4 threads and requires +4 GB memory to perform whole genome vcf file.
