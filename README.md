# VARIANT BENCHMARKING TOOLS (VBT)

VBT is set of features that is used for aligner + variant calling benchmarking. Currently there are two tools under vbt :

    1. VARIANT COMPARISON
    2. MENDELIAN VIOLATION DETECTOR


### Quick Install:
```
git clone https://gitlab.sbgdinc.com/Variants/vbt.git
cd vbt
make
```
**Note:** If you are unable to compile program or getting runtime error, it is mostly because the project is still under heavy development. Please shoot me an email(berke.toptas@sbgdinc.com) for any problems, suggestions etc.

#### Parameter format:
```
./vbt varComp [PARAMETERS]
./vbt mendelian [PARAMETERS]
```
#### Warning: #### Parameters are **case sensitive**

## 1. VARIANT COMPARISON

Variant Comparison tool is written using the algorithm of [vcfeval](https://github.com/RealTimeGenomics/rtg-tools "vcfeval github page") using [this](http://biorxiv.org/content/early/2015/08/02/023754) paper. Please see the command line parameters below.

#### Parameter Format:
```
./vbt varcomp -called <called_vcf_file> -base <baseline_vcf> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

## Command line parameters:


### --help
A **single** parameter mode which prints all parameter options to the console. (./vbt varComp --help)


### -base <Baseline_vcf_path>

A **required** parameter which specifies the baseline vcf file path. It supports both bcf and vcf file formats.


### -called <Query_vcf_path>

A **required** parameter which specifies the query vcf file path. It supports both bcf and vcf file formats.


### -ref <Reference_fasta_path>

A **required** parameter which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.


### -outDir <Output_Directory_path>

A **required** parameter which specifies the output directory for program/error logs and ga4gh output vcf file. In current version, in order not to damage multiplatform capability, **directory should be created by user**.

### -output-mode <Output_Mode>

An **optional** parameter which specifies the way output will generated. **SPLIT** mode creates 4 vcf files(TPbase, TPcalled, FP, FN). **GA4GH** mode creates a single merged vcf file. Default valus is **SPLIT**.

### --allele-match

An **optinal** parameter which changes comparison mode to Allele Matching

### -sampleBase <Baseline_vcf_sample_name>

An **optional** parameter which is used to select sample name from baseline vcf file. By default, first sample is selected.


### -sampleCalled <Called_vcf_sample_name>

An **optional** parameter which is used to select sample name from called vcf file. By default, first sample is selected.

### -filter <Filter_name>

An **optional** parameter which is used to filter variants with given filter name. Filter name should be same in baseline and called variants. By default, **PASS filtering** is applied to the variants. In order to disable filtering, **'-filter none'** should be used.


### --snp_only

An **optional** parameter which is used to eliminate all of INDELs and SVs from input vcf files.


### --indel_only

An **optional** parameter which is used to eliminate all of SNPs and SVs from input vcf files.


### --ref_overlap [Almost Finished]

An **optional** parameter which tolerates reference overlaps to include more variants to the truth set.

### --platform-mode

An **optional** parameter which maximizes the thread count (to 25) assuming that there is no memory and processor limit. The default execution uses up to 4 threads and requires +4 GB memory to perform whole genome vcf file.

### -thread-count <[1-25]>
An **optional** parameter to specify number of threads. Default value is 2


## 2. MENDELIAN VIOLATION DETECTOR

A tool for checking all mendelian violations in parent-child data.

#### Parameter Format:
```
./vbt mendelian -father <father_vcf_path> -mother <mother_vcf_path> -child <child_vcf_path> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

## Command line parameters:


### --help
A **single** parameter mode which prints all parameter options to the console. (./vbt mendelian --help)


### -father <father_vcf_path>

A **required** parameter which specifies the father vcf file path. It supports both bcf and vcf file formats.


### -mother <mother_vcf_path>

A **required** parameter which specifies the mother vcf file path. It supports both bcf and vcf file formats.

### -child <child_vcf_path>

A **required** parameter which specifies the mother vcf file path. It supports both bcf and vcf file formats.

### -ref <Reference_fasta_path>

A **required** parameter which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.

### -no-call <No Call Mode>

An **optional** parameter which sets the way program process no call variants. Default mode is none.

**implicit :** Marks both implicit and explicit noCall variants as NoCall
**explicit :** Marks only explicit noCall variants as NoCall, process implicit noCall sites as 0/0 HomRef
**none     :** Process both implicit and explicit noCall sites as 0/0 HoRef

### -outDir <Output_Directory_path>

A **required** parameter which specifies the output directory for program/error logs and ga4gh output vcf file. In current version, in order not to damage multiplatform capability, **directory should be created by user**.


### -sampleChild <child_vcf_sample_name>

An **optional** parameter which is used to select sample name from child vcf file. By default, first sample is selected.


### -sampleMother <mother_vcf_sample_name>

An **optional** parameter which is used to select sample name from mother vcf file. By default, first sample is selected.


### -sampleFather <father_vcf_sample_name>

An **optional** parameter which is used to select sample name from father vcf file. By default, first sample is selected.

### -filter <Filter_name>

An **optional** parameter which is used to filter variants with given filter name. Filter name should be same in baseline and called variants. By default, **PASS filtering** is applied to the variants. In order to disable filtering, **'-filter none'** should be used.


