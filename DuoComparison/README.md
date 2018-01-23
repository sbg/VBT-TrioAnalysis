
# VARIANT COMPARISON
**Author:** Berke Cagkan Toptas

Variant Comparison tool is written using the algorithm of [vcfeval](https://github.com/RealTimeGenomics/rtg-tools "vcfeval github page") using [this](http://biorxiv.org/content/early/2015/08/02/023754) paper. Please see command line parameters below.

### Parameter Format:

```
./vbt varcomp -called <called_vcf_file> -base <baseline_vcf> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

## Command line parameters:


### --help
A **single** parameter mode which prints all parameter options to the console. (./vbt varComp --help)


### -base Baseline_vcf_path

A **required** parameter which specifies the baseline vcf file path. It supports both vcf.gz and vcf file formats.


### -called Query_vcf_path

A **required** parameter which specifies the query vcf file path. It supports both vcf.gz and vcf file formats.


### -ref Reference_fasta_path

A **required** parameter which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.


### -outDir Output_Directory_path

A **required** parameter which specifies the output directory for program/error logs and ga4gh output vcf file. In current version, in order not to damage multiplatform capability, **directory should be created by user**.

### -output-mode Output_Mode

An **optional** parameter which specifies the way output will generated. **SPLIT** mode creates 4 vcf files(TPbase, TPcalled, FP, FN). **GA4GH** mode creates a single merged vcf file. Default value is **SPLIT**.

### --allele-match

An **optional** parameter which changes comparison mode to Allele Matching

### -sample-base Baseline_vcf_sample_name

An **optional** parameter which is used to select sample name from baseline vcf file. By default, first sample is selected.


### -sample-called Called_vcf_sample_name

An **optional** parameter which is used to select sample name from called vcf file. By default, first sample is selected.

### -filter Filter_name

An **optional** parameter which is used to filter variants with given filter name. Filter name should be same in baseline and called variants. By default, **PASS** filtering is applied to the variants. In order to disable filtering, **'-filter none'** should be used.

### -bed BED_file_path
An **optional** parameter which is used to select regions from VCF file.

### --snp-only

An **optional** parameter which is used to eliminate all of INDELs and SVs from input vcf files.

### --indel-only

An **optional** parameter which is used to eliminate all of SNPs and SVs from input vcf files.

### --disable-ref-overlap

An **optional** parameter which disables reference overlapping in variant comparison. Only one of overlapping variants will be placed to the True Positive variant set.

### --generate-sync-points

An **optional** parameter that prints the sync point list of two vcf file.

### --trim-endings-first
An **optional** parameter where in reference overlapping mode variants will be trimmed starting from the longest suffix first. By default, VBT is trimming prefix of alleles first.

### -thread-count [1-25]
An **optional** parameter to specify number of threads. Default value is 2

### -max-path-size Unsigned_Integer

An **optional** parameter to specify the maximum size of path that core algorithm can store inside. Default value is 150,000.

### -max-iteration-count Unsigned_Integer
An **optional** parameter to specify the maximum iteration count that core algorithm can decide to include/exclude variant. Default value is 10,000,000

### -max-bp-length Unsigned_Integer
An **optional** parameter to specify the maximum base pair length of variant to process. Default value is 1000. Variants larger than the base pair are filtered out.

