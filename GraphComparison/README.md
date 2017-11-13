
# GRAPH VARIANT COMPARISON
**Author:** Berke Cagkan Toptas

Graph Variant Comparison tool is written using the algorithm of [vcfeval](https://github.com/RealTimeGenomics/rtg-tools "vcfeval github page") using [this](http://biorxiv.org/content/early/2015/08/02/023754) paper. VBT uses [htslib](http://www.htslib.org/) to read/write vcf files and read fasta files. Please see command line parameters below.

### Parameter Format:

```
./vbt graphcomp -called <called_vcf_file> -base <baseline_vcf> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

## Command line parameters:


### --help
A **single** parameter mode which prints all parameter options to the console. (./vbt varComp --help)


### -base <Baseline_vcf_path>
A **required** parameter which specifies the baseline graph vcf file path. It supports both bcf and vcf file formats.


### -called <Query_vcf_path>
A **required** parameter which specifies the query graph vcf file path. It supports both bcf and vcf file formats.


### -ref <Reference_fasta_path>
A **required** parameter which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.


### -outDir <Output_Directory_path>
A **required** parameter which specifies the output directory for program/error logs 4 vcf outputs (TPbase, TPcalled, FP and FN vcfs. In current version, in order not to damage multiplatform capability, **directory should be created by user**.

### --pass-filter
An **optional** parameter which is used to filter out variants rather than **PASS**. By default all variants are included in comparison.

### -bed <BED_file_path>
An **optional** parameter which is used to select regions from graph VCF file.

### -max-path-size <Unsigned_Integer>
An **optional** parameter to specify the maximum size of path that core algorithm can store inside. Default value is 150,000.

### -max-iteration-count <Unsigned_Integer>
An **optional** parameter to specify the maximum iteration count that core algorithm can decide to include/exclude variant. Default value is 10,000,000

### -max-bp-length <Unsigned_Integer>
An **optional** parameter to specify the maximum base pair length of variant to process. Default value is 1000. Variants larger than the base pair are filtered out.

