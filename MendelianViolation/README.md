# MENDELIAN VIOLATION DETECTOR
**Author:** Berke Cagkan Toptas

A tool for checking all mendelian violations in parent-child data.

### Parameter Format:

```
./vbt mendelian -father <father_vcf_path> -mother <mother_vcf_path> -child <child_vcf_path> -ref <reference_fasta> -outDir <output_directory> [OPTIONAL PARAMETERS]
```

## Command line parameters:
**Warning:** Parameters are **case sensitive**.


### --help
A **single** parameter which prints all parameter options to the console. (./vbt mendelian --help)


### -father <father_vcf_path>

A **required** parameter which specifies the father vcf file path. It supports both bcf and vcf file formats. If given vcf contain multiple samples, first sample will be taken by default unless a Pedigree file is provided (see -pedigree command)


### -mother <mother_vcf_path>

A **required** parameter which specifies the mother vcf file path. It supports both bcf and vcf file formats. If given vcf contain multiple samples, first sample will be taken by default unless a Pedigree file is provided (see -pedigree command)

### -child <child_vcf_path>

A **required** parameter which specifies the mother vcf file path. It supports both bcf and vcf file formats. If given vcf contain multiple samples, first sample will be taken by default unless a Pedigree file is provided (see -pedigree command)

### -ref <Reference_fasta_path>

A **required** parameter which specifies the reference FASTA file path. It should be in FASTA (.fa) format. A FASTA index file is not mandatory. The tool will automatically generate a FASTA index file (.fai) if it does not exist.

### -no-call <No Call Mode>

An **optional** parameter which sets the way program process no call variants. Default mode is none.

**implicit :** Marks both implicit and explicit noCall variants as NoCall  
**explicit :** Marks only explicit noCall variants as NoCall, process implicit noCall sites as 0/0 HomRef  
**none     :** Process both implicit and explicit noCall sites as 0/0 HoRef

### -outDir <Output_Directory_path>

A **required** parameter which specifies the output directory for program/error logs and ga4gh output vcf file. In current version, in order not to damage multiplatform capability, **directory should be created by user**.


### -pedigree <PED_File_path>


An **optional** parameter which is used to determine all 3 samples from a PED file. Is useful where input vcf files contain more than 1 sample

### --autosome-only

An **optional** parameter which is used to process only autosomes [chr1:chr22] By default, it is set to FALSE. Note that in order to use this parameter, the chromosome naming should be either 1,2,3 etc. or chr1, chr2, chr3. For other chromosome namings, please use -bed command to specify autosomes.

### -bed <BED_file_path>

An **optional** parameter which is used to select regions from VCF file.


### -filter <Filter_name>

An **optional** parameter which is used to filter variants with given filter name. Filter name should be same in baseline and called variants. By default, **PASS filtering** is applied to the variants. In order to disable filtering, **'-filter none'** should be used.


