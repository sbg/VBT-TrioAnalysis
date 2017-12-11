# VARIANT BENCHMARKING TOOLS (VBT)

**Author:** Berke Cagkan Toptas (berke.toptas@sbgdinc.com) <br /> 
**Project Language:** C++11 <br />
<br />
VBT provides a set of tools that is used for aligner/variant calling benchmarking. VBT uses [vcfeval](https://github.com/RealTimeGenomics/rtg-tools) as its core variant comparison library and [htslib](https://htslib.org) to read/write VCF and FASTA files. 
<br />
<br />
**Preprint of VBT on Biorxiv :**  TO_BE_FILLED <br/>
<br />
VBT contains 2 tool:

   1. [VARIANT COMPARISON](DuoComparison/README.md)   
   2. [MENDELIAN VIOLATION DETECTOR](MendelianViolation/README.md)

## Installing VBT:


#### Using makefile:

VBT uses htslib as an external dependency and it needs to be installed before VBT. Please download the latest htslib version from [htslib official page](https://htslib.org). After you download and extract htslib, you can install it using the following steps: (following steps are tested with htslib 1.6)

```
cd htslib-1.x
./configure
make
make install
```

You can download the latest version of VBT using the following git command:

```
git clone https://gitlab.sbgsbgd.com/cagkantoptas/VBT-TrioAnalysis   **[TO BE EDITED TO SBG GITHUB]**
```

Once you have successfully installed htslib and download VBT, Then you can use following commands to compile VBT:

```
cp -R (HTSLIB_PATH)/htslib (VBT_PATH)          //Copy htslib folder to VBT directory
cp (HTSLIB_PATH)/libhts.a (VBT_PATH)/lib       //Copy htslib binaries to VBT/lib 
cp (HTSLIB_PATH)/libhts.so (VBT_PATH)/lib      //Copy htslib binaries to VBT/lib
cp (HTSLIB_PATH)/libhts.so.2 (VBT_PATH)/lib    //Copy htslib binaries to VBT/lib
cd (VBT_PATH)                                  //Go inside the VBT folder where makefile is
make all                                       //Compile VBT
```
#### Using dockerfile:

You can run VBT as a docker image. A dockerfile is added to the project that is running on Ubuntu 14.04 instance.

```
cd (VBT_PATH)
docker build -t vbtapp:v1.0 .
docker run -ti vbtapp:v1.0
```

**Note:** If you are unable to compile program or getting runtime error, please contact me. **(berke.toptas@sbgdinc.com)**

## Parameter format:

```
./vbt varComp [PARAMETERS]
./vbt mendelian [PARAMETERS]
```
For additional help, please use the following commands:

```
./vbt --help
./vbt [module_name] --help
```

## Sample Execution:

Example folder contains a minimal example for testing VBT. We placed chr21 of merged BWA+UG CEPH trio. For the FASTA file, we extract chr21 from GRCh37 reference. In order to run the sample test, following command can be executed :

```
./vbt -mother Example/UG_CEU_merged_cleaned_chr21.vcf -father Example/UG_CEU_merged_cleaned_chr21.vcf -child Example/UG_CEU_merged_cleaned_chr21.vcf -ped Example/ceu.ped --ref-overlap --outDir Example/output -out-prefix ceu_21_sample
```
This sample execution will produce:

	1. ceu_21_sample_BestPathLogs.txt              //Logs of Best Path Algorithm
	2. ceu_21_sample_ChildReportLog.txt            //Mendelian decision for Non-Ref child variants only
	3. ceu_21_sample_DetailedLogs.txt              //Detailed human-readable logs of trio concordance analysis
	4. ceu_21_sample_tab_delim_detailed_log.tsv    //Tab delimited version of detailed logs
	5. ceu_21_sample_trio.vcf                      //Merged output trio (Mendelian decisions are annotated for each record)
	
For additional parameters, please check paramater specifications.

## Source Code Documentation:

A [Doxygen html documentation](Doxygen) is available for VBT under Doxygen folder.

