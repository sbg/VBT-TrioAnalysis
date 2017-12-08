# VARIANT BENCHMARKING TOOLS (VBT)

**Author:** Berke Cagkan Toptas  
**Project Language:** C++11

VBT provides a set of tools that is used for aligner/variant calling benchmarking. VBT uses an improved version of [vcfeval](https://github.com/RealTimeGenomics/rtg-tools) as its core variant comparison library and uses [htslib](https://htslib.org) to read/write VCF and FASTA files. 


VBT contains 2 tool:

    1. VARIANT COMPARISON
    2. MENDELIAN VIOLATION DETECTOR

### Quick Install:

```
git clone https://gitlab.sbgdinc.com/Variants/vbt.git
cd vbt
make
```

**Note:** If you are unable to compile program or getting runtime error, please contact me. **(berke.toptas@sbgdinc.com)**

### Parameter format:

```
./vbt varComp [PARAMETERS]
./vbt mendelian [PARAMETERS]
```
For additional help, please use the following commands:

```
./vbt --help
./vbt [module_name] --help
```
Additional information can be found at each tools folder:

[Variant Comparison](DuoComparison/README.md)

[Mendelian Violation Detector](MendelianViolation/README.md)
