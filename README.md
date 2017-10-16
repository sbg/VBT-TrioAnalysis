# VARIANT BENCHMARKING TOOLS (VBT)

**Author:** Berke Cagkan Toptas  
**Project Language:** C++11

VBT provides set of tools that is used for aligner/variant calling benchmarking. Currently there are three tools implemented:

    1. VARIANT COMPARISON
    2. MENDELIAN VIOLATION DETECTOR
    3. GRAPH VCF COMPARISON

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
./vbt graphComp [PARAMETERS]
```
For additional help, please use the following commands:

```
./vbt --help
./vbt [module_name] --help
```
Links to parameter structure of tools:

[Variant Comparison](DuoComparison/README.md)

[Mendelian Violation Detector](MendelianViolation/README.md)

[Graph Vcf Comparison](GraphComparison/README.md)

