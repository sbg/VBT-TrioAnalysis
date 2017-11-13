# VBT CORE LIBRARY

This folder contains comparison library originated from [RTG Tools vcfeval](https://github.com/RealTimeGenomics/rtg-tools/tree/master/src/com/rtg/vcf/eval). We have made slight changes on the original implementation to improve performance and accuracy.

```

CHR      POS              REF                    ALT                BASELINE   Query
1        4582981          ATCTATCATCTCTCTCTC     A           GT       0/1       0/1
1        4582988          ATCTCTCTC              ATCTC       GT       0/1       0/1
1        4582994          C                      A           GT       0/1       0/1
```


For example, in ref overlap mode, 3 True Positive would expected for above variants between baseline and query variants. However the third variant is excluded by vcfeval.

```
4582988      ATCTCTCTC   ATCTC              Range : 4582988 - 4582997
4582988      *****TCTC   *****      Trimmed Range : 4582993 - 4582997
```
In -ref-overlap mode, vcfeval try to trim variants in a lenient way to solve reference overlapping. After vcfeval does the trimming as seen above, third variant still resides in range of second variant (4582994). But if vcfeval would do the trimming as below, the third variant would be included in TP set without a problem:

```
4582988      *TCTC****  *****      Trimmed Range : 4582989 - 4582993
or
4582988      **CTCT***  *****      Trimmed Range : 4582990 - 4582994
```

We implemented a smart trimming method to catch these variants.

