---
title: "Chromosomes and their substructures 1: Reference genomes"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---





# Genomic sequence, reference builds

## Human

The genomic sequence for humans has recently
been revised.  We can use the most recent major
revision as follows:


```r
library(BSgenome.Hsapiens.NCBI.GRCh38)
Hsapiens
```

```
## Human genome
## | 
## | organism: Homo sapiens (Human)
## | provider: NCBI
## | provider version: GRCh38
## | release date: 2013-12-17
## | release name: Genome Reference Consortium Human Build 38
## | 455 sequences:
## |   1                                  2                                 
## |   3                                  4                                 
## |   5                                  6                                 
## |   7                                  8                                 
## |   9                                  10                                
## |   ...                                ...                               
## |   HSCHR19KIR_FH05_B_HAP_CTG3_1       HSCHR19KIR_FH06_A_HAP_CTG3_1      
## |   HSCHR19KIR_FH06_BA1_HAP_CTG3_1     HSCHR19KIR_FH08_A_HAP_CTG3_1      
## |   HSCHR19KIR_FH08_BAX_HAP_CTG3_1     HSCHR19KIR_FH13_A_HAP_CTG3_1      
## |   HSCHR19KIR_FH13_BA2_HAP_CTG3_1     HSCHR19KIR_FH15_A_HAP_CTG3_1      
## |   HSCHR19KIR_RP5_B_HAP_CTG3_1                                          
## | (use 'seqnames()' to see all the sequence names, use the '$' or '[['
## | operator to access a given sequence)
```

```r
h38 = Hsapiens # for later
```

Notice the number of sequences reported, and their names.  We can
get the sequence for a chromosome by using list-like
syntax with `Hsapiens`.


```r
h38$"22"
```

```
##   50818468-letter "DNAString" instance
## seq: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

This shows that the starting and ending bases are indeterminate.
We can obtain the overall nucleotide frequencies as


```r
alphabetFrequency(Hsapiens$"22")
```

```
##        A        C        G        T        M        R        W        S 
## 10382214  9160652  9246186 10370725        0        2        1        0 
##        Y        K        V        H        D        B        N        - 
##        2        0        0        0        0        0 11658686        0 
##        +        . 
##        0        0
```

A great deal of reference data in use are annotated to 
build hg19 (also known as GRCh37).


```r
library(BSgenome.Hsapiens.UCSC.hg19)
```

```
## 
## Attaching package: 'BSgenome.Hsapiens.UCSC.hg19'
## 
## The following object is masked from 'package:BSgenome.Hsapiens.NCBI.GRCh38':
## 
##     Hsapiens
```

```r
Hsapiens
```

```
## Human genome
## | 
## | organism: Homo sapiens (Human)
## | provider: UCSC
## | provider version: hg19
## | release date: Feb. 2009
## | release name: Genome Reference Consortium GRCh37
## | 93 sequences:
## |   chr1                  chr2                  chr3                 
## |   chr4                  chr5                  chr6                 
## |   chr7                  chr8                  chr9                 
## |   chr10                 chr11                 chr12                
## |   chr13                 chr14                 chr15                
## |   ...                   ...                   ...                  
## |   chrUn_gl000235        chrUn_gl000236        chrUn_gl000237       
## |   chrUn_gl000238        chrUn_gl000239        chrUn_gl000240       
## |   chrUn_gl000241        chrUn_gl000242        chrUn_gl000243       
## |   chrUn_gl000244        chrUn_gl000245        chrUn_gl000246       
## |   chrUn_gl000247        chrUn_gl000248        chrUn_gl000249       
## | (use 'seqnames()' to see all the sequence names, use the '$' or '[['
## | operator to access a given sequence)
```

```r
h19 = Hsapiens
```

Note that there is a different sequence naming convention
and a different number of sequences managed in this build.


## Other organisms

If you have an internet connection, the `available.genomes` function
will list packages that contain reference sequences.


```r
available.genomes()
```
 
For organisms not covered at present by the project, tools
for building compatible packages are available in the
BSgenome package (see the BSgenomeForge vignette).


