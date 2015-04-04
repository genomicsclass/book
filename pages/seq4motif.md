---
title: "Genomic sequence -- utility for motif checking"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---





# Overview

In this document we'll show how to look for occurrences
of a binding motif in genomic sequence underlying binding peaks.

Recall that we have the ER binding peaks for two cell
lines in the ERBS package.  We'll focus on HepG2


```r
library(ERBS)
data(HepG2)
HepG2
```

```
## GRanges object with 303 ranges and 7 metadata columns:
##         seqnames                 ranges strand   |      name     score
##            <Rle>              <IRanges>  <Rle>   | <numeric> <integer>
##     [1]     chr2 [ 20335378,  20335787]      *   |      <NA>         0
##     [2]    chr20 [   328285,    329145]      *   |      <NA>         0
##     [3]     chr6 [168135432, 168136587]      *   |      <NA>         0
##     [4]    chr19 [  1244419,   1245304]      *   |      <NA>         0
##     [5]    chr11 [ 64071828,  64073069]      *   |      <NA>         0
##     ...      ...                    ...    ... ...       ...       ...
##   [299]     chr4 [  1797182,   1797852]      *   |      <NA>         0
##   [300]     chr1 [110198573, 110199126]      *   |      <NA>         0
##   [301]    chr17 [ 17734052,  17734469]      *   |      <NA>         0
##   [302]     chr1 [ 48306453,  48306908]      *   |      <NA>         0
##   [303]    chr12 [123867207, 123867554]      *   |      <NA>         0
##               col signalValue    pValue       qValue      peak
##         <logical>   <numeric> <numeric>    <numeric> <integer>
##     [1]      <NA>      58.251    75.899 6.143712e-72       195
##     [2]      <NA>      10.808    69.685 5.028065e-66       321
##     [3]      <NA>      17.103    54.311 7.930665e-51       930
##     [4]      <NA>      12.427    43.855 1.359756e-40       604
##     [5]      <NA>       10.85    40.977 7.333863e-38       492
##     ...       ...         ...       ...          ...       ...
##   [299]      <NA>       9.681    10.057 1.423343e-08       402
##   [300]      <NA>       7.929    10.047 1.442076e-08       197
##   [301]      <NA>       5.864      9.99 1.638918e-08       227
##   [302]      <NA>        5.66     9.948 1.799414e-08       211
##   [303]      <NA>      13.211     9.918 1.921805e-08       163
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
```

We'd like to look at the genomic sequence underneath the peaks
and inspect it for the binding motif "TCAAGGTCA".  This is
easy to do with the Biostrings and BSGenome infrastructure.

# Reference genomic sequence for humans

We'll work with hg19.  The BSgenome... package will
create variable `Hsapiens` on attachment.
This variable gives a metadata report.


```r
library(BSgenome.Hsapiens.UCSC.hg19)
```

```
## Loading required package: BSgenome
## Loading required package: Biostrings
## Loading required package: XVector
## Loading required package: rtracklayer
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

The reference sequence for a chromosome can be obtained
with the $$ operator.


```r
Hsapiens$chr17
```

```
##   81195210-letter "DNAString" instance
## seq: AAGCTTCTCACCCTGTTCCTGCATAGATAATTGC...GGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGT
```

# Targeted retrieval of reference sequence

The getSeq function obtains sequence
corresponding to addresses listed in GRanges.
We'll obtain the sequence under the peaks as
`hepseq`, and a set of control sequences of similar
lengths obtained by shifting the binding peak intervals
by 2500 bases and obtaining the reference sequence in the
shifted intervals.


```r
hepseq = getSeq(Hsapiens, HepG2)
rhepseq = getSeq(Hsapiens, shift(HepG2,2500))
hepseq
```

```
##   A DNAStringSet instance of length 303
##       width seq
##   [1]   410 GAGACAGGGTTTCACCATGTTGGCCAGGCT...TCCAGGAAGCAGAAATGTTCAAGGACTCTC
##   [2]   861 TGGGAAGGACACAACTGAATGAGGCTGTGC...AGAACCTCCAACCGTGTGTGTGTGTGTGTG
##   [3]  1156 GACACCTGCCACCCCGGACCCCACAGAATG...TCGTGTCTGCTTTCTTATGTGTTTTTGTTT
##   [4]   886 GTGAAGGCCCTGGAGTAGGCGGTGCGTACC...GTTTTTGGCACCTCCGTGGGCACCTAGGCT
##   [5]  1242 CATCCTCCACCTTAACACTCAGCACCCTTA...TGTGTCCTACAAGCAGCCGGCGGCGCCGCC
##   ...   ... ...
## [299]   671 CACTGGAGCTGGTGAAACAGGTAGTGAGTT...ATCTAGGGAGGCATGCAGCCCTCACCTGAG
## [300]   554 TCCGGAGAAGAAGAAACGGGGGAAGAACTT...GCCGAGCGGCTGGGGACCGGCTCTAGGGAC
## [301]   418 CCACACCTGGAGCCAGTCTCAATGGCTCCC...TGAATGGTTGGAGACCAGGGGAGTTCTGTG
## [302]   456 TGATGACATTTCTCAAGGATTAAGAAAAAG...TGCACCCATTTTGGTTTTGCTGTAGGGCCT
## [303]   348 TCCAAAGCAGACACTCCAGGACACCTGATT...CTTTTTTTGAGACGGAGTCTCGTTCTGTCG
```

# Counting motif occurrences

We count the occurrences of the ESRRA
binding motif
"TCAAGGTCA" in the bound intervals (and their reverse complement
representation).  This is compared to the frequency of occurrence in the
control sequences.  We'll use the `vcountPattern` function of the
Biostrings package to carry this out.


```r
sum(vcountPattern("TCAAGGTCA", hepseq))+sum(vcountPattern("TCAAGGTCA", 
   reverseComplement(hepseq)))
```

```
## [1] 55
```

```r
sum(vcountPattern("TCAAGGTCA", rhepseq))+sum(vcountPattern("TCAAGGTCA", 
   reverseComplement(rhepseq)))
```

```
## [1] 6
```

We see a 9-fold increase in occupancy in the bound regions compared
to the shifted regions.  This is not the way one assesses motif occurrences.
First, the motif is generally represented as a model and not a string.
The model is typically expressed as a position weight matrix (PWM).
Second, the most common software tools for evaluating motif enrichment are
MEME and FIMO; matchPWM of the Biostrings package can perform similar analyses.
package can also 
