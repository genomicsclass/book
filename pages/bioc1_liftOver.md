---
title: "Translating addresses between genome builds"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---






# Translating addresses between genome builds: liftOver

The rtracklayer package includes an interface to the
liftOver utilities developed for the UCSC genome browser.
The idea is that a collection of local alignments
can be defined and used to remap coordinates from
one reference build to another.

We can illustrate this with gene addresses created for hg38,
the current reference build.  We want to translate them
for comparison to addresses asserted for hg19.

We need a "chain file", uncompressed.  You can
get it from the following URL, and use gunzip on your
system to uncompress in your home dir, if you would
like to emulate the commands below.

"ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"


```r
library(rtracklayer)
ch = import.chain("~/hg38ToHg19.over.chain")
ch
```

```
## Chain of length 25
## names(25): chr22 chr21 chr19 chr20 chrY chr18 ... chr5 chr4 chr3 chr2 chr1
```

```r
str(ch[[1]])
```

```
## Formal class 'ChainBlock' [package "rtracklayer"] with 6 slots
##   ..@ ranges  :Formal class 'IRanges' [package "IRanges"] with 6 slots
##   .. .. ..@ start          : int [1:6842] 16367189 16386933 16386970 16387001 16387128 16395491 16395528 16395841 16395860 16395956 ...
##   .. .. ..@ width          : int [1:6842] 19744 36 31 112 8362 36 312 18 95 33 ...
##   .. .. ..@ NAMES          : NULL
##   .. .. ..@ elementType    : chr "integer"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
##   ..@ offset  : int [1:6842] -480662 -480702 -480702 -480726 -480726 -480726 -480726 -480726 -480726 -480726 ...
##   ..@ score   : int [1:1168] -1063867308 68830488 21156147 20814926 7358950 3927744 2928210 991419 880681 802146 ...
##   ..@ space   : chr [1:1168] "chr22" "chr14" "chr22" "chr21" ...
##   ..@ reversed: logi [1:1168] FALSE FALSE FALSE FALSE FALSE FALSE ...
##   ..@ length  : int [1:1168] 1124 1280 173 465 398 110 43 173 342 84 ...
```

Let's get the addresses for genes on chromosome 1
in hg38.


```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx38 = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(tx38, force=TRUE) = "chr1"
g1_38 = genes(tx38)
```

```
## Warning: 'elementLengths' is deprecated.
## Use 'elementNROWS' instead.
## See help("Deprecated")

## Warning: 'elementLengths' is deprecated.
## Use 'elementNROWS' instead.
## See help("Deprecated")
```

Now execute the liftOver:


```r
g1_19L = liftOver(g1_38, ch)
```

The result is a list of GRanges, one for
each translation event.


```r
g1_19L
```

```
## Warning: 'elementLengths' is deprecated.
## Use 'elementNROWS' instead.
## See help("Deprecated")
```

```
## GRangesList object of length 2478:
## $$10000 
## GRanges object with 1 range and 1 metadata column:
##       seqnames                 ranges strand |     gene_id
##          <Rle>              <IRanges>  <Rle> | <character>
##   [1]     chr1 [243651535, 244006886]      - |       10000
## 
## $$100034743 
## GRanges object with 1 range and 1 metadata column:
##       seqnames                 ranges strand |   gene_id
##   [1]     chr1 [147466094, 147476214]      - | 100034743
## 
## $$100126331 
## GRanges object with 1 range and 1 metadata column:
##       seqnames                 ranges strand |   gene_id
##   [1]     chr1 [117637265, 117637350]      + | 100126331
## 
## ...
## <2475 more elements>
## -------
## seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

Verification of accuracy of translation is covered in exercises.
