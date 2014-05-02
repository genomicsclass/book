---
layout: page
title: Basic Bioconductor infrastructure
---




## IRanges


```r
library(IRanges)
```

```
## Loading required package: methods
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
```

```r
ir <- IRanges(5, 10)
ir
```

```
## IRanges of length 1
##     start end width
## [1]     5  10     6
```

```r
start(ir)
```

```
## [1] 5
```

```r
end(ir)
```

```
## [1] 10
```

```r
width(ir)
```

```
## [1] 6
```

```r
# ?IRanges
```



```r
ir <- IRanges(start = c(3, 5, 17), end = c(10, 8, 20))
ir
```

```
## IRanges of length 3
##     start end width
## [1]     3  10     8
## [2]     5   8     4
## [3]    17  20     4
```

```r
ir <- IRanges(5, 10)
```



```r
# ?'intra-range-methods'
shift(ir, -2)
```

```
## IRanges of length 1
##     start end width
## [1]     3   8     6
```


Remeber, all of these commands can work on more than one range at once. Here we show the effects of the different methods using a single range:


```r
shift(ir, -2)
```

```
## IRanges of length 1
##     start end width
## [1]     3   8     6
```

```r
narrow(ir, start = 2)
```

```
## IRanges of length 1
##     start end width
## [1]     6  10     5
```

```r
narrow(ir, end = 5)
```

```
## IRanges of length 1
##     start end width
## [1]     5   9     5
```

```r
flank(ir, width = 3, start = TRUE, both = FALSE)
```

```
## IRanges of length 1
##     start end width
## [1]     2   4     3
```

```r
flank(ir, width = 3, start = FALSE, both = FALSE)
```

```
## IRanges of length 1
##     start end width
## [1]    11  13     3
```

```r
flank(ir, width = 3, start = TRUE, both = TRUE)
```

```
## IRanges of length 1
##     start end width
## [1]     2   7     6
```

```r
ir * 2
```

```
## IRanges of length 1
##     start end width
## [1]     6   8     3
```

```r
ir + 2
```

```
## IRanges of length 1
##     start end width
## [1]     3  12    10
```

```r
ir - 2
```

```
## IRanges of length 1
##     start end width
## [1]     7   8     2
```




```r
# set up a plotting window so we can look at range operations
plotir <- function(ir, i) {
    arrows(start(ir) - 0.5, i, end(ir) + 0.5, i, code = 3, angle = 90, lwd = 3)
}
plot(0, 0, xlim = c(0, 15), ylim = c(0, 11), type = "n", xlab = "", ylab = "", 
    xaxt = "n")
axis(1, 0:15)
abline(v = 0:30 + 0.5, col = rgb(0, 0, 0, 0.5))

# plot the original IRange
plotir(ir, 1)

# draw a red shadow for the original IRange
polygon(c(start(ir) - 0.5, start(ir) - 0.5, end(ir) + 0.5, end(ir) + 0.5), c(-1, 
    12, 12, -1), col = rgb(1, 0, 0, 0.2), border = NA)
plotir(shift(ir, -2), 2)
plotir(narrow(ir, start = 2), 3)
plotir(narrow(ir, end = 5), 4)
plotir(flank(ir, width = 3, start = TRUE, both = FALSE), 5)
plotir(flank(ir, width = 3, start = FALSE, both = FALSE), 6)
plotir(flank(ir, width = 3, start = TRUE, both = TRUE), 7)
plotir(ir * 2, 8)
plotir(ir + 2, 9)
plotir(ir - 2, 10)
```

![plot of chunk unnamed-chunk-5](figure/basic_Bioconductor_infrastructure-unnamed-chunk-5.png) 



```r
# ?'inter-range-methods'
ir <- IRanges(start = c(3, 5, 17), end = c(10, 8, 20))
range(ir)
```

```
## IRanges of length 1
##     start end width
## [1]     3  20    18
```

```r
reduce(ir)
```

```
## IRanges of length 2
##     start end width
## [1]     3  10     8
## [2]    17  20     4
```

```r
gaps(ir)
```

```
## IRanges of length 1
##     start end width
## [1]    11  16     6
```

```r
disjoin(ir)
```

```
## IRanges of length 4
##     start end width
## [1]     3   4     2
## [2]     5   8     4
## [3]     9  10     2
## [4]    17  20     4
```


## GRanges and GRangesList

### GRanges


```r
library(GenomicRanges)
```

```
## Loading required package: GenomeInfoDb
```

```r
gr <- GRanges("chrZ", IRanges(start = c(5, 10), end = c(35, 45)), strand = "+", 
    seqlengths = c(chrZ = 100L))
gr
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [ 5, 35]      +
##   [2]     chrZ  [10, 45]      +
##   ---
##   seqlengths:
##    chrZ
##     100
```

```r
shift(gr, 10)
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [15, 45]      +
##   [2]     chrZ  [20, 55]      +
##   ---
##   seqlengths:
##    chrZ
##     100
```

```r
shift(gr, 80)
```

```
## Warning: 'ranges' contains values outside of sequence bounds. See ?trim to
## subset ranges.
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ [85, 115]      +
##   [2]     chrZ [90, 125]      +
##   ---
##   seqlengths:
##    chrZ
##     100
```

```r
trim(shift(gr, 80))
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ [85, 100]      +
##   [2]     chrZ [90, 100]      +
##   ---
##   seqlengths:
##    chrZ
##     100
```

```r
mcols(gr)
```

```
## DataFrame with 2 rows and 0 columns
```

```r
mcols(gr)$value <- c(-1, 4)
gr
```

```
## GRanges with 2 ranges and 1 metadata column:
##       seqnames    ranges strand |     value
##          <Rle> <IRanges>  <Rle> | <numeric>
##   [1]     chrZ  [ 5, 35]      + |        -1
##   [2]     chrZ  [10, 45]      + |         4
##   ---
##   seqlengths:
##    chrZ
##     100
```


### GRangesList


```r
gr2 <- GRanges("chrZ", IRanges(11:13, 51:53))
mcols(gr)$value <- NULL
grl <- GRangesList(gr, gr2)
grl
```

```
## GRangesList of length 2:
## [[1]] 
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [ 5, 35]      +
##   [2]     chrZ  [10, 45]      +
## 
## [[2]] 
## GRanges with 3 ranges and 0 metadata columns:
##       seqnames   ranges strand
##   [1]     chrZ [11, 51]      *
##   [2]     chrZ [12, 52]      *
##   [3]     chrZ [13, 53]      *
## 
## ---
## seqlengths:
##  chrZ
##   100
```

```r
length(grl)
```

```
## [1] 2
```

```r
grl[[1]]
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [ 5, 35]      +
##   [2]     chrZ  [10, 45]      +
##   ---
##   seqlengths:
##    chrZ
##     100
```

```r
mcols(grl)$value <- c(5, 7)
grl
```

```
## GRangesList of length 2:
## [[1]] 
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [ 5, 35]      +
##   [2]     chrZ  [10, 45]      +
## 
## [[2]] 
## GRanges with 3 ranges and 0 metadata columns:
##       seqnames   ranges strand
##   [1]     chrZ [11, 51]      *
##   [2]     chrZ [12, 52]      *
##   [3]     chrZ [13, 53]      *
## 
## ---
## seqlengths:
##  chrZ
##   100
```

```r
mcols(grl)
```

```
## DataFrame with 2 rows and 1 column
##       value
##   <numeric>
## 1         5
## 2         7
```


### findOverlaps and %over%


```r
gr1 <- GRanges("chrZ", IRanges(c(1, 11, 21, 31, 41), width = 5))
gr2 <- GRanges("chrZ", IRanges(c(19, 33), c(38, 35)))
gr1
```

```
## GRanges with 5 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [ 1,  5]      *
##   [2]     chrZ  [11, 15]      *
##   [3]     chrZ  [21, 25]      *
##   [4]     chrZ  [31, 35]      *
##   [5]     chrZ  [41, 45]      *
##   ---
##   seqlengths:
##    chrZ
##      NA
```

```r
gr2
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [19, 38]      *
##   [2]     chrZ  [33, 35]      *
##   ---
##   seqlengths:
##    chrZ
##      NA
```

```r
fo <- findOverlaps(gr1, gr2)
fo
```

```
## Hits of length 3
## queryLength: 5
## subjectLength: 2
##   queryHits subjectHits 
##    <integer>   <integer> 
##  1         3           1 
##  2         4           1 
##  3         4           2
```

```r
queryHits(fo)
```

```
## [1] 3 4 4
```

```r
subjectHits(fo)
```

```
## [1] 1 1 2
```

```r
gr1 %over% gr2
```

```
## [1] FALSE FALSE  TRUE  TRUE FALSE
```

```r
gr1[gr1 %over% gr2]
```

```
## GRanges with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrZ  [21, 25]      *
##   [2]     chrZ  [31, 35]      *
##   ---
##   seqlengths:
##    chrZ
##      NA
```


### Rle and Views


```r
r <- Rle(c(1, 1, 1, 0, 0, -2, -2, -2, rep(-1, 20)))
r
```

```
## numeric-Rle of length 28 with 4 runs
##   Lengths:  3  2  3 20
##   Values :  1  0 -2 -1
```

```r
str(r)
```

```
## Formal class 'Rle' [package "IRanges"] with 4 slots
##   ..@ values         : num [1:4] 1 0 -2 -1
##   ..@ lengths        : int [1:4] 3 2 3 20
##   ..@ elementMetadata: NULL
##   ..@ metadata       : list()
```

```r
as.numeric(r)
```

```
##  [1]  1  1  1  0  0 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
## [24] -1 -1 -1 -1 -1
```

```r
Views(r, start = c(4, 2), end = c(7, 6))
```

```
## Views on a 28-length Rle subject
## 
## views:
##     start end width
## [1]     4   7     4 [ 0  0 -2 -2]
## [2]     2   6     5 [ 1  1  0  0 -2]
```




## ExpressionSet and SummarizedExperiment


```r
library(Biobase)
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
library(GEOquery)
```

```
## Error: there is no package called 'GEOquery'
```

```r
geoq <- getGEO("GSE9514")
```

```
## Error: could not find function "getGEO"
```

```r
names(geoq)
```

```
## Error: object 'geoq' not found
```

```r
e <- geoq[[1]]
```

```
## Error: object 'geoq' not found
```


### ExpressionSet


```r
dim(e)
```

```
## Error: object 'e' not found
```

```r
exprs(e)[1:3, 1:3]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'e' not found
```

```r
dim(exprs(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'e' not found
```

```r

pData(e)[1:3, 1:6]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
dim(pData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
names(pData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
pData(e)$characteristics_ch1
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r

fData(e)[1:3, 1:3]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
dim(fData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
names(fData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
head(fData(e)$"Gene Symbol")
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in fData(e) : 
##   error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
head(rownames(e))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in rownames(e) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'e' not found
```

```r

experimentData(e)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'experimentData': Error: object 'e' not found
```

```r
annotation(e)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'annotation': Error: object 'e' not found
```


### Summarized Experiment


```r
library(parathyroidSE)
```

```
## Error: there is no package called 'parathyroidSE'
```

```r
data(parathyroidGenesSE)
```

```
## Warning: data set 'parathyroidGenesSE' not found
```

```r
se <- parathyroidGenesSE
```

```
## Error: object 'parathyroidGenesSE' not found
```

```r
se
```

```
## Error: object 'se' not found
```




```r
dim(se)
```

```
## Error: object 'se' not found
```

```r
assay(se)[1:3, 1:3]
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'se' not found
```

```r
dim(assay(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'se' not found
```

```r

colData(se)[1:3, 1:6]
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'colData': Error: object 'se' not found
```

```r
dim(colData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'colData': Error: object 'se' not found
```

```r
names(colData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'colData': Error: object 'se' not found
```

```r
colData(se)$treatment
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'colData': Error: object 'se' not found
```

```r

rowData(se)[1]
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'rowData': Error: object 'se' not found
```

```r
class(rowData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'rowData': Error: object 'se' not found
```

```r
length(rowData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'rowData': Error: object 'se' not found
```

```r
head(rownames(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in rownames(se) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'se' not found
```

```r
metadata(rowData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'metadata': Error in rowData(se) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rowData': Error: object 'se' not found
```

```r

exptData(se)$MIAME
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'exptData': Error: object 'se' not found
```

```r
abstract(exptData(se)$MIAME)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'abstract': Error in exptData(se) : 
##   error in evaluating the argument 'x' in selecting a method for function 'exptData': Error: object 'se' not found
```

