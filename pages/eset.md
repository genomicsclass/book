---
title: "The ExpressionSet container"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
---



# Overview

We'll work with the basic representation of expression experiments
in Bioconductor.  An example is in package Biobase.


```r
library(Biobase)
```

```
## Loading required package: BiocGenerics
## Loading required package: methods
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
##     table, tapply, union, unique, unlist, unsplit
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
data(sample.ExpressionSet)
sample.ExpressionSet
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 500 features, 26 samples 
##   element names: exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: A B ... Z (26 total)
##   varLabels: sex type score
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: hgu95av2
```

We'll abbreviate the name:


```r
samp = sample.ExpressionSet
```

# Queries and extractors


```r
dim(samp)
```

```
## Features  Samples 
##      500       26
```

```r
exprs(samp)[1:5,1:6]  # extract expression values
```

```
##                        A         B        C        D        E       F
## AFFX-MurIL2_at  192.7420  85.75330 176.7570 135.5750 64.49390 76.3569
## AFFX-MurIL10_at  97.1370 126.19600  77.9216  93.3713 24.39860 85.5088
## AFFX-MurIL4_at   45.8192   8.83135  33.0632  28.7072  5.94492 28.2925
## AFFX-MurFAS_at   22.5445   3.60093  14.6883  12.3397 36.86630 11.2568
## AFFX-BioB-5_at   96.7875  30.43800  46.1271  70.9319 56.17440 42.6756
```

```r
pData(samp)    # extract sample level data
```

```
##      sex    type score
## A Female Control  0.75
## B   Male    Case  0.40
## C   Male Control  0.73
## D   Male    Case  0.42
## E Female    Case  0.93
## F   Male Control  0.22
## G   Male    Case  0.96
## H   Male    Case  0.79
## I Female    Case  0.37
## J   Male Control  0.63
## K   Male    Case  0.26
## L Female Control  0.36
## M   Male    Case  0.41
## N   Male    Case  0.80
## O Female    Case  0.10
## P Female Control  0.41
## Q Female    Case  0.16
## R   Male Control  0.72
## S   Male    Case  0.17
## T Female    Case  0.74
## U   Male Control  0.35
## V Female Control  0.77
## W   Male Control  0.27
## X   Male Control  0.98
## Y Female    Case  0.94
## Z Female    Case  0.32
```

```r
experimentData(samp)
```

```
## Experiment data
##   Experimenter name: Pierre Fermat 
##   Laboratory: Francis Galton Lab 
##   Contact information: pfermat@lab.not.exist 
##   Title: Smoking-Cancer Experiment 
##   URL: www.lab.not.exist 
##   PMIDs:  
## 
##   Abstract: A 8 word abstract is available. Use 'abstract' method.
##   notes:
##    notes:     
##       An example object of expression set (exprSet) class
```

```r
abstract(samp)  # special accessor
```

```
## [1] "An example object of expression set (ExpressionSet) class"
```

Have a look at annotation package pmid2MIAME function to see
how to extract abstracts of papers from pubmed.  These can be
bound into ExpressionSets with experimentData().

# Matrix-like subscripting

We can use matrix-like syntax directly to restrict the
ExpressionSet, getting back a new ExpressionSet

```r
samp[1:4,3:20]
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 4 features, 18 samples 
##   element names: exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: C D ... T (18 total)
##   varLabels: sex type score
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: hgu95av2
```

```r
samp[, samp$sex=="Male"]
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 500 features, 15 samples 
##   element names: exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: B C ... X (15 total)
##   varLabels: sex type score
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: hgu95av2
```


