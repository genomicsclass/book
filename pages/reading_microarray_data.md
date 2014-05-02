---
layout: page
title: Reading in microarray data
---




## Affymterix CEL files

We start by reading in the sample information table. This is usually created by the person who performed the experiment. 


```r
library(affy)
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
##     table, tapply, union, unique, unlist
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
basedir <- "celfiles"
tab <- read.delim(file.path(basedir, "sampleinfo.txt"), check.names = FALSE, 
    as.is = TRUE)
fns <- list.celfiles(basedir)
fns %in% tab[, 1]  ##check
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
ab <- ReadAffy(filenames = file.path(basedir, tab[, 1]), phenoData = tab)
```

```
## Warning: Mismatched phenoData and celfile names!
## 
## Please note that the row.names of your phenoData object should be identical to what you get from list.celfiles()!
## Otherwise you are responsible for ensuring that the ordering of your phenoData object conforms to the ordering of the celfiles as they are read into the AffyBatch!
## If not, errors may result from using the phenoData for subsetting or creating linear models, etc.
```


This creates an AffyBatch object which object contains infomration you need.


```r
dim(pm(ab))
```

```
## 
## The downloaded source packages are in
## 	'/private/var/folders/6d/d_8pbllx7318htlp5wv_rm580000gn/T/RtmpjcHJ2q/downloaded_packages'
```

```
## Warning: replacing previous import by 'utils::head' when loading 'hgu95acdf'
## Warning: replacing previous import by 'utils::tail' when loading 'hgu95acdf'
```

```
## 
```

```
## [1] 201807      6
```

```r
dim(pData(ab))
```

```
## [1]  6 17
```

```r
annotation(ab)
```

```
## [1] "hgu95a"
```


Note, this object You can then preprocess RMA

```r
e <- rma(ab)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```


If you are not interested in probe level data you could can use this function

```r
e <- justRMA(filenames = tab[, 1], celfile.path = basedir, phenoData = tab)
```

```
## Warning: Mismatched phenoData and celfile names!
## 
## Please note that the row.names of your phenoData object should be identical to what you get from list.celfiles()!
## Otherwise you are responsible for ensuring that the ordering of your phenoData object conforms to the ordering of the celfiles as they are read into the AffyBatch!
## If not, errors may result from using the phenoData for subsetting or creating linear models, etc.
```


##Agilent data


```r
library(limma)
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
basedir <- "agilent"
targets <- readTargets(file.path(basedir, "TargetBeta7.txt"))
RG <- read.maimages(targets$FileName, source = "genepix", path = basedir)
```

```
## Warning: Name partially matched in data frame
```

```
## Read agilent/6Hs.195.1.gpr 
## Read agilent/6Hs.168.gpr 
## Read agilent/6Hs.166.gpr 
## Read agilent/6Hs.187.1.gpr 
## Read agilent/6Hs.194.gpr 
## Read agilent/6Hs.243.1.gpr
```

```r
MA <- MA.RG(RG, bc.method = "none")
mypar(1, 1)
imageplot(MA$M[, 2], RG$printer, zlim = c(-3, 3))
```

![plot of chunk unnamed-chunk-5](figure/reading_microarray_data-unnamed-chunk-5.png) 

```r
dev.off()
```

```
## null device 
##           1
```






## oligo
We can also use oligo to read affy arrays


```r
detach("package:affy")
library(oligo)
```

```
## Loading required package: oligoClasses
## Welcome to oligoClasses version 1.26.0
## Loading required package: Biostrings
## Loading required package: IRanges
## Loading required package: XVector
## ===========================================================================
## Welcome to oligo version 1.28.0
## ===========================================================================
## 
## Attaching package: 'oligo'
## 
## The following object is masked from 'package:limma':
## 
##     backgroundCorrect
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     normalize
```

```r
basedir <- "celfiles"
tab <- read.delim(file.path(basedir, "sampleinfo.txt"), check.names = FALSE, 
    as.is = TRUE)
fns <- list.celfiles(basedir, listGzipped = TRUE)
fns %in% tab[, 1]  ##check
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
pd <- as(tab, "AnnotatedDataFrame")
efs <- read.celfiles(filenames = file.path(basedir, tab[, 1]), phenoData = pd, 
    sampleNames = sampleNames(pd))
```

```
## Loading required package: pd.hg.u95a
```

```
## Warning: there is no package called 'pd.hg.u95a'
```

```
## Attempting to obtain 'pd.hg.u95a' from BioConductor website.
## Checking to see if your internet connection works...
```

```
## 
## The downloaded source packages are in
## 	'/private/var/folders/6d/d_8pbllx7318htlp5wv_rm580000gn/T/RtmpjcHJ2q/downloaded_packages'
```

```
## Loading required package: pd.hg.u95a
## Loading required package: RSQLite
## Loading required package: DBI
## Platform design info loaded.
```

```
## Reading in : celfiles/1521a99hpp_av06.CEL.gz
## Reading in : celfiles/1532a99hpp_av04.CEL.gz
## Reading in : celfiles/2353a99hpp_av08.CEL.gz
## Reading in : celfiles/1521b99hpp_av06.CEL.gz
## Reading in : celfiles/1532b99hpp_av04.CEL.gz
## Reading in : celfiles/2353b99hpp_av08r.CEL.gz
```

```
## Warning: 'channel' automatically added to varMetadata in phenoData.
```



```r
e <- rma(efs)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```

