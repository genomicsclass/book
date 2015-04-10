---
layout: page
title: "Architecture: considerations on high performance computing with Bioconductor"
---




## Overview of performance enhancements

There are two main approaches to achieving scalability of
analytical workflows in Bioconductor.

- Shared memory parallelism.  The R process is forked an
arbitrary number of times with full copies of the memory
image and computations proceed independently for each image.
Selected results are returned to the master process.  This
is often effective in multicore environments.

- Distributed parallelism.  Independent processors, potentially
running different operating systems, can run compatible instances
of R.  Job control can be carried out by R or by a cluster scheduler.

For tasks that are "embarrassingly parallel", that do not require
communication between processes beyond starting, stopping, and
returning results, either of these approaches can be used
reasonably simply in R.

### Simple illustration

We can demonstrate the shared memory approach with our laptop.


```r
system.time( lapply(1:8, function(x)Sys.sleep(1) ) )
```

```
##    user  system elapsed 
##   0.001   0.000   8.007
```

```r
library(parallel)
detectCores()
```

```
## [1] 4
```

```r
options(mc.cores=4)
system.time( mclapply(1:8, function(x)Sys.sleep(1) ) )
```

```
##    user  system elapsed 
##   0.011   0.027   2.024
```

For this meaningless computation, we achieved linear speedup:
we cut the time for serial computation by a factor of four.

### Detour: BAM files from an RNA-seq experiment

We will be working with a set of BAM files
created in the study of [Zarnack et al. 2013](http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3629564&tool=pmcentrez&rendertype=abstract).
The hypothesis of this study is that a class of proteins,
the heterogeneous nuclear ribonucleoproteins C1 and C2,
referred to as hnRNP C, are responsible for preventing
incorporation of Alu elements into mature RNA transcripts.
Such incorporation may lead to protein dysfunction that could
cause disease; see the references of the Zarnack paper for
further background.

The package that we will work with has 8 BAM files with
reads aligned to chr14.  Four of the files are reads obtained from
from HeLa cells (negative control) and four are obtained from
HeLa cells in which hnRNP C has been "knocked down" with siRNA.

```r
library(RNAseqData.HNRNPC.bam.chr14)
dir(system.file("extdata", package="RNAseqData.HNRNPC.bam.chr14"))
```

```
##  [1] "ERR127302_chr14.bam"     "ERR127302_chr14.bam.bai"
##  [3] "ERR127303_chr14.bam"     "ERR127303_chr14.bam.bai"
##  [5] "ERR127304_chr14.bam"     "ERR127304_chr14.bam.bai"
##  [7] "ERR127305_chr14.bam"     "ERR127305_chr14.bam.bai"
##  [9] "ERR127306_chr14.bam"     "ERR127306_chr14.bam.bai"
## [11] "ERR127307_chr14.bam"     "ERR127307_chr14.bam.bai"
## [13] "ERR127308_chr14.bam"     "ERR127308_chr14.bam.bai"
## [15] "ERR127309_chr14.bam"     "ERR127309_chr14.bam.bai"
## [17] "tophat2_out"
```
These are Tabix-indexed BAM files.  In the alphabetic ordering, the
first four are two pairs of replicates of HeLa cells
that have undergone hnRNP C knockdown, and the second four are
two pairs of control replicates.


### Implicit parallelism through BiocParallel

To foster harmonious development of reliably performant procedures,
Bioconductor introduced the BiocParallel package.  Users
benefit from autonomous (but optionally controllable) pursuit
of parallel computing when feasible.  Consider the following
example: we will count reads into bins defined by exon addresses
in the HNRNPC example dataset.  
The RNAseqData.HNRNPC.bam.chr14 package includes a vector
of absolute pathnames of the BAM files, which we assign to `fns`.


```r
library(RNAseqData.HNRNPC.bam.chr14)
fns = RNAseqData.HNRNPC.bam.chr14_BAMFILES
length(fns)
```

```
## [1] 8
```

Here we establish the exon bins into which we will count reads.

```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) = "chr14"
ebg = exonsBy(txdb, by="gene")
```

Now we use the `summarizeOverlaps` function from
the GenomicAlignments package to count reads into the exon bins.
We'll time the counting from a single file, and then time
the counting from all files at once.

```r
library(GenomicAlignments)
# summarizeOverlaps uses bplapply to iterate over files
s1 = system.time(i1 <- summarizeOverlaps( ebg, fns[3] ))
s1
```

```
##    user  system elapsed 
##   7.565   0.266   7.872
```

```r
# show implicit config
BiocParallel::bpparam()
```

```
## class: MulticoreParam; bpisup: TRUE; bpworkers: 4; catch.errors: TRUE
## setSeed: TRUE; recursive: TRUE; cleanup: TRUE; cleanupSignal: 15;
##   verbose: FALSE
```

```r
s2 = system.time(i2 <- summarizeOverlaps( ebg, fns ))
s2
```

```
##    user  system elapsed 
##  36.426   3.659  23.955
```
This is not a thorough way of measuring speedup but it
shows reasonable enhancement.  
In the second computation, we did approximately 8 times as
much computation, but the clock time elapsed increased only
by a factor of (3.04).  
We did nothing by way of configuration
or request.

What happened?  The `summarizeOverlaps` function will iterate
over the files using `bplapply` from the BiocParallel package.
That function will check the R session 
for specific parallelization configuration information,
and if it finds none, will check for multiple cores
and make arrangements to use them if present.
The "check" occurs via the function `bpparam`.

The default situation on a
MacBook Air running MacOSX 10.9.5 with
an Intel core i7 processor,
(two physical cores with two logical cores each, allowing
for four concurrent threads)
as follows.

```r
library(BiocParallel)
bpparam()
```

```
## class: MulticoreParam; bpisup: TRUE; bpworkers: 4; catch.errors: TRUE
## setSeed: TRUE; recursive: TRUE; cleanup: TRUE; cleanupSignal: 15;
##   verbose: FALSE
```
This identifies an object called `MulticoreParam` which is
used to configure the behavior of `bplapply` and other utilities
of the BiocParallel package.  There are various configuration
object classes that can be used.

```
‘SnowParam’: distributed memory computing

‘MulticoreParam’: shared memory computing

‘BatchJobsParam’: scheduled cluster computing

‘DoparParam’: foreach computing

‘SerialParam’: non-parallel execution
```
We need to use `register` to determine the type of
concurrent computation that will be performed.  

If process size is large, we may want to leave
some cores idle.  We can accomplish that by using `register`.

```r
library(BiocParallel)
register(MulticoreParam(workers=2))
system.time(i3 <- summarizeOverlaps( ebg, fns ))
```

```
##    user  system elapsed 
##  34.044   2.690  23.263
```

```r
all.equal(i3,i2)  # check that the results do not change
```

```
## [1] TRUE
```
Note here that by reducing the number of CPUs by a factor of 2, we
do not double run time.  This is because there are communication
costs that are reduced by reducing the number of CPUs.

In summary, it is very easy to perform embarrassignly parallel
tasks with R, and this carries over to genomic data analysis
thanks to BiocParallel.  There are some strategic considerations
concerning control of memory consumption and communication costs,
and full mastery of the topic involves attention to profiling and
benchmarking methods that can be addressed in an advanced software
development course.


