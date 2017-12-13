---
layout: page
title: "Combining S4 with NoSQL (mongodb) to query ENCODE bedfiles"
---




## Introduction

For an application involving many thousands of files,
we have found that [NoSQL](https://en.wikipedia.org/wiki/NoSQL)
strategies may be effective.  The [TxRegInfra](https://github.com/vjcitn/TxRegInfra) is under development and illustrates use of mongodb
with a small collection of BED files obtained from the
ENCODE project.  In this section we sketch the most
basic aspects of wrapping a mongodb connection in S4 and
implementing subsetByOverlaps to query the data store.

To carry out the tasks in this section, you will need mongod
(the database managing daemon) running on your system.  

## Basic considerations



Our long term goal is
to define a package,
TxRegQuery, to exploration of transcriptional regulatory networks
by integrating data on eQTL, digital genomic footprinting (DGF), DnaseI
hypersensitivity binding data (DHS), and transcription
factor binding site (TFBS) data.  Owing to the volume of emerging tissue-specific
data, special data modalities are used.  In this document
we'll focus on DHS.

## Managing bed file content with mongodb

### Importing and querying documents

The package comes with a small number of bed files to demonstrate
import utilities.

```r
# ENCODE
f1 = dir(system.file("bedfiles", package="TxRegInfra"), full=TRUE, patt="ENCFF971VCD")
cat(readLines(f1, n=3), sep="\n")
```

```
## chr1	181415	181565	.	0	.	26	-1	-1	75
## chr1	629855	630005	.	0	.	45	-1	-1	75
## chr1	629895	630045	.	0	.	44.0	-1	-1	75
```

```r
# ChromHMM
f2 = dir(system.file("bedfiles", package="TxRegInfra"), full=TRUE, patt="E096_imp12")
cat(readLines(f2, n=3), sep="\n")
```

```
## chr17	0	14400	25_Quies
## chr17	14400	14600	19_DNase
## chr17	14600	28200	25_Quies
```

The function `importBedToMongo` uses system() to run
mongodb.
There is a `bedType` parameter that indicates what fields are available; it
defaults to `broadPeak`.

The following code imports a broadPeak and chromHMM document.
We deal with metadata about these documents below.
We assume a database called 'txregnet' has been established
for a running mongodb server.

```r
importBedToMongo(f1, "vjc1", db="txregnet")
```

```
## [1] TRUE
```

```r
importBedToMongo(f2, "vjc2", db="txregnet", bedType="chromHMM")
```

```
## [1] TRUE
```

Now that the documents are imported, we can query for
information in an interval specified by a GRanges instance.

```r
library(RMongo)
```

```
## Loading required package: rJava
```

```r
con = mongoDbConnect("txregnet") # defaults for local server
queryBedInMongo(con, "vjc1", GRanges("chr1", IRanges(1, 800000)), skip=0, limit=5)
```

```
##   score chromStart chromEnd strand name pValue                     X_id
## 1     0     181415   181565      .    .     -1 5a2188f454be7d9ced2569e3
## 2     0     633955   634105      .    .     -1 5a2188f454be7d9ced2569e4
## 3     0     633995   634145      .    .     -1 5a2188f454be7d9ced2569e5
## 4     0     778635   778785      .    .     -1 5a2188f454be7d9ced2569e6
## 5     0     779135   779285      .    .     -1 5a2188f454be7d9ced2569e7
##   peak qValue signalValue chrom
## 1   75     -1          26  chr1
## 2   75     -1          63  chr1
## 3   75     -1          62  chr1
## 4   75     -1          76  chr1
## 5   75     -1          11  chr1
```

```r
queryBedInMongo(con, "vjc2", GRanges("chr17", IRanges(1, 800000)), skip=0, limit=5)
```

```
##   chromStart chromEnd                     X_id      state chrom
## 1      28200    29600 5a2188f554be7d9ced26a57a  24_ReprPC chr17
## 2      29600    30800 5a2188f554be7d9ced26a57b 23_PromBiv chr17
## 3      30800    31600 5a2188f554be7d9ced26a57c   22_PromP chr17
## 4      14600    28200 5a2188f554be7d9ced26a57d   25_Quies chr17
## 5      31600    32000 5a2188f554be7d9ced26a57e    2_PromU chr17
```

## An integrative container

We need to bind the metadata and information about the mongodb.

### BED file metadata

The BED files are extracted from a few different places.  We have
metadata on 10 of them:

```r
data(hsFiles_subset) # holds hsFiles
hsFiles[1:3,1:6]
```

```
##      File.accession    File.format Output.type Experiment.accession
## 1916    ENCFF001UUW bed narrowPeak       peaks          ENCSR000EIX
## 859     ENCFF936ICC  bed broadPeak    hotspots          ENCSR000EPK
## 1462    ENCFF460FCG bed narrowPeak       peaks          ENCSR000EOI
##          Assay Biosample.term.id
## 1916 DNase-seq        CL:0002620
## 859  DNase-seq        CL:0001054
## 1462 DNase-seq        CL:2000017
```
We added an additional four.  This will become colData for an
instance of an extended RaggedExperiment class to be defined.


```r
cd = DataFrame(rbind(hsFiles, rbind(e072, e073, e088, e096)))
cd[1:4,1:6]
```

```
## DataFrame with 4 rows and 6 columns
##   File.accession    File.format Output.type Experiment.accession
##      <character>    <character> <character>          <character>
## 1    ENCFF001UUW bed narrowPeak       peaks          ENCSR000EIX
## 2    ENCFF936ICC  bed broadPeak    hotspots          ENCSR000EPK
## 3    ENCFF460FCG bed narrowPeak       peaks          ENCSR000EOI
## 4    ENCFF879FWO  bed broadPeak    hotspots          ENCSR782XFY
##         Assay Biosample.term.id
##   <character>       <character>
## 1   DNase-seq        CL:0002620
## 2   DNase-seq        CL:0001054
## 3   DNase-seq        CL:2000017
## 4   DNase-seq    UBERON:0000970
```








