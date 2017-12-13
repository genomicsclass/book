---
layout: page
title: "Combining S4 with NoSQL (mongodb) to query ENCODE bedfiles"
---




## Introduction

For an application involving many thousands of files,
we have found that [NoSQL](https://en.wikipedia.org/wiki/NoSQL)
strategies may be effective.  The [TxRegInfra](https://github.com/vjcitn/TxRegInfra) package is under development in github 
and illustrates use of mongodb
with a small collection of BED files obtained from the
ENCODE project.  In this section we sketch the most
basic aspects of wrapping a mongodb connection in S4 and
implementing subsetByOverlaps to query the data store.

To carry out the tasks in this section, you will need mongod
(the database managing daemon) running on your system.  
The [community server edition](https://www.mongodb.com/download-center?jmp=nav#community) should be easy to install.

After you get mongod running, you can install TxRegInfra using
`library(BiocInstaller); biocLite("vjcitn/TxRegInfra")`.

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

### S4: Extending the RaggedExperiment class

(From the RaggedExperiment vignette:)
The *[RaggedExperiment](http://bioconductor.org/packages/RaggedExperiment)* package provides a flexible data
representation for copy number, mutation and other ragged array schema for 
genomic location data. It aims to provide a framework for a set of samples
that have differing numbers of genomic ranges.

In TxRegInfra, we extend the RaggedExperiment class to deal
with external data managed by mongodb.  We've created
a database 'txregnet' and we connect this to 
the extended RaggedExperiment
'rme1', an instance of `RaggedMongoExpt`.


```r
okdf = DataFrame(hsFiles)
rownames(okdf) = hsFiles[,1]
loccon = localMongolite(db="txregnet")
rme1 = RaggedMongoExpt(loccon, colData=okdf)
rme1
```

```
## class: RaggedMongoExpt 
## dim: 10 10 
## assays(0):
## rownames(10): ENCFF001UUW ENCFF936ICC ... ENCFF971VCD ENCFF350ZQV
## colnames(10): ENCFF001UUW ENCFF936ICC ... ENCFF971VCD ENCFF350ZQV
## colData names(10): File.accession File.format ...
##   Biosample.life.stage Biosample.sex
```

## The upshot: peak densities by tissue type

In the following, we produce a table of number of
peaks by tissue type, in a small region of chromosome 1.


```r
brp = which(colData(rme1)$File.format == "bed broadPeak")
allst = subsetByOverlaps(rme1[,brp], 
               GRanges("chr1", IRanges(1,8e5))) 
data.frame(tiss=colData(rme1)[brp, "Biosample.term.name"], 
             num.peaks=sapply(allst,nrow))
```

```
##                                    tiss num.peaks
## ENCFF936ICC      CD14-positive monocyte         9
## ENCFF879FWO                         eye        19
## ENCFF970PRR   renal cortex interstitium         7
## ENCFF001WQH stromal cell of bone marrow        18
## ENCFF942PVS                     GM12864         5
```

## Some additional details

Ultimately we would like to make use of the RaggedExperiment
infrastructure directly.  To do this we need to bind
a GRangesList to the assay data; once this is done,
we can use the sparseAssay, compactAssay, and qreduceAssay
methods.  Longer term utility of this approach will be
demonstrated in the TxRegQuery package, under development.


```r
badn = c("seqnames", "ranges", "strand", "seqlevels", 
   "seqlengths", "isCircular", "start", "end", "width", "element")
cleanCols = function(x) setdiff(colnames(x), badn)
grl = GRangesList(lapply(allst, function(x) {
     ans = GRanges(x$chrom, IRanges(x$chromStart, x$chromEnd)); mcols(ans) = x[,cleanCols(x)]; ans
     }))
re = RaggedExperiment(grl, colData=colData(rme1[,brp])) 
dim(sparseAssay(re))
```

```
## [1] 58  5
```

```r
dim(compactAssay(re))
```

```
## [1] 51  5
```

To conclude, we peek at the details of the mongodb connection
established by *[mongolite](https://CRAN.R-project.org/package=mongolite)*.  It includes a
variety of hints concerning the R interface.


```r
rme1@con
```

```
## mongolite connection for db txregnet, coll. test
## URL:  mongodb://127.0.0.1
```

```r
rme1@con@con
```

```
## <Mongo collection> 'test' 
##  $$aggregate(pipeline = "{}", options = "{\"allowDiskUse\":true}", handler = NULL, pagesize = 1000) 
##  $$count(query = "{}") 
##  $$distinct(key, query = "{}") 
##  $$drop() 
##  $$export(con = stdout(), bson = FALSE) 
##  $$find(query = "{}", fields = "{\"_id\":0}", sort = "{}", skip = 0, limit = 0, handler = NULL, pagesize = 1000) 
##  $$import(con, bson = FALSE) 
##  $$index(add = NULL, remove = NULL) 
##  $$info() 
##  $$insert(data, pagesize = 1000, ...) 
##  $$iterate(query = "{}", fields = "{\"_id\":0}", sort = "{}", skip = 0, limit = 0) 
##  $$mapreduce(map, reduce, query = "{}", sort = "{}", limit = 0, out = NULL, scope = NULL) 
##  $$remove(query, just_one = FALSE) 
##  $$rename(name, db = NULL) 
##  $$update(query, update = "{\"$set\":{}}", upsert = FALSE, multiple = FALSE)
```
