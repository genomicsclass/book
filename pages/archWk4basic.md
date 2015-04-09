---
layout: page
title: "Architecture: Overview of last of the four As"
---



## Introduction to architectural concepts for Bioconductor

The basic objective is to support an efficient and reliable flow of experimental
and reference data.

Start with

- Assay outputs bound to sample-level data

Pass to

- Algorithms for preprocessing to remove technical artifacts

Combine clean assay outputs with

- Annotation on genome structure and function and on experimental design

Continue with

- Algorithms for inference on biological hypotheses

Conclude with

- Efficient and appropriate reporting, visualization and export

As noted previously, your experiments and analyses may serve
as data and annotation for future experiments in other labs.

In this subunit we want to clarify some of the architectural
principles underlying Bioconductor so that this objective of
efficient and reliable data flow can be
achieved *making good use of community collaboration and
total commitment to open source*.

Our key topics will be:

- How to create R packages

- How to support integrative analysis of multiple assay types

- How to streamline access to curated institutional archives like GEO

- How to make good use of parallel computing concepts on laptops and clusters

- How Bioconductor ensures reliable interoperability of project software and data assets


## What is an R package?

Conceptually, an R package is a collection of functions, data
objects, and documentation that coherently support a family
of related data analysis operations.

Concretely, an R package is a structured collection of folders,
organized and populated according to the rules of
[Writing R Extensions](http://cran.r-project.org/doc/manuals/r-release/R-exts.html).

### A new software package with `package.skeleton`

We can create our own packages using `package.skeleton`.  We'll illustrate that now
with an enhancement to the ERBS package that was created for the course.
We'll create a new package that utilizes the peak data, defining
a function `juxta` that allows us to compare binding peak patterns for the two cell
types on a selected chromosome.

Here's a definition of `juxta`.  Add it to your R session.

```r
juxta = function (chrname="chr22", ...) 
{
    require(ERBS)
    data(HepG2)
    data(GM12878)
    require(ggbio)
    require(GenomicRanges)  # "subset" is overused, need import detail
    ap1 = autoplot(GenomicRanges::subset(HepG2, seqnames==chrname))
    ap2 = autoplot(GenomicRanges::subset(GM12878, seqnames==chrname))
    tracks(HepG2 = ap1, Bcell = ap2, ...)
}
```

Now demonstrate it as follows.


```r
library(ERBS)
juxta("chr22", main="ESRRA binding peaks on chr22")
```

```
## Loading required package: ggbio
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
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: ggplot2
## Need specific help about ggbio? try mailing 
##  the maintainer or visit http://tengfei.github.com/ggbio/
## 
## Attaching package: 'ggbio'
## 
## The following objects are masked from 'package:ggplot2':
## 
##     geom_bar, geom_rect, geom_segment, ggsave, stat_bin,
##     stat_identity, xlim
## 
## Loading required package: GenomicRanges
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
```

![plot of chunk doj](figure/archWk4basic-doj-1.png) 

In the video we will show how to use `package.skeleton` and the Rstudio
editor to generate, document, and install this new package!  We will
streamline the code in `juxta` to make use of inter-package
symbol transfer by properly writing the DESCRIPTION and NAMESPACE
files for the package.

### A new annotation package with OrganismDbi

We have found the `Homo.sapiens` package to be quite convenient.
We can get gene models, symbol to GO mappings, and so on, without
remembering any more than `keytypes`, `columns`, `keys`, and `select`.
At present there is no similar resource for *S. cerevisiae*.
We can make one, following the OrganismDbi vignette.  This is
a very lightweight integrative package.


```r
library(OrganismDbi)
```

```
## Loading required package: AnnotationDbi
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## 
## Attaching package: 'AnnotationDbi'
## 
## The following object is masked from 'package:GenomeInfoDb':
## 
##     species
## 
## Loading required package: GenomicFeatures
```

```r
gd = list( join1 = c(GO.db="GOID", org.Sc.sgd.db="GO"),
           join2 = c(org.Sc.sgd.db="ENTREZID",
              TxDb.Scerevisiae.UCSC.sacCer3.sgdGene="GENEID"))
if (!file.exists("Sac.cer3")) # don't do twice...
makeOrganismPackage(pkgname="Sac.cer3",  # simplify typing!
  graphData=gd, organism="Saccharomyces cerevisiae",
  version="1.0.0", maintainer="Student <ph525x@harvardx.edu>",
  author="Student <ph525x@harvardx.edu>",
  destDir=".",
  license="Artistic-2.0")
```

At this point we have a folder structure in our
working folder that can support an installation.

```r
install.packages("Sac.cer3", repos=NULL, type="source")
library(Sac.cer3)
```

```
## Loading required package: GO.db
## Loading required package: DBI
## 
## Loading required package: org.Sc.sgd.db
## 
## Loading required package: TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
```

```r
Sac.cer3
```

```
## class: OrganismDb 
## Annotation resources:
## [1] "GO.db"                                
## [2] "org.Sc.sgd.db"                        
## [3] "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"
## Annotation relationships:
##      xDbs            yDbs                                    xKeys     
## [1,] "GO.db"         "org.Sc.sgd.db"                         "GOID"    
## [2,] "org.Sc.sgd.db" "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" "ENTREZID"
##      yKeys   
## [1,] "GO"    
## [2,] "GENEID"
## For more details, please see the show methods for the component objects listed above.
```

```r
columns(Sac.cer3)
```

```
##  [1] "GOID"         "TERM"         "ONTOLOGY"     "DEFINITION"  
##  [5] "ORF"          "ENTREZID"     "PFAM"         "ENZYME"      
##  [9] "PATH"         "PMID"         "REFSEQ"       "ENSEMBL"     
## [13] "ENSEMBLPROT"  "ENSEMBLTRANS" "UNIPROT"      "GO"          
## [17] "EVIDENCE"     "GOALL"        "EVIDENCEALL"  "ONTOLOGYALL" 
## [21] "DESCRIPTION"  "COMMON"       "INTERPRO"     "SMART"       
## [25] "SGD"          "ALIAS"        "CHRLOC"       "CHRLOCEND"   
## [29] "GENENAME"     "CHR"          "CDSID"        "CDSNAME"     
## [33] "CDSCHROM"     "CDSSTRAND"    "CDSSTART"     "CDSEND"      
## [37] "EXONID"       "EXONNAME"     "EXONCHROM"    "EXONSTRAND"  
## [41] "EXONSTART"    "EXONEND"      "GENEID"       "TXID"        
## [45] "EXONRANK"     "TXNAME"       "TXCHROM"      "TXSTRAND"    
## [49] "TXSTART"      "TXEND"
```

```r
genes(Sac.cer3)
```

```
## GRanges object with 6534 ranges and 1 metadata column:
##             seqnames           ranges strand   |        GENEID
##                <Rle>        <IRanges>  <Rle>   | <IntegerList>
##       Q0010     chrM   [ 3952,  4415]      +   |         Q0010
##       Q0032     chrM   [11667, 11957]      +   |         Q0032
##       Q0055     chrM   [13818, 26701]      +   |         Q0055
##       Q0075     chrM   [24156, 25255]      +   |         Q0075
##       Q0080     chrM   [27666, 27812]      +   |         Q0080
##         ...      ...              ...    ... ...           ...
##     YPR200C   chrXVI [939279, 939671]      -   |       YPR200C
##     YPR201W   chrXVI [939922, 941136]      +   |       YPR201W
##     YPR202W   chrXVI [943032, 944188]      +   |       YPR202W
##   YPR204C-A   chrXVI [946856, 947338]      -   |     YPR204C-A
##     YPR204W   chrXVI [944603, 947701]      +   |       YPR204W
##   -------
##   seqinfo: 17 sequences (1 circular) from sacCer3 genome
```

## Integrative analysis concepts
 
### TF binding and expression co-regulation

An example of integrative analysis was given in the introductory
lecture, in connection with the regulatory program of the yeast 
cell cycle.  There are two key experimental components:

- Protein binding patterns: based on ChIP-chip experiments, we can determine
the gene promoter regions to which transcription factors bind.

- Expression patterns: based on timed observations of gene expression in a yeast colony
we can identify times at which groups of genes reach maximal expression.


The diagram that we looked at indicated that the Mbp1 transcription
factor played a role in regulating expression in the transition
from G1 to S phases of the cell cycle.  The ChIP-chip data is
in the `harbChIP` package.


```r
library(harbChIP)
```

```
## Loading required package: tools
## Loading required package: Biostrings
## 
## Attaching package: 'harbChIP'
## 
## The following object is masked from 'package:OrganismDbi':
## 
##     keys
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     keys
## 
## The following object is masked from 'package:GenomeInfoDb':
## 
##     organism
```

```r
data(harbChIP)
harbChIP
```

```
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 6230 features, 204 samples 
##   element names: exprs, se.exprs 
## protocolData: none
## phenoData
##   sampleNames: A1 (MATA1) ABF1 ... ZMS1 (204 total)
##   varLabels: txFac
##   varMetadata: labelDescription
## featureData
##   featureNames: YAL001C YAL002W ... MRH1 (6230 total)
##   fvarLabels: ID PLATE ... REV_SEQ (12 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
##   pubMedIds: 15343339 
## Annotation:
```
This is a well-documented data object, and we can read the abstract
of the paper directly.


```r
abstract(harbChIP)
```

```
## [1] "DNA-binding transcriptional regulators interpret the genome's regulatory code by binding to specific sequences to induce or repress gene expression. Comparative genomics has recently been used to identify potential cis-regulatory sequences within the yeast genome on the basis of phylogenetic conservation, but this information alone does not reveal if or when transcriptional regulators occupy these binding sites. We have constructed an initial map of yeast's transcriptional regulatory code by identifying the sequence elements that are bound by regulators under various conditions and that are conserved among Saccharomyces species. The organization of regulatory elements in promoters and the environment-dependent use of these elements by regulators are discussed. We find that environment-specific use of regulatory elements predicts mechanistic models for the function of a large population of yeast's transcriptional regulators."
```

Let's find MBP1 and assess the distribution of reported binding affinity
measures.  The sample names of the ExpressionSet (structure
used for convenience
even though the data are not expression data)
are the names of the proteins "chipped" onto the yeast
promoter array.


```r
mind = which(sampleNames(harbChIP)=="MBP1")
qqnorm(exprs(harbChIP)[,mind], main="MBP1 binding")
```

![plot of chunk lkm](figure/archWk4basic-lkm-1.png) 

The shape of the qq-normal plot is indicative of
a strong
departure from Gaussianity in the distribution
of binding scores, with a very long right tail.
We'll focus on the top five genes.


```r
topb = featureNames(harbChIP)[ order(
  exprs(harbChIP)[,mind], decreasing=TRUE)[1:5] ]
topb
```

```
## [1] "YDL101C" "YPR075C" "YCR064C" "YCR065W" "YGR109C"
```

```r
library(org.Sc.sgd.db)
select(org.Sc.sgd.db, keys=topb, keytype="ORF",
  columns="COMMON")
```

```
##       ORF        SGD COMMON
## 1 YDL101C S000002259   DUN1
## 2 YPR075C S000006279   OPY2
## 3 YCR064C S000000660   <NA>
## 4 YCR065W S000000661   HCM1
## 5 YGR109C S000003341   CLB6
```

Our conjecture is that these genes will exhibit
similar expression trajectories, peaking well
within the first half of the 66 minute cell cycle
for the yeast strain studied.

We will subset the cell cycle expression data from
the `yeastCC` package to a colony whose cycling was
synchronized using alpha pheromone.


```r
library(yeastCC)
data(spYCCES)
alp = spYCCES[, spYCCES$syncmeth=="alpha"]
par(mfrow=c(1,1))
plot(exprs(alp)[ topb[1], ]~alp$time, lty=1,
   type="l", ylim=c(-1.5,1.5), lwd=2, ylab="Expression",
    xlab="Minutes elapsed")
for (i in 2:5) lines(exprs(alp)[topb[i],]~alp$time, lty=i, lwd=2)
legend(75,-.5, lty=1:10, legend=topb, lwd=2, cex=.6, seg.len=4)
```

![plot of chunk doalp](figure/archWk4basic-doalp-1.png) 

We have the impression that at least three of these
genes reach peak expression roughly together near times
20 and 80 minutes.  There is considerable variability.
A data filtering and visualization pattern is emerging
by which genes bound by a given transcription factor
can be assessed for coregulation of expression.  We
have not entered into the assessment of statistical
significance, but have focused on how the data
types are brought together.

Consider how the Bioconductor architecture facilitated
this analysis.

- Easily installed and highly self-descriptive packages provide the key experimental data

- Conventional containers (ExpressionSets) are used
for assay plus sample-level data (even when the experiment
does not assess expression) so that it is easy
to quickly isolate features and samples of interest

- Immediate access to R's visualization and statistical
analysis functions makes appraisal and inference
very convenient.

### Harvesting GEO for families of microarray archives

The NCBI Gene Expression Omnibus is a basic resource for
integrative bioinformatics.  The Bioconductor GEOmetadb
package helps with discovery and retrieval of tools for
analyzing GEO datasets.

The GEOmetadb database is a 240MB download that decompresses to 3.6 GB
of SQLite.  Once you have acquired the GEOmetadb.sqlite file using
the `getSQLiteFile` function, you can create a connection
and start interrogating.


```r
library(RSQLite)
lcon = dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListTables(lcon)
```

```
##  [1] "gds"               "gds_subset"        "geoConvert"       
##  [4] "geodb_column_desc" "gpl"               "gse"              
##  [7] "gse_gpl"           "gse_gsm"           "gsm"              
## [10] "metaInfo"          "sMatrix"
```

We will build a query that returns all the GEO GSE entries
that have the phrase "pancreatic cancer" in their titles.
Because GEO uses uninformative labels for array platforms,
we will retrieve a field that records the Bioconductor array
annotation package name so that we know what technology was
in use.  We'll tabulate the various platforms used.


```r
vbls = "gse.gse, gse.title, gpl.gpl, gpl.bioc_package"
req1 = " from gse join gse_gpl on gse.gse=gse_gpl.gse"
req2 = " join gpl on gse_gpl.gpl=gpl.gpl"
goal = " where gse.title like '%pancreatic%cancer%'"
quer = paste0("select ", vbls, req1, req2, goal)
lkpc = dbGetQuery(lcon, quer)
dim(lkpc)
```

```
## [1] 137   4
```

```r
table(lkpc$bioc_package)
```

```
## 
##                     hgu133a                    hgu133a2 
##                           3                           3 
##                     hgu133b                 hgu133plus2 
##                           3                          27 
##                   hgug4110b       HsAgilentDesign026652 
##                           1                           1 
## hugene10sttranscriptcluster IlluminaHumanMethylation27k 
##                           1                           2 
##             illuminaHumanv4                  mouse430a2 
##                           7                           7
```

