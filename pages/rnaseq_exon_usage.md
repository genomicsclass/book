---
layout: page
title: RNA-seq differential exon usage
---



The [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) package offers differential testing of exon usage within each gene. Here we will explore the R code used in a *DEXSeq* analysis. We omit the python calls for preparing the annotation and count tables, but these can be found in the vignette at the above link. The python calls are generally along the lines of:

```
python dexseq_prepare_annotation.py gtffile.gtf dexseq.gff
python dexseq_count.py dexseq.gff sample1.sam sample1.txt
```

Once we have repeated the `dexseq_count` script for each sample, we can read the data into R using the code chunks below. As we are working with pre-prepared data, we first point to these files which live within the *pasilla* package. 

The *pasilla* package contains counts from an experiment by [Brooks et al](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3032923/)

We will run DEXSeq on a subset of the genes, for demonstration purposes.


```r
library("pasilla")
inDir = system.file("extdata", package="pasilla", mustWork=TRUE)
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
genesForSubset = read.table(file.path(inDir, "geneIDsinsubset.txt"),
  stringsAsFactors=FALSE)[[1]]
```

As in *DESeq2* we use a `sampleTable` to define the samples:


```r
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3",
    "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",
    "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end",
    "single-end", "single-end", "paired-end", "paired-end" ) )
sampleTable
```

```
##            condition    libType
## treated1   knockdown single-end
## treated2   knockdown paired-end
## treated3   knockdown paired-end
## untreated1   control single-end
## untreated2   control single-end
## untreated3   control paired-end
## untreated4   control paired-end
```

We now read the data into a `DEXSeqDataSet` object:


```r
library("DEXSeq")
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
```

Subset the genes, for demonstration purposes:


```r
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]
```

Now we run the estimation and testing functions:


```r
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
```

The following code extracts a results table, makes an MA-plot, and draws the expression levels over the exons to highlight differential exon usage:


```r
dxr = DEXSeqResults( dxd )
plotMA( dxr, cex=0.8 )
```

![plot of chunk unnamed-chunk-6](figure/rnaseq_exon_usage-unnamed-chunk-6-1.png) 

```r
plotDEXSeq( dxr, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
```

![plot of chunk unnamed-chunk-6](figure/rnaseq_exon_usage-unnamed-chunk-6-2.png) 

Again, drawing the expression levels, now showing the annotated transcripts below:


```r
plotDEXSeq( dxr, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE,
              cex.axis=1.2, cex=1.3, lwd=2 )
```

![plot of chunk unnamed-chunk-7](figure/rnaseq_exon_usage-unnamed-chunk-7-1.png) 

For more details on the *DEXSeq* software, see the vignette and the paper, which is linked from the vignette page:


```r
browseVignettes("DEXSeq")
```
