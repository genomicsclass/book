# Exploring Cufflinks output with cummeRbund

Here we show the exploratory plots offered by the [cummeRbund](http://www.bioconductor.org/packages/release/bioc/html/cummeRbund.html) package. These plots require loading in a directory in which results from a [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) analysis has been run. Follow the vignette in the above link in order in order to perform a Cufflinks gene- and isoform-level analysis. From the vignette:

> CummeRbund begins by re-organizing output files of a cuffdiff analysis, and storing these data in a local SQLite database. CummeRbund indexes the data to speed up access to specific feature data (genes, isoforms, TSS, CDS, etc.), and preserves the various relationships between these features. 


```r
library(cummeRbund)
myDir <- system.file("extdata", package="cummeRbund") 
gtfFile <- system.file("extdata/chr1_snippet.gtf",package="cummeRbund")
```

Read in the prepared Cufflinks files from the directory:


```r
cuff <- readCufflinks(dir=myDir,gtfFile=gtfFile,genome="hg19",rebuild=TRUE)
```

Boxplots of expression (FPKM) at the gene and isoform level:


```r
csBoxplot(genes(cuff))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```r
csBoxplot(genes(cuff),replicates=TRUE)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) 

```r
csBoxplot(isoforms(cuff),replicates=TRUE)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png) 

Scatterplot matrix of gene and isoform level expression:


```r
csScatterMatrix(genes(cuff))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
csScatterMatrix(isoforms(cuff))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) 

Sample dendrograms using Jensen-Shannon distances:


```r
csDendro(genes(cuff),replicates=TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

```
## 'dendrogram' with 2 branches and 6 members total, at height 0.2685017
```

```r
csDendro(isoforms(cuff),replicates=TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png) 

```
## 'dendrogram' with 2 branches and 6 members total, at height 0.4377249
```

MA-plot comparing two conditions:


```r
MAplot(genes(cuff),"hESC","Fibroblasts")
```

```
## Warning in loop_apply(n, do.ply): Removed 54 rows containing missing
## values (geom_point).
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

```r
MAplot(isoforms(cuff),"hESC","Fibroblasts")
```

```
## Warning in loop_apply(n, do.ply): Removed 187 rows containing missing
## values (geom_point).
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png) 

A "volcano plot" matrix. Each volcano plot is the -log10(p-value) over the log fold change.


```r
csVolcanoMatrix(genes(cuff))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

```r
csVolcanoMatrix(isoforms(cuff))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png) 

For all of these functions, see the help pages in the *cummeRbund* package for more details, and check the vignette for a sample workflow. The [Cufflinks homepage](http://cole-trapnell-lab.github.io/cufflinks/) has details about running the pipeline upstream of producing these figures.


```r
browseVignettes("cummeRbund")
```
