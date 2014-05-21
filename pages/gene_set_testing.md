---
layout: page
title: Gene set testing
---




# Gene set testing

In the following unit, we will explore software for testing differential expression in a set of genes. These tests differ from the gene-by-gene tests we saw previously. Again, the gene set testing software we will use lives in the `limma` package.

We download an experiment from the GEO website, using the `getGEO` function from the `GEOquery` package:


```r
library(GEOquery)
```

```
## Loading required package: methods
## Loading required package: Biobase
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
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Setting options('download.file.method.GEOquery'='auto')
```

```r
g <- getGEO("GSE34313")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34313/matrix/
## Found 1 file(s)
## GSE34313_series_matrix.txt.gz
## File stored at: 
## /var/folders/6d/d_8pbllx7318htlp5wv_rm580000gn/T//RtmpY5CFbp/GPL6480.soft
```

```r
e <- g[[1]]
```


This dataset is hosted by GEO at the following link: <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34313>

The experiment is described in the paper by [Masuno 2011](#foot).

Briefly, the investigators applied a glucocorticoid hormone to cultured human airway smooth muscle. The glucocorticoid hormone is used to treat asthma, as it reduces the inflammation response, however it has many other effects throughout the different tissues of the body.

The groups are defined in the `characteristics_ch1.2` variable:


```r
e$condition <- e$characteristics_ch1.2
levels(e$condition) <- c("dex24", "dex4", "control")
table(e$condition)
```

```
## 
##   dex24    dex4 control 
##       3       3       4
```


By examining boxplots, we can guess that the data has already been normalized somehow, and on the GEO site the investigators report that they normalized using Agilent software.

We will subset to the control samples and the samples treated with dexamethasone (the hormone) after 4 hours.


```r
boxplot(exprs(e), range = 0)
```

![plot of chunk unnamed-chunk-3](figure/gene_set_testing-unnamed-chunk-3.png) 

```r
names(fData(e))
```

```
##  [1] "ID"                   "SPOT_ID"              "CONTROL_TYPE"        
##  [4] "REFSEQ"               "GB_ACC"               "GENE"                
##  [7] "GENE_SYMBOL"          "GENE_NAME"            "UNIGENE_ID"          
## [10] "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
## [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"         
## [16] "GO_ID"                "SEQUENCE"
```

```r
lvls <- c("control", "dex4")
es <- e[, e$condition %in% lvls]
es$condition <- factor(es$condition, levels = lvls)
```


The following lines run the linear model in `limma`. We note that the top genes are common immune-response genes (CSF2, LIF, CCL2, IL6). Also present is FKBP5, a gene which regulates and is regulated by the protein which receives the glucocorticoid hormone.


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
design <- model.matrix(~es$condition)
fit <- lmFit(es, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, genelist = fData(es)$GENE_SYMBOL)
tt
```

```
##                    ID  logFC AveExpr      t   P.Value adj.P.Val     B
## A_23_P133408     CSF2 -4.725   8.188 -69.36 2.382e-12 6.249e-08 15.93
## A_24_P122137      LIF -6.690  11.053 -67.24 3.048e-12 6.249e-08 15.85
## A_32_P5276   ARHGEF26  3.411   7.374  54.48 1.619e-11 1.745e-07 15.17
## A_23_P89431      CCL2 -3.619  13.763 -54.14 1.702e-11 1.745e-07 15.15
## A_23_P42257      IER3 -4.809  13.417 -48.58 4.018e-11 3.295e-07 14.72
## A_23_P71037       IL6 -4.568  11.595 -45.79 6.422e-11 4.213e-07 14.47
## A_24_P20327     KLF15  4.013   6.825  43.91 8.952e-11 4.213e-07 14.28
## A_24_P38081     FKBP5  3.567   9.310  43.77 9.176e-11 4.213e-07 14.26
## A_23_P213944    HBEGF -2.902  10.269 -43.73 9.247e-11 4.213e-07 14.26
## A_24_P250922    PTGS2 -3.518  10.477 -42.53 1.153e-10 4.728e-07 14.13
```



We will use the [ROAST method](#foot) for gene set testing. We can test a single gene set by looking up the genes which contain a certain GO ID, and providing this to the `roast` function. We will show how to get such lists of genes associated with a GO ID in the next chunk.

The roast function performs an advanced statistical technique, *rotation of residuals*, in order to generate a sense of the null distribution for the test statistic. The test statistics in this case is the summary of the scores from each gene. The tests are *self-contained* because only the summary for a single set is used, whereas other gene set tests might compare a set to all the other genes in the dataset, e.g., a *competitive* gene set test.

The result here tells us that the *immune response* genes are significantly down-regulated, and additionally, mixed up and down.


```r
# Immune response
idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)
```

```
## [1] 504
```

```r
r1 <- roast(es, idx, design)
# ?roast
r1
```

```
##       Active.Prop P.Value
## Down      0.16270   0.016
## Up        0.09325   0.985
## Mixed     0.25595   0.009
```



## Testing multiple gene sets

We can also use the `mroast` function to perform multiple roast tests. First we need to create a list, which contains the indices of genes in the ExpressionSet for each of a number of gene sets. We will use the `org.Hs.eg.db` package to gather the gene set information.


```r
# biocLite('org.Hs.eg.db')
library(org.Hs.eg.db)
```

```
## Loading required package: AnnotationDbi
## Loading required package: GenomeInfoDb
## Loading required package: DBI
```

```r
org.Hs.egGO2EG
```

```
## GO2EG map for Human (object of class "Go3AnnDbBimap")
```

```r
go2eg <- as.list(org.Hs.egGO2EG)
head(go2eg)
```

```
## $`GO:0000002`
##     TAS     TAS     ISS     IMP     NAS     IMP     IEA     IMP 
##   "291"  "1890"  "4205"  "4358"  "9361" "10000" "80119" "92667" 
## 
## $`GO:0000003`
##    IEP 
## "8510" 
## 
## $`GO:0000012`
##         IDA         IDA         IEA         IMP         IDA         IDA 
##      "3981"      "7141"      "7515"     "23411"     "54840"     "55775" 
##         IMP         IMP         IEA 
##     "55775"    "200558" "100133315" 
## 
## $`GO:0000018`
##     TAS     TAS     TAS     IMP     IMP     IEP 
##  "3575"  "3836"  "3838"  "9984" "10189" "56916" 
## 
## $`GO:0000019`
##     TAS     IDA 
##  "4361" "10111" 
## 
## $`GO:0000022`
##    TAS    TAS 
## "9055" "9493"
```


The following code unlists the list, then gets matches for each Entrez gene ID to the index in the ExpressionSet. Finally, we rebuild the list.


```r
govector <- unlist(go2eg)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)
```

```
## [1] "400451" "10239"  "9899"   "348093" "57099"  "146050"
```

```r
idxvector <- match(govector, fData(es)$GENE)
table(is.na(idxvector))
```

```
## 
##  FALSE   TRUE 
## 206613   6488
```

```r
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
```

```
##     TAS     TAS     ISS     IMP     NAS     IMP     IEA     IMP 
##   "291"  "1890"  "4205"  "4358"  "9361" "10000" "80119" "92667"
```

```r
fData(es)$GENE[idx[[1]]]
```

```
## [1] "291"   "1890"  "4205"  "4358"  "9361"  "10000" "80119" "92667"
```


We need to clean this list such that there are no `NA` values. We also clean it to remove gene sets which have less than 10 genes.


```r
idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)
```

```
## [1] 2739
```


The following line of code runs the multiple ROAST test. This can take about 3 minutes.


```r
r2 <- mroast(es, idxsub, design)
head(r2)
```

```
##            NGenes PropDown  PropUp Direction PValue     FDR PValue.Mixed
## GO:0006954    297   0.2290 0.09428      Down  0.002 0.01497        0.001
## GO:0008083    167   0.2874 0.07784      Down  0.002 0.01497        0.001
## GO:0005125    166   0.2651 0.04819      Down  0.002 0.01497        0.001
## GO:0043433     61   0.2787 0.11475      Down  0.002 0.01497        0.001
## GO:0006959     52   0.2500 0.09615      Down  0.002 0.01497        0.001
## GO:0042102     48   0.1875 0.06250      Down  0.002 0.01497        0.001
##            FDR.Mixed
## GO:0006954  0.002365
## GO:0008083  0.002365
## GO:0005125  0.002365
## GO:0043433  0.002365
## GO:0006959  0.002365
## GO:0042102  0.002365
```

```r
r2 <- r2[order(r2$PValue.Mixed), ]
```


We can use the `GO.db` annotation package to extract the GO terms for the top results, by the *mixed* test.


```r
# biocLite('GO.db')
library(GO.db)
```

```
## 
```

```r
columns(GO.db)
```

```
## [1] "GOID"       "TERM"       "ONTOLOGY"   "DEFINITION"
```

```r
keytypes(GO.db)
```

```
## [1] "GOID"       "TERM"       "ONTOLOGY"   "DEFINITION"
```

```r
GOTERM[[rownames(r2)[1]]]
```

```
## GOID: GO:0006954
## Term: inflammatory response
## Ontology: BP
## Definition: The immediate defensive reaction (by vertebrate
##     tissue) to infection or injury caused by chemical or physical
##     agents. The process is characterized by local vasodilation,
##     extravasation of plasma into intercellular spaces and
##     accumulation of white blood cells and macrophages.
```

```r
r2tab <- select(GO.db, keys = rownames(r2)[1:10], columns = c("GOID", "TERM", 
    "DEFINITION"), keytype = "GOID")
r2tab[, 1:2]
```

```
##          GOID
## 1  GO:0006954
## 2  GO:0008083
## 3  GO:0005125
## 4  GO:0043433
## 5  GO:0006959
## 6  GO:0042102
## 7  GO:0051781
## 8  GO:0007623
## 9  GO:0048661
## 10 GO:0043525
##                                                                                  TERM
## 1                                                               inflammatory response
## 2                                                              growth factor activity
## 3                                                                   cytokine activity
## 4  negative regulation of sequence-specific DNA binding transcription factor activity
## 5                                                             humoral immune response
## 6                                         positive regulation of T cell proliferation
## 7                                                positive regulation of cell division
## 8                                                                    circadian rhythm
## 9                             positive regulation of smooth muscle cell proliferation
## 10                                    positive regulation of neuron apoptotic process
```


We can also look for the top results using the standard p-value and in the *up* direction.


```r
r2 <- r2[order(r2$PValue), ]
r2tab <- select(GO.db, keys = rownames(r2)[r2$Direction == "Up"][1:10], columns = c("GOID", 
    "TERM", "DEFINITION"), keytype = "GOID")
r2tab[, 1:2]
```

```
##          GOID
## 1  GO:0042594
## 2  GO:0030032
## 3  GO:0003950
## 4  GO:0042813
## 5  GO:0046847
## 6  GO:0006471
## 7  GO:0005528
## 8  GO:0071889
## 9  GO:0008209
## 10 GO:0008631
##                                                                     TERM
## 1                                                 response to starvation
## 2                                                 lamellipodium assembly
## 3                                   NAD+ ADP-ribosyltransferase activity
## 4                                        Wnt-activated receptor activity
## 5                                                    filopodium assembly
## 6                                               protein ADP-ribosylation
## 7                                                          FK506 binding
## 8                                                 14-3-3 protein binding
## 9                                             androgen metabolic process
## 10 intrinsic apoptotic signaling pathway in response to oxidative stress
```


Again but for the *down* direction.


```r
r2tab <- select(GO.db, keys = rownames(r2)[r2$Direction == "Down"][1:5], columns = c("GOID", 
    "TERM", "DEFINITION"), keytype = "GOID")
r2tab[, 1:2]
```

```
##         GOID
## 1 GO:0006954
## 2 GO:0008083
## 3 GO:0005125
## 4 GO:0043433
## 5 GO:0006959
##                                                                                 TERM
## 1                                                              inflammatory response
## 2                                                             growth factor activity
## 3                                                                  cytokine activity
## 4 negative regulation of sequence-specific DNA binding transcription factor activity
## 5                                                            humoral immune response
```



## Footnotes <a name="foot"></a>

### Methods within the limma package

Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK. "ROAST: rotation gene set tests for complex microarray experiments". Bioinformatics. 2010.
<http://www.ncbi.nlm.nih.gov/pubmed/20610611>

Di Wu and Gordon K. Smyth, "Camera: a competitive gene set test accounting for inter-gene correlation" Nucleic Acids Research, 2012.
<http://nar.oxfordjournals.org/content/40/17/e133>

### GSEA

Subramanian A1, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP, "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles" Proc Natl Acad Sci U S A. 2005.
<http://www.ncbi.nlm.nih.gov/pubmed/16199517>

### Correlation within gene sets

William T. Barry, Andrew B. Nobel, and Fred A. Wright, "A statistical framework for testing functional categories in microarray data" Ann. Appl. Stat, 2008.
<http://projecteuclid.org/euclid.aoas/1206367822>

William Barry has a package `safe` in Bioconductor for gene set testing with resampling.
<http://www.bioconductor.org/packages/release/bioc/html/safe.html>

Daniel M Gatti, William T Barry, Andrew B Nobel, Ivan Rusyn and Fred A Wright, "Heading Down the Wrong Pathway: on the Influence of Correlation within Gene Sets", BMC Genomics, 2010.
<http://www.biomedcentral.com/1471-2164/11/574#B24>

### Gene sets and power

The following article points out an issue with gene set testing: the power to detect differential expression for an individual gene depends on the number of NGS reads which align to that gene, which depends on the transcript length among other factors.

Alicia Oshlack* and Matthew J Wakefield, "Transcript length bias in RNA-seq data confounds systems biology", Biology Direct, 2009.
<http://www.biologydirect.com/content/4/1/14>

### The dataset used in this lab

Masuno K, Haldar SM, Jeyaraj D, Mailloux CM, Huang X, Panettieri RA Jr, Jain MK, Gerber AN., "Expression profiling identifies Klf15 as a glucocorticoid target that regulates airway hyperresponsiveness". Am J Respir Cell Mol Biol. 2011.
<http://www.ncbi.nlm.nih.gov/pubmed/21257922>
