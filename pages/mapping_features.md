---
layout: page
title: Mapping features to genes
---




## Using Bioconductor annotation packages

This unit will focus on mapping features to genes, i.e., getting annotation information from one format to another. We start by loading in the `maPooling` dataset from previous lectures.


```r
# library(devtools) install_github('dagdata','genomicsclass')
library(dagdata)
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
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
data(maPooling)
e <- maPooling
head(rownames(e))
```

```
## [1] "1367452_at" "1367453_at" "1367454_at" "1367455_at" "1367456_at"
## [6] "1367457_at"
```

```r
annotation(e)
```

```
## [1] "rae230a"
```


The annotation for this ExpressionSet is *rae230a*. Many platforms will have database annotation packages already existing on Bioconductor. We can access these, first by installing, and then loading the library. We will use the `AnnotationDbi` package to query the information in the library.

While in this unit we will use a microarray annotation package as an example, the same commands can be used for an organism package, such as the homo sapiens annotation package `org.Hs.eg.db`, which let's one query from one kind of gene annotation to another.


```r
# biocLite(paste0(annotation(e),'.db'))
library(rae230a.db)
```

```
## Loading required package: AnnotationDbi
## Loading required package: GenomeInfoDb
## Loading required package: org.Rn.eg.db
## Loading required package: DBI
```

```r
# biocLite('AnnotationDbi')
library(AnnotationDbi)
```


Annotation packages have *columns*, some of which may be *keys*. You can query the database using a *key*, and ask for one or more *columns* in return. We will use the rownames of the ExpressionSet as keys.


```r
columns(rae230a.db)
```

```
##  [1] "PROBEID"      "ENTREZID"     "PFAM"         "IPI"         
##  [5] "PROSITE"      "ACCNUM"       "ALIAS"        "CHR"         
##  [9] "CHRLOC"       "CHRLOCEND"    "ENZYME"       "PATH"        
## [13] "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"     
## [17] "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"    
## [21] "UNIPROT"      "GO"           "EVIDENCE"     "ONTOLOGY"    
## [25] "GOALL"        "EVIDENCEALL"  "ONTOLOGYALL"
```

```r
keytypes(rae230a.db)
```

```
##  [1] "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
##  [5] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"      
##  [9] "CHRLOCEND"    "ENZYME"       "PATH"         "PMID"        
## [13] "REFSEQ"       "SYMBOL"       "UNIGENE"      "ENSEMBL"     
## [17] "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"     "UNIPROT"     
## [21] "GO"           "EVIDENCE"     "ONTOLOGY"     "GOALL"       
## [25] "EVIDENCEALL"  "ONTOLOGYALL"  "PROBEID"
```

```r
head(keys(rae230a.db, keytype = "PROBEID"))
```

```
## [1] "1367452_at" "1367453_at" "1367454_at" "1367455_at" "1367456_at"
## [6] "1367457_at"
```

```r
head(rownames(e))
```

```
## [1] "1367452_at" "1367453_at" "1367454_at" "1367455_at" "1367456_at"
## [6] "1367457_at"
```


The following `select` call will return the Entrez ID, ENSEMBL ID, and gene symbol for each Probe ID, which are the rownames of the ExpressionSet.


```r
res <- select(rae230a.db, keys = rownames(e), columns = c("ENTREZID", "ENSEMBL", 
    "SYMBOL"), keytype = "PROBEID")
```

```
## Warning: 'select' resulted in 1:many mapping between keys and return rows
```

```r
head(res)
```

```
##      PROBEID ENTREZID            ENSEMBL SYMBOL
## 1 1367452_at   690244 ENSRNOG00000032840  Sumo2
## 2 1367453_at   114562 ENSRNOG00000033426  Cdc37
## 3 1367454_at    60384 ENSRNOG00000045871  Copb2
## 4 1367455_at   116643 ENSRNOG00000034242    Vcp
## 5 1367456_at    81920 ENSRNOG00000013741 Ube2d3
## 6 1367456_at   641452 ENSRNOG00000013741 Ube2d2
```

```r
idx <- match(rownames(e), res$PROBEID)
```


We need to align the `res` object so that we pull out, in order, one row for each row of the ExpressionSet.


```r
head(rownames(e))
```

```
## [1] "1367452_at" "1367453_at" "1367454_at" "1367455_at" "1367456_at"
## [6] "1367457_at"
```

```r
head(res$PROBEID, 7)
```

```
## [1] "1367452_at" "1367453_at" "1367454_at" "1367455_at" "1367456_at"
## [6] "1367456_at" "1367457_at"
```

```r
head(idx)
```

```
## [1] 1 2 3 4 5 7
```


Here we add the new information to the `fData` of `e`. If there were already information in `fData`, we would have used `cbind` to add the new columns. Note here that, since we have a one-to-many mapping, the `match` function gave us the first match that it found. You could also collapse all possible matches of the Probe ID to the Genes using `split` and `paste` with the `collapse` argument. However, here we keep it simple and just take the first match in the `res` object.


```r
fData(e) <- res[idx, ]
head(fData(e), 10)
```

```
##       PROBEID ENTREZID            ENSEMBL SYMBOL
## 1  1367452_at   690244 ENSRNOG00000032840  Sumo2
## 2  1367453_at   114562 ENSRNOG00000033426  Cdc37
## 3  1367454_at    60384 ENSRNOG00000045871  Copb2
## 4  1367455_at   116643 ENSRNOG00000034242    Vcp
## 5  1367456_at    81920 ENSRNOG00000013741 Ube2d3
## 7  1367457_at   114558 ENSRNOG00000020513  Becn1
## 8  1367458_at    83510 ENSRNOG00000010067 Lypla2
## 9  1367459_at    64310 ENSRNOG00000030591   Arf1
## 10 1367460_at    29662 ENSRNOG00000018091   Gdi2
## 11 1367461_at   114023 ENSRNOG00000012000  Copb1
```

```r
all.equal(fData(e)$PROBEID, rownames(e))
```

```
## [1] TRUE
```


## Using Biomart

An alternate way to map from one annotation to another is using the `biomaRt` package. For more information on which Biomarts are available and how to access them, see the `biomaRt` vignette.


```r
# biocLite('biomaRt')
library(biomaRt)
# vignette('biomaRt')
m <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
map <- getBM(mart = m, attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", 
    values = fData(e)$ENSEMBL)
head(map)
```

```
##      ensembl_gene_id entrezgene
## 1 ENSRNOG00000000007      24379
## 2 ENSRNOG00000000010     498922
## 3 ENSRNOG00000000024     362454
## 4 ENSRNOG00000000028      83836
## 5 ENSRNOG00000000029      24582
## 6 ENSRNOG00000000033     305095
```


Finally, we need to align the new information with the old information using the `match` function as before, again picking the first match from a one-to-many mapping. We see that for the most part the new and the old Entrez IDs are the same, though some differences occur when we pick one from the one-to-many mappings that exist.



```r
idx <- match(fData(e)$ENSEMBL, map$ensembl_gene_id)
fData(e)$NEW_ENTREZID <- map$entrezgene[idx]
head(fData(e))
```

```
##      PROBEID ENTREZID            ENSEMBL SYMBOL NEW_ENTREZID
## 1 1367452_at   690244 ENSRNOG00000032840  Sumo2       682787
## 2 1367453_at   114562 ENSRNOG00000033426  Cdc37       114562
## 3 1367454_at    60384 ENSRNOG00000045871  Copb2        60384
## 4 1367455_at   116643 ENSRNOG00000034242    Vcp       116643
## 5 1367456_at    81920 ENSRNOG00000013741 Ube2d3        81920
## 7 1367457_at   114558 ENSRNOG00000020513  Becn1       114558
```

```r
mean(fData(e)$ENTREZID == fData(e)$NEW_ENTREZID, na.rm = TRUE)
```

```
## [1] 0.993
```


