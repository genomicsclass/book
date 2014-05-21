---
layout: page
title: Using limma for microarray analysis
---




## First, simple t-tests

In this unit, we will show the difference between using the simple t-test and doing differential expression with the `limma` hierarchical model. The reference is [Smyth 2004](#foot), listed in the footnotes.

Here we also show the basic steps for performing a `limma` analysis. Note that the `limma` package is very powerful, and has hundreds of pages of documentation which we cannot cover in this course, however we recommend that users wanting to explore further should check out this guide.

We start by loading the spike-in data which was introduced in lecture, which has already been normalized.


```r
# biocLite('SpikeInSubset')
library(SpikeInSubset)
```

```
## Loading required package: Biobase
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
## 
## Loading required package: affy
```

```r
data(rma95)
fac <- factor(rep(1:2, each = 3))
```


We can now perform simple t-tests using the `rowttests` function in the `genefilter` package:


```r
library(genefilter)
```

```
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
rtt <- rowttests(exprs(rma95), fac)
```


We will define colors depending on whether the p-value is small, the absolute difference in means is large, and whether the feature is a spike-in value.


```r
mask <- with(rtt, abs(dm) < 0.2 & p.value < 0.01)
spike <- rownames(rma95) %in% colnames(pData(rma95))
cols <- ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))
```


We now plot the results, using the colors defined above. We multiply the `dm` by -1, because we are interested in the difference from the second group to the first (this is the difference used by `lm` and the `limma` package by default). The spike-in genes are in blue, which have mostly small p-value and large difference in means. The red points indicate genes which have small p-values but also small differences in means. We will see how these points change after using `limma`.


```r
with(rtt, plot(-dm, -log10(p.value), cex = 0.8, pch = 16, xlim = c(-1, 1), ylim = c(0, 
    5), xlab = "difference in means", col = cols))
abline(h = 2, v = c(-0.2, 0.2), lty = 2)
```

![plot of chunk unnamed-chunk-4](figure/using_limma-unnamed-chunk-4.png) 


Note that the red genes have mostly low estimates of standard deviation.


```r
rtt$s <- apply(exprs(rma95), 1, function(row) sqrt(0.5 * (var(row[1:3]) + var(row[4:6]))))
with(rtt, plot(s, -log10(p.value), cex = 0.8, pch = 16, log = "x", xlab = "estimate of standard deviation", 
    col = cols))
```

![plot of chunk unnamed-chunk-5](figure/using_limma-unnamed-chunk-5.png) 


## limma steps

The following three steps perform the basic `limma` analysis. We specify `coef=2` because we are interested in the difference between groups, not the intercept.


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
fit <- lmFit(rma95, design = model.matrix(~fac))
colnames(coef(fit))
```

```
## [1] "(Intercept)" "fac2"
```

```r
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2)
tt
```

```
##            logFC AveExpr       t   P.Value adj.P.Val     B
## 1708_at  -7.0611   7.945 -73.529 7.816e-17 9.869e-13 8.647
## 36202_at  0.8526   9.373   9.975 4.936e-07 3.116e-03 4.588
## 36311_at  0.8318   8.564   8.363 3.017e-06 1.270e-02 3.568
## 33264_at  0.7119   4.919   7.435 9.666e-06 2.707e-02 2.836
## 32660_at  0.6554   8.680   7.356 1.072e-05 2.707e-02 2.768
## 38734_at  0.7467   6.256   7.185 1.345e-05 2.831e-02 2.618
## 1024_at   0.8427   9.697   6.731 2.503e-05 4.400e-02 2.196
## 36085_at  0.6449  12.193   6.654 2.788e-05 4.400e-02 2.121
## 33818_at  0.5322  12.286   6.455 3.699e-05 5.190e-02 1.923
## 39058_at  0.6091   7.535   6.279 4.768e-05 5.688e-02 1.743
```


`topTable` will return the top genes ranked by whichever value you define. You can also ask topTable to return all the values, sorted by `"none"`. Note that a column automatically is included which gives the *adjusted p-values* for each gene. By default the method of Benjamini-Hochberg is used, by calling the `p.adjust` function.


```r
# ?topTable
dim(topTable(fit, coef = 2, number = Inf, sort.by = "none"))
```

```
## [1] 12626     6
```

```r
# ?p.adjust
```


Here we will compare the previous volcano plot with the `limma` results. Note that the red points are now all under the line where `-log10(p.value)` is equal to 2. Also, the blue points which represent real differences have p-values which are even higher than before.


```r
limmares <- data.frame(dm = coef(fit)[, "fac2"], p.value = fit$p.value[, "fac2"])
with(limmares, plot(dm, -log10(p.value), cex = 0.8, pch = 16, col = cols, xlab = "difference in means", 
    xlim = c(-1, 1), ylim = c(0, 5)))
abline(h = 2, v = c(-0.2, 0.2), lty = 2)
```

![plot of chunk unnamed-chunk-8](figure/using_limma-unnamed-chunk-8.png) 


Finally, we will construct a plot which shows how `limma` shrinks the variance estimates towards a common value, eliminating false positives which might arise from too-low estimates of variance.

Here we pick, for each of 40 bins of different variance estimates, a single gene which falls in that bin. We remove bins which do not have any such genes.


```r
n <- 40
qs <- seq(from = 0, to = 0.2, length = n)
idx <- sapply(seq_len(n), function(i) which(as.integer(cut(rtt$s^2, qs)) == 
    i)[1])
idx <- idx[!is.na(idx)]
```


Now we will plot a line, from the initial estimate of variance for these genes to the estimate after running `limma`.


```r
par(mar = c(5, 5, 2, 2))
plot(1, 1, xlim = c(0, 0.21), ylim = c(0, 1), type = "n", xlab = "variance estimates", 
    ylab = "", yaxt = "n")
axis(2, at = c(0.1, 0.9), c("before", "after"), las = 2)
segments((rtt$s^2)[idx], rep(0.1, n), fit$s2.post[idx], rep(0.9, n))
```

![plot of chunk unnamed-chunk-10](figure/using_limma-unnamed-chunk-10.png) 


## Footnotes <a name="foot"></a>

Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>
