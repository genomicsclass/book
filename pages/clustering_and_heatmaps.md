---
layout: page
title: Clustering lab
---



# Introduction

We load the tissue gene expression data


```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
```

Let's act as if we don't know these are different tissues and are interested in clustering. The first step is to compute the distance between each sample


```r
d <- dist(t(e))
```

<a name="hierarchical"></a>

# Hierarchical clustering

We can perform hierarchical clustering based on the
distances defined above, using the `hclust` function.
The `plot` method will make a plot of the tree that results from `hclust`. 


```r
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
mypar2(1,1)
hc <- hclust(d)
hc
```

```
## 
## Call:
## hclust(d = d)
## 
## Cluster method   : complete 
## Distance         : euclidean 
## Number of objects: 189
```

```r
plot(hc,labels=tissue)
```

![plot of chunk unnamed-chunk-3](figure/clustering_and_heatmaps-unnamed-chunk-3-1.png) 

Does this technique "discover" the clusters defined by the different tissues? In this case it is not easy to see the different tissues so we add colors.
 

```r
myplclust(hc, labels=tissue, lab.col=as.fumeric(tissue))
```

![plot of chunk unnamed-chunk-4](figure/clustering_and_heatmaps-unnamed-chunk-4-1.png) 

Note that hierarchical clustering does not define specific clusters but rather defines a dendrogram. To define clusters we need to "cut the tree" at some height and group all samples that are within that split into different groups below that line. We use 120 as an example:


```r
myplclust(hc, labels=tissue, lab.col=as.fumeric(tissue))
abline(h=120)
```

![plot of chunk unnamed-chunk-5](figure/clustering_and_heatmaps-unnamed-chunk-5-1.png) 

If we use the line above to cut the tree into clusters, we can examine how the clusters overlap with the actual tissues:


```r
hclusters <- cutree(hc, h=120)
table(true=tissue, cluster=hclusters)
```

```
##              cluster
## true           1  2  3  4  5  6  7  8  9 10 11 12 13 14
##   cerebellum   0  0  0  0 31  0  0  0  2  0  0  5  0  0
##   colon        0  0  0  0  0  0 34  0  0  0  0  0  0  0
##   endometrium  0  0  0  0  0  0  0  0  0  0 15  0  0  0
##   hippocampus  0  0 12 19  0  0  0  0  0  0  0  0  0  0
##   kidney       9 18  0  0  0 10  0  0  2  0  0  0  0  0
##   liver        0  0  0  0  0  0  0 24  0  2  0  0  0  0
##   placenta     0  0  0  0  0  0  0  0  0  0  0  0  2  4
```

<a name="kmeans"></a>

# K-means

We can also cluster with the `kmeans` function to perform k-means clustering. As an example, let's run k-means on the samples in the space of the first two genes:


```r
plot(e[1,], e[2,])
```

![plot of chunk unnamed-chunk-7](figure/clustering_and_heatmaps-unnamed-chunk-7-1.png) 

```r
set.seed(1)
km <- kmeans(t(e[1:2,]), centers=7)
names(km)
```

```
## [1] "cluster"      "centers"      "totss"        "withinss"    
## [5] "tot.withinss" "betweenss"    "size"         "iter"        
## [9] "ifault"
```

```r
plot(e[1,], e[2,], col=km$cluster, pch=16)
```

![plot of chunk unnamed-chunk-7](figure/clustering_and_heatmaps-unnamed-chunk-7-2.png) 

```r
plot(e[1,], e[2,], col=as.fumeric(tissue), pch=16)
```

![plot of chunk unnamed-chunk-7](figure/clustering_and_heatmaps-unnamed-chunk-7-3.png) 

```r
table(true=tissue,cluster=km$cluster)
```

```
##              cluster
## true           1  2  3  4  5  6  7
##   cerebellum   0  1  8  0  6  0 23
##   colon        2 11  2 15  4  0  0
##   endometrium  0  3  4  0  0  0  8
##   hippocampus 19  0  2  0 10  0  0
##   kidney       7  8 20  0  0  0  4
##   liver        0  0  0  0  0 18  8
##   placenta     0  4  0  0  0  0  2
```

We can instead perform k-means clustering using all of the genes. And to visualize this, we can use an MDS plot



```r
mds <- cmdscale(d)
plot(mds[,1], mds[,2]) 
```

![plot of chunk unnamed-chunk-8](figure/clustering_and_heatmaps-unnamed-chunk-8-1.png) 

```r
km <- kmeans(t(e), centers=7)
plot(mds[,1], mds[,2], col=km$cluster, pch=16)
```

![plot of chunk unnamed-chunk-8](figure/clustering_and_heatmaps-unnamed-chunk-8-2.png) 

```r
table(true=tissue,cluster=km$cluster)
```

```
##              cluster
## true           1  2  3  4  5  6  7
##   cerebellum   0  0  5  0 31  2  0
##   colon        0 34  0  0  0  0  0
##   endometrium  0 15  0  0  0  0  0
##   hippocampus  0  0 31  0  0  0  0
##   kidney       0 37  0  0  0  2  0
##   liver        2  0  0  0  0  0 24
##   placenta     0  0  0  6  0  0  0
```


<a name="heatmap"></a>

# Heatmaps

Heatmaps are useful plots for visualizing the measurements  for a subset of rows over all the samples. A *dendrogram* is added on top and on the side is a hierarchical clustering as we saw before. First we will
use the `heatmap` available in base R. First define a color palette.


```r
# install.packages("RColorBrewer")
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
```

Now, pick the genes with the top variance over all samples:


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
rv <- rowVars(e)
idx <- order(-rv)[1:40]
```

Now we can plot a heatmap of these genes:


```r
heatmap(e[idx,], col=hmcol)
```

![plot of chunk unnamed-chunk-11](figure/clustering_and_heatmaps-unnamed-chunk-11-1.png) 

The `heatmap.2` function in the `gplots` package on CRAN is a bit more
customization. For example, it stretches to fill the window. Here we add colors to indicate the tissue on the top:


```r
# install.packages("gplots")
library(gplots)
```

```
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]
cbind(colnames(e),cols)
```

```
##                           cols     
##   [1,] "GSM11805.CEL.gz"  "#1B9E77"
##   [2,] "GSM11814.CEL.gz"  "#1B9E77"
##   [3,] "GSM11823.CEL.gz"  "#1B9E77"
##   [4,] "GSM11830.CEL.gz"  "#1B9E77"
##   [5,] "GSM12067.CEL.gz"  "#1B9E77"
##   [6,] "GSM12075.CEL.gz"  "#1B9E77"
##   [7,] "GSM12079.CEL.gz"  "#1B9E77"
##   [8,] "GSM12098.CEL.gz"  "#1B9E77"
##   [9,] "GSM12105.CEL.gz"  "#1B9E77"
##  [10,] "GSM12268.CEL.gz"  "#1B9E77"
##  [11,] "GSM12270.CEL.gz"  "#1B9E77"
##  [12,] "GSM12283.CEL.gz"  "#1B9E77"
##  [13,] "GSM12298.CEL.gz"  "#1B9E77"
##  [14,] "GSM12300.CEL.gz"  "#1B9E77"
##  [15,] "GSM12399.CEL.gz"  "#1B9E77"
##  [16,] "GSM12444.CEL.gz"  "#1B9E77"
##  [17,] "GSM21203.cel.gz"  "#D95F02"
##  [18,] "GSM21204.cel.gz"  "#D95F02"
##  [19,] "GSM21205.cel.gz"  "#D95F02"
##  [20,] "GSM21206.cel.gz"  "#D95F02"
##  [21,] "GSM21207.cel.gz"  "#D95F02"
##  [22,] "GSM21208.cel.gz"  "#D95F02"
##  [23,] "GSM21209.cel.gz"  "#D95F02"
##  [24,] "GSM21210.cel.gz"  "#D95F02"
##  [25,] "GSM21212.cel.gz"  "#D95F02"
##  [26,] "GSM21213.cel.gz"  "#D95F02"
##  [27,] "GSM21214.cel.gz"  "#D95F02"
##  [28,] "GSM21215.cel.gz"  "#D95F02"
##  [29,] "GSM21216.cel.gz"  "#D95F02"
##  [30,] "GSM21217.cel.gz"  "#D95F02"
##  [31,] "GSM21218.cel.gz"  "#D95F02"
##  [32,] "GSM21219.cel.gz"  "#D95F02"
##  [33,] "GSM21220.cel.gz"  "#D95F02"
##  [34,] "GSM21221.cel.gz"  "#D95F02"
##  [35,] "GSM21222.cel.gz"  "#D95F02"
##  [36,] "GSM21223.cel.gz"  "#D95F02"
##  [37,] "GSM21224.cel.gz"  "#D95F02"
##  [38,] "GSM21225.cel.gz"  "#D95F02"
##  [39,] "GSM21226.cel.gz"  "#D95F02"
##  [40,] "GSM21227.cel.gz"  "#D95F02"
##  [41,] "GSM21228.cel.gz"  "#D95F02"
##  [42,] "GSM21229.cel.gz"  "#D95F02"
##  [43,] "GSM21230.cel.gz"  "#D95F02"
##  [44,] "GSM21231.cel.gz"  "#D95F02"
##  [45,] "GSM21232.cel.gz"  "#D95F02"
##  [46,] "GSM21233.cel.gz"  "#D95F02"
##  [47,] "GSM87058.cel.gz"  "#7570B3"
##  [48,] "GSM87064.cel.gz"  "#7570B3"
##  [49,] "GSM87065.cel.gz"  "#7570B3"
##  [50,] "GSM87066.cel.gz"  "#7570B3"
##  [51,] "GSM87071.cel.gz"  "#7570B3"
##  [52,] "GSM87072.cel.gz"  "#7570B3"
##  [53,] "GSM87073.cel.gz"  "#7570B3"
##  [54,] "GSM87074.cel.gz"  "#7570B3"
##  [55,] "GSM87075.cel.gz"  "#7570B3"
##  [56,] "GSM87085.cel.gz"  "#7570B3"
##  [57,] "GSM87086.cel.gz"  "#7570B3"
##  [58,] "GSM87087.cel.gz"  "#7570B3"
##  [59,] "GSM87088.cel.gz"  "#7570B3"
##  [60,] "GSM87089.cel.gz"  "#7570B3"
##  [61,] "GSM87090.CEL.gz"  "#7570B3"
##  [62,] "GSM87091.cel.gz"  "#7570B3"
##  [63,] "GSM87092.cel.gz"  "#7570B3"
##  [64,] "GSM87093.cel.gz"  "#7570B3"
##  [65,] "GSM87094.cel.gz"  "#7570B3"
##  [66,] "GSM87095.cel.gz"  "#7570B3"
##  [67,] "GSM87096.cel.gz"  "#7570B3"
##  [68,] "GSM87097.cel.gz"  "#7570B3"
##  [69,] "GSM87098.cel.gz"  "#7570B3"
##  [70,] "GSM87099.cel.gz"  "#7570B3"
##  [71,] "GSM87100.cel.gz"  "#7570B3"
##  [72,] "GSM87101.cel.gz"  "#7570B3"
##  [73,] "GSM87102.cel.gz"  "#7570B3"
##  [74,] "GSM146778.CEL.gz" "#1B9E77"
##  [75,] "GSM146779.CEL.gz" "#1B9E77"
##  [76,] "GSM146781.CEL.gz" "#1B9E77"
##  [77,] "GSM146782.CEL.gz" "#1B9E77"
##  [78,] "GSM146783.CEL.gz" "#1B9E77"
##  [79,] "GSM146785.CEL.gz" "#1B9E77"
##  [80,] "GSM146787.CEL.gz" "#1B9E77"
##  [81,] "GSM146788.CEL.gz" "#1B9E77"
##  [82,] "GSM146789.CEL.gz" "#1B9E77"
##  [83,] "GSM146791.CEL.gz" "#1B9E77"
##  [84,] "GSM146793.CEL.gz" "#1B9E77"
##  [85,] "GSM146795.CEL.gz" "#1B9E77"
##  [86,] "GSM146797.CEL.gz" "#1B9E77"
##  [87,] "GSM92240.CEL.gz"  "#E7298A"
##  [88,] "GSM92241.CEL.gz"  "#E7298A"
##  [89,] "GSM92242.CEL.gz"  "#E7298A"
##  [90,] "GSM92243.CEL.gz"  "#E7298A"
##  [91,] "GSM92244.CEL.gz"  "#E7298A"
##  [92,] "GSM92245.CEL.gz"  "#E7298A"
##  [93,] "GSM92247.CEL.gz"  "#E7298A"
##  [94,] "GSM92248.CEL.gz"  "#E7298A"
##  [95,] "GSM92249.CEL.gz"  "#E7298A"
##  [96,] "GSM92250.CEL.gz"  "#E7298A"
##  [97,] "GSM92253.CEL.gz"  "#E7298A"
##  [98,] "GSM92254.CEL.gz"  "#E7298A"
##  [99,] "GSM92255.CEL.gz"  "#E7298A"
## [100,] "GSM92256.CEL.gz"  "#E7298A"
## [101,] "GSM92257.CEL.gz"  "#E7298A"
## [102,] "GSM92258.CEL.gz"  "#E7298A"
## [103,] "GSM92259.CEL.gz"  "#E7298A"
## [104,] "GSM92260.CEL.gz"  "#E7298A"
## [105,] "GSM92261.CEL.gz"  "#E7298A"
## [106,] "GSM92262.CEL.gz"  "#E7298A"
## [107,] "GSM92263.CEL.gz"  "#E7298A"
## [108,] "GSM92264.CEL.gz"  "#E7298A"
## [109,] "GSM92265.CEL.gz"  "#E7298A"
## [110,] "GSM92266.CEL.gz"  "#E7298A"
## [111,] "GSM92267.CEL.gz"  "#E7298A"
## [112,] "GSM92268.CEL.gz"  "#E7298A"
## [113,] "GSM92269.CEL.gz"  "#E7298A"
## [114,] "GSM92270.CEL.gz"  "#E7298A"
## [115,] "GSM92271.CEL.gz"  "#E7298A"
## [116,] "GSM92272.CEL.gz"  "#E7298A"
## [117,] "GSM92273.CEL.gz"  "#E7298A"
## [118,] "GSM92274.CEL.gz"  "#E7298A"
## [119,] "GSM92275.CEL.gz"  "#E7298A"
## [120,] "GSM92276.CEL.gz"  "#E7298A"
## [121,] "GSM35979.cel.gz"  "#1B9E77"
## [122,] "GSM35980.cel.gz"  "#1B9E77"
## [123,] "GSM35981.cel.gz"  "#1B9E77"
## [124,] "GSM35982.cel.gz"  "#66A61E"
## [125,] "GSM35983.cel.gz"  "#66A61E"
## [126,] "GSM35984.cel.gz"  "#66A61E"
## [127,] "GSM35991.cel.gz"  "#1B9E77"
## [128,] "GSM35992.cel.gz"  "#1B9E77"
## [129,] "GSM35993.cel.gz"  "#1B9E77"
## [130,] "GSM35994.cel.gz"  "#66A61E"
## [131,] "GSM35995.cel.gz"  "#66A61E"
## [132,] "GSM35996.cel.gz"  "#66A61E"
## [133,] "GSM36003.cel.gz"  "#1B9E77"
## [134,] "GSM44675.CEL.gz"  "#1B9E77"
## [135,] "GSM44689.CEL.gz"  "#7570B3"
## [136,] "GSM44697.CEL.gz"  "#D95F02"
## [137,] "GSM44702.CEL.gz"  "#66A61E"
## [138,] "GSM181429.CEL.gz" "#66A61E"
## [139,] "GSM181430.CEL.gz" "#66A61E"
## [140,] "GSM181431.CEL.gz" "#66A61E"
## [141,] "GSM181432.CEL.gz" "#66A61E"
## [142,] "GSM181433.CEL.gz" "#66A61E"
## [143,] "GSM18917.CEL.gz"  "#7570B3"
## [144,] "GSM18918.CEL.gz"  "#7570B3"
## [145,] "GSM18953.CEL.gz"  "#66A61E"
## [146,] "GSM18954.CEL.gz"  "#66A61E"
## [147,] "GSM18955.CEL.gz"  "#1B9E77"
## [148,] "GSM18956.CEL.gz"  "#1B9E77"
## [149,] "GSM296875.CEL.gz" "#E6AB02"
## [150,] "GSM296876.CEL.gz" "#E6AB02"
## [151,] "GSM296878.CEL.gz" "#E6AB02"
## [152,] "GSM296879.CEL.gz" "#E6AB02"
## [153,] "GSM296880.CEL.gz" "#E6AB02"
## [154,] "GSM296881.CEL.gz" "#E6AB02"
## [155,] "GSM296882.CEL.gz" "#E6AB02"
## [156,] "GSM296883.CEL.gz" "#E6AB02"
## [157,] "GSM296886.CEL.gz" "#E6AB02"
## [158,] "GSM296887.CEL.gz" "#E6AB02"
## [159,] "GSM296888.CEL.gz" "#E6AB02"
## [160,] "GSM296889.CEL.gz" "#E6AB02"
## [161,] "GSM296890.CEL.gz" "#E6AB02"
## [162,] "GSM296891.CEL.gz" "#E6AB02"
## [163,] "GSM296892.CEL.gz" "#E6AB02"
## [164,] "GSM298747.CEL.gz" "#66A61E"
## [165,] "GSM298748.CEL.gz" "#66A61E"
## [166,] "GSM298749.CEL.gz" "#66A61E"
## [167,] "GSM298750.CEL.gz" "#66A61E"
## [168,] "GSM299110.CEL.gz" "#66A61E"
## [169,] "GSM299111.CEL.gz" "#66A61E"
## [170,] "GSM299112.CEL.gz" "#66A61E"
## [171,] "GSM299113.CEL.gz" "#66A61E"
## [172,] "GSM299244.CEL.gz" "#66A61E"
## [173,] "GSM299245.CEL.gz" "#66A61E"
## [174,] "GSM299246.CEL.gz" "#66A61E"
## [175,] "GSM299247.CEL.gz" "#66A61E"
## [176,] "GSM322969.CEL.gz" "#7570B3"
## [177,] "GSM323054.CEL.gz" "#7570B3"
## [178,] "GSM323523.CEL.gz" "#7570B3"
## [179,] "GSM323524.CEL.gz" "#7570B3"
## [180,] "GSM323527.CEL.gz" "#7570B3"
## [181,] "GSM323565.CEL.gz" "#7570B3"
## [182,] "GSM323566.CEL.gz" "#7570B3"
## [183,] "GSM323567.CEL.gz" "#7570B3"
## [184,] "GSM246492.CEL.gz" "#A6761D"
## [185,] "GSM246493.CEL.gz" "#A6761D"
## [186,] "GSM246494.CEL.gz" "#A6761D"
## [187,] "GSM307639.CEL.gz" "#A6761D"
## [188,] "GSM307640.CEL.gz" "#A6761D"
## [189,] "GSM307641.CEL.gz" "#A6761D"
```

```r
heatmap.2(e[idx,], labCol=tissue,
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)
```

![plot of chunk unnamed-chunk-12](figure/clustering_and_heatmaps-unnamed-chunk-12-1.png) 


