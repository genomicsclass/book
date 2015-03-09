---
layout: page
title: Distance lecture
---



# Introduction

The concept of distance can be generalized from  physical distance. For example, we cluster animals into groups. When we do this, we put animals that "close" in the same group:

<img src="images/animals.jpg" align="middle" width=300>

Any time we cluster individuals into separate groups we are, explicitely or implicitely computing a distance. 

Do create _heatmaps_ a distance is computed explicitely. Heatmaps are widely used in genomics and other highthroughput fields:

<img src="images/Heatmap.png" align="middle" width=300>
Image Source: Heatmap, Gaeddal, 01.28.2007, http://commons.wikimedia.org/wiki/File:Heatmap.png, PD

In these plots the measurements, which are stored ina matrix, are represented with colors after the columns and rows have been clustered. Here we will learn the necessary mathematics and computing skill to understand and create heatmaps. We start by reviewing the mathematical definition of distance. 


# Euclidean Distance

As a review, let's define the distance between two points, $$A$$ and $$B$$, on a cartesian plane.


```
## Loading required package: RColorBrewer
```

<img src="figure/distance-unnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />

The euclidean distance between A and B is simply

$$\sqrt{ (A_x-B_x)^2 + (A_y-B_y)^2}$$

# High dimensional Data

In this chapter we focus on high-dimensional data. We introduce a data set with gene expression measurements for 22215 genes from 189 samples. The R ojects can be downloaded like this:


```r
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
```

The data represent RNA expression levels for eight tissues, each with several individuals.


```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
table(tissue)
```

```
## tissue
##  cerebellum       colon endometrium hippocampus      kidney       liver 
##          38          34          15          31          39          26 
##    placenta 
##           6
```

# Distance in high dimensions

We are interested in describing distance in the context of this dataset. We might also be interested in finding genes that _behave similarly_ across samples.

To define distance we need to know what points are since distance is computed between points. With high dimensional data, points are no longer on the cartesian plan. Instead they are in higher dimensions. For exampe, sample $$i$$ is defined by the point in 22215 dimesions $$(Y_{1,i},\dots,Y_{22215,i})'$$. Feature $$g$$ feature $$g$$ is defined by the point in 189 dimensions $$(Y_{g,189},\dots,Y_{g,189})'$$

Once we define points, the Euclidean distance is defined in a very similar way as it is defined for two dimensions. For example, the  distance between two samples $$i$$ and $$j$$ is

$$
d(i,j) = \sqrt{ \sum_{g=1}^{22215} (Y_{g,i}-Y_{g,j })^2 }
$$

and the distance between two features $$h$$ and $$g$$ as:
$$
d(h,g) = \sqrt{ \sum_{i=1}^{189} (Y_{h,i}-Y_{g,i})^2 }
$$


# Distance with Matrix Algebra

The distance between samples $$i$$ and $$j$$ can be written as


$$ d(i,j) = (\mathbf{Y}_i - \mathbf{Y}_j)^\top(\mathbf{Y}_i - \mathbf{Y}_j)$$

With $$\mathbf{Y}_i$$ and $$\mathbf{Y}_j$$ coliumns $$i$$ and $$j$$


# Examples

We can now use the formulas above to compute distance. Let's compute distance between samples 1 and 2, both kidneys, and then to 87, a colon.


```r
x <- e[,1]
y <- e[,2]
z <- e[,87]
sqrt(sum((x-y)^2))
```

```
## [1] 85.8546
```

```r
sqrt(sum((x-z)^2))
```

```
## [1] 122.8919
```

As expeceted the kidneys are closer to each other. A faster way to compute this is using matrix algebra


```r
sqrt( crossprod(x-y) )
```

```
##         [,1]
## [1,] 85.8546
```

```r
sqrt( crossprod(x-z) )
```

```
##          [,1]
## [1,] 122.8919
```

Now to compute all the distances at once we have the function `dist`. Because it computes the distance between each row, and here we are interested in the distance between samples we transpose the matrix


```r
d <- dist(t(e))
class(d)
```

```
## [1] "dist"
```

Note that this produces the an object of class `dist` and to access to entries we need to coerce into a matrix:


```r
as.matrix(d)[1,2]
```

```
## [1] 85.8546
```

```r
as.matrix(d)[1,87]
```

```
## [1] 122.8919
```

It is important to keep in mind that if we run `dist` on `e` it will compute all pairwise distances between genes. This will try to create a $$22215 \times 22215$$ matrix that may kill crash your R sessions.


