---
layout: page
Title: EDA for Highthroughput Exercises
---

{pagebreak} 

## Exercises

We will be using a handful of Bioconductor packages. These are installed using the function `biocLite` which you can source from the web:


```r
source("http://www.bioconductor.org/biocLite.R")
```

or you can run the `bioc_install` in the `rafalib` package.


```r
library(rafalib)
bioc_install()
```


Download and install the Bioconductor package `SpikeInSubset` and then load the library and the `mas133` data:


```r
library(rafalib)
install_bioc("SpikeInSubset")
library(SpikeInSubset)
data(mas133)
```

Now make the following plot of the first two samples and compute the correlation:

```r
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")
```

1. What proportion of the points are inside the box?




2. Now make the sample plot with log:

    
    ```r
    plot(log2(e[,1]),log2(e[,2]))
    k <- log2(3000)
    b <- log2(0.5)
    polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")
    ```

    What is an advantage of taking the log?
    
    - A) The tails do not dominate the plot: 95% of data is not in a tiny section of plot.
    - B) There are less points.
    - C) There is exponential growth.
    - D) We always take logs.



3. Make an MA-plot:

    
    ```r
    e <- log2(exprs(mas133))
    plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)
    ```

    The two samples we are plotting are replicates (they are random samples from the same batch of RNA). The correlation of the data was 0.997 in the original scale and 0.96 in the log-scale. High correlations are sometimes confused with evidence of replication. However, replication implies we get very small differences between the observations, which is better measured with distance or differences.

    What is the standard deviation of the log ratios for this comparison? 



4. How many fold changes above 2 do we see?



