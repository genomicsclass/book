---
Title: EDA for Highthroughput Exercises
---

A>## Exercises
A>
A>We will be using a handful of Bioconductor packages. These are installed using the function `biocLite` which you can source from the web
A>
A>
```r
source("http://www.bioconductor.org/biocLite.R")
```
A>
A>or you can run the `bioc_install` in the `rafalib` package.
A>
A>
```r
library(rafalib)
bioc_install()
```
A>
A>
A>Download and install the Bioconductor package SpikeInSubset and then load the library and the mas133 data:
A>
A>
```r
library(rafalib)
install_bioc("SpikeInSubset")
library(SpikeInSubset)
data(mas133)
```
A>
A>Now make the following plot of the first two samples and compute the correlation:
A>
```r
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")
```
A>
A>1. What proportion of the points are inside the box?
A>
A>
A>
A>
A>2. Now make the sample plot with log:
A>
A>    
    ```r
    plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
    k <- log2(3000)
    b <- log2(0.5)
    polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")
    ```
A>
A>    What is an advantage of taking the log?
A>      - A) The tails do not dominate the plot: 95% of data is not in a tiny section of plot.
A>      - B) There are less points.
A>      - C) There is exponential growth.
A>      - D) We always take logs
A>
A>
A>
A>3. Make an MA-plot
A>
A>    
    ```r
    e <- log2(exprs(mas133))
    plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)
    ```
A>
A>    The two samples we are plotting are replicates (they random samples fro the same batch of RNA) The correlation of the data was 0.997 in the original scale, 0.96 in the log-scale. High correlations are sometimes confused with evidence of replication. But replication implies we get very small difference between the observations which is better measured with distance or differences.
A>
A>    What is the standard deviation of the log ratios for this comparison? 
A>
A>
A>
A>4. How many fold changes above 2 do we see?
A>
A>
A>
