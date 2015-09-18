---
Title: Introduction to high throughput data exercises
---

## Exercises

For the remaining parts of this book we will be downloading larger datasets than those we have been using. Most of thes datasets are not avaialbe as part of the standard R installation or packages such as `UsingR`. For some of these packages we have created pacakges and offer them via GitHub. To download these you will need to install the `devtools` package. Once you do this you can install packages such as the `GSE5859Subset` which we will be using here:


```r
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)
```

This package loads three tables:  `geneAnnotation`, `geneExpression`, and `sampleInfo`. Answer the following questions to familiarize yourself with the data set:


1. How many samples where processed on 2005-06-27?



2. Question: How many of the genes represented in this particular technology are on chromosome Y? 



3.  What is the log expression value of the for gene ARPC1A
on the one subject that we measured on 2005-06-10 ?


