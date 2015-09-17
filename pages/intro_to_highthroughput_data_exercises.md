---
Title: Introduction to high throughput data exercises
---

A>## Exercises
A>
A>For the remaining parts of this book we will be downloading larger datasets than those we have been using. Most of thes datasets are not avaialbe as part of the standard R installation or packages such as `UsingR`. For some of these packages we have created pacakges and offer them via GitHub. To download these you will need to install the `devtools` package. Once you do this you can install packages such as the `GSE5859Subset` which we will be using here:
A>
A>
```r
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)
```
A>
A>This package loads three tables:  `geneAnnotation`, `geneExpression`, and `sampleInfo`. Answer the following questions to familiarize yourself with the data set:
A>
A>
A>1. How many samples where processed on 2005-06-27?
A>
A>
A>
A>2. Question: How many of the genes represented in this particular technology are on chromosome Y? 
A>
A>
A>
A>3.  What is the log expression value of the for gene ARPC1A
A>on the one subject that we measured on 2005-06-10 ?
A>
A>
