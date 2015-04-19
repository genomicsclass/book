---
layout: page
title: Inference with bioc
---




# Introduction

In this section we will cover inference in the context of genomics experiments. We apply some of the concepts we have covered in previous sections including t-tests, multiple comparisons and standard deviation estimates from hierarchical models. 

We start by loading the pooling experiment data 



```r
library(Biobase)
library(maPooling)
data(maPooling)
individuals=which(rowSums(pd)==1)
```

```
## Error in is.data.frame(x): object 'pd' not found
```

And extracting the individual mice as well as their strain


```r
individuals=which(rowSums(pd)==1)
```

```
## Error in is.data.frame(x): object 'pd' not found
```

```r
individuals=individuals[-grep("tr",names(individuals))]
```

```
## Error in eval(expr, envir, enclos): object 'individuals' not found
```

```r
y=exprs(maPooling)[,individuals]
```

```
## Error in eval(expr, envir, enclos): object 'individuals' not found
```

```r
g=factor(as.numeric(grepl("b",names(individuals))))
```

```
## Error in grepl("b", names(individuals)): object 'individuals' not found
```


# T-test

We can now apply a t-test to each gene using the `rowttest` function in the `genefilter` package


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
tt=rowttests(y,g)
```

```
## Error in rowttests(y, g): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'y' not found
```

Now which genes do we report as statistically significant? For somewhat arbitrary reasons, in science p-values of 0.01 and 0.05 are used as cutoff. In this particular example we get 


```r
sum(tt$p.value<0.01)
```

```
## Error in eval(expr, envir, enclos): object 'tt' not found
```

```r
sum(tt$p.value<0.05)
```

```
## Error in eval(expr, envir, enclos): object 'tt' not found
```


# Multiple testing
We described multiple testing in detail in course 3. Here we provide a quick summary.

Do we report all these genes? Let's explore what happens if we split the first group into two, forcing the null hypothesis to be true


```r
set.seed(0)
shuffledIndex <- factor(sample(c(0,1),sum(g==0),replace=TRUE ))
```

```
## Error in sample.int(length(x), size, replace, prob): object 'g' not found
```

```r
nulltt <- rowttests(y[,g==0],shuffledIndex)
```

```
## Error in rowttests(y[, g == 0], shuffledIndex): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'y' not found
```

```r
sum(nulltt$p.value<0.01)
```

```
## Error in eval(expr, envir, enclos): object 'nulltt' not found
```

```r
sum(nulltt$p.value<0.05)
```

```
## Error in eval(expr, envir, enclos): object 'nulltt' not found
```

If we use the 0.05 cutoff we will be reporting 840 false positives. We have described several ways to adjust for this include the `qvalue` method available in the `qvalue` package. After this adjustment we include a smaller list of genes.


```r
library(qvalues)
```

```
## Error in library(qvalues): there is no package called 'qvalues'
```

```r
qvals = qvalue(tt$p.value)$qvalue
```

```
## Error in eval(expr, envir, enclos): could not find function "qvalue"
```

```r
sum(qvals<0.05)
```

```
## Error in eval(expr, envir, enclos): object 'qvals' not found
```

```r
sum(qvals<0.01)
```

```
## Error in eval(expr, envir, enclos): object 'qvals' not found
```

And now the null case generates fewer false positives:


```r
library(qvalues)
```

```
## Error in library(qvalues): there is no package called 'qvalues'
```

```r
nullqvals = qvalue(nulltt$p.value)$qvalue
```

```
## Error in eval(expr, envir, enclos): could not find function "qvalue"
```

```r
sum(nullqvals<0.05)
```

```
## Error in eval(expr, envir, enclos): object 'nullqvals' not found
```

```r
sum(nullqvals<0.01)
```

```
## Error in eval(expr, envir, enclos): object 'nullqvals' not found
```

