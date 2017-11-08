---
layout: page
title: "Inference: t-tests, multiple comparisons"
---




# Introduction

In the previous section, we focused on a pair of genes to
illustrate two aspects of variation.  One of the genes appeared to
have high between-mouse variation that was hidden in the act
of pooling samples within strain.  When strains were compared on
the basis of the pooled data, there was an appearance of a significant
strain
effect for this gene ($$p < 10^{-6}$$), but when individual-level data were used to
perform the comparison, the strain effect was found to be very
weak at best ($$p = 0.089$$).  The lesson is to recognize that the
most scientifically compelling questions concern biological variation,
which can only be directly measured with good experimental design.  Accurate
interpretation of origin and size of biological variation requires
appropriate statistical analysis.

In this section we will cover inference in the context of genome-scale experiments.  There are several serious conceptual problems:

- there are many tests, often at least one test for each one of tens of thousands of features
- each feature (typically a gene) exhibits its own technical and biological variability
- there may be unmeasured or unreported sources of biological variation (such as time of day)
- many features are inherently interrelated, so the tests are not independent

We will apply some of the concepts we have covered in previous 
sections including t-tests and multiple comparisons; later we will
compute standard deviation estimates from hierarchical models. 

We start by loading the pooling experiment data 



```r
library(Biobase)
library(maPooling)
data(maPooling)
pd=pData(maPooling)
individuals=which(rowSums(pd)==1)
```

And extracting the individual mice as well as their strain


```r
individuals=which(rowSums(pd)==1)
individuals=individuals[-grep("tr",names(individuals))]
y=exprs(maPooling)[,individuals]
g=factor(as.numeric(grepl("b",names(individuals))))
```

<a name="rowWiseT"></a>

# T-tests

We can now apply a t-test to each gene using the `rowttest` function in the `genefilter` package


```r
library(genefilter)
tt=rowttests(y,g)
```

<a name="naive"></a>

Now which genes do we report as statistically significant? For somewhat arbitrary reasons, in science p-values of 0.01 and 0.05 are used as cutoff. In this particular example we get 


```r
NsigAt01 = sum(tt$p.value<0.01)
NsigAt01
```

```
## [1] 1578
```

```r
NsigAt05 = sum(tt$p.value<0.05)
NsigAt05
```

```
## [1] 2908
```

<a name="sham"></a>

# Multiple testing
We described multiple testing in detail [in course 3](http://genomicsclass.github.io/book/pages/multiple_testing.html). Here we provide a quick summary.

Do we report all the nominally significant
genes identified above? Let's explore what happens if we split the first group into two, forcing the null hypothesis to be true


```r
set.seed(0)
shuffledIndex <- factor(sample(c(0,1),sum(g==0),replace=TRUE ))
nulltt <- rowttests(y[,g==0],shuffledIndex)
NfalselySigAt01 = sum(nulltt$p.value<0.01)
NfalselySigAt01 
```

```
## [1] 79
```

```r
NfalselySigAt05 = sum(nulltt$p.value<0.05)
NfalselySigAt05
```

```
## [1] 840
```





