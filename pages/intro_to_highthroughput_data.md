---
title: "Introduction to Advanced Statistics for the Life Sciences"
author: "Rafa"
date: "January 31, 2015"
output: html_document
layout: page
---



# Inference For High Dimensional Data

## Introduction

The R markdown document for this section is available [here](https://github.com/genomicsclass/labs/tree/master/advinference/intro_to_highthroughput_data.Rmd).

High-throughput technologies have changed basic biology and the biomedical sciences from data poor disciplines to data intensive ones. A specific example comes from research fields interested in understanding gene expression. Gene expression is the process in which DNA, the blueprint for life, is copied into RNA, the templates for the synthesis of proteins, the building blocks for life. In the 1990s, the analysis of gene expression data amounted to spotting black dots on a piece of paper or extracting a few numbers from standard curves. With high-throughput technologies, such as microarrays, this suddenly changed to sifting through tens of thousands of numbers. More recently, RNA-Sequencing has further increased data complexity. Biologists went from using their eyes or simple summaries to categorize results to having thousands (and now millions) of measurements per sample to analyze. In this chapter we will focus on statistical inference. Specifically, we focus on the problem of detecting differences in groups using statistical tests and quantifying uncertainty in a meaningful way. In later chapters we will study the statistics behind clustering, machine learning, factor analysis and multi-level modeling. 

Since there is a vast number of available public datasets, we use several gene expression examples. Nonetheless, the statistical techniques you will learn have also proven useful in other fields that make use of high-throughput technologies. Technologies such as microarrays, next generation sequencing, fRMI, and mass spectrometry all produce data to answer questions for which what we learn here will be indispensable. The specific topics we will focus on are inference in the context of high-throughput data, distance and clustering, dimension reduction, machine learning, modeling including Bayesian/hierarchical models and advanced exploratory data analysis. Because there is an interplay between these topics, we will cover each separately. 


```r
library(rafalib)
```

<a name="threetables"></a>

#### Data packages
Several of the  examples we are going to use in the following sections are best obtained through R packages. These are available from GitHub and can be installed using the `install_github` function from the `devtools` package. Microsoft Windows users might need to follow [these instructions](https://github.com/genomicsclass/windows) to properly install `devtools`. 

Once `devtools` is installed, you can then install the data packages like this:


```r
library(devtools)
install_github("genomicsclass/GSE5859Subset")
```

#### The three tables

Most of the data we use as examples in this class are created with high-throughput technologies. These technologies measure thousands of _features_. Examples of feature are genes, single base locations of the genome, genomic regions, or image pixel intensities. Each specific measurement product is defined by a specific set of features. For example, a specific gene expression microarray product is defined by the set of genes that it measures. 

A specific study will typically use one product to make measurements on several experimental units, such as individuals. The most common experimental unit will be the individual, but they can also be defined by other entities, for example different parts of a tumor. We often call the experimental units _samples_ following because experimental jargon. It is important that these are not confused with samples as referred to in previous chapters, for example "random sample". 

So a high-throughput experiment is usually defined by three tables: one with the high-throughput measurements and two tables with information about the columns and rows of this first table respectively.

Because a dataset is typically defined by a set of experimental units and a product defines a fixed set of features, the high-throughput measurements can be stored in an {$$}n \times m{/$$} matrix with {$$}n{/$$} the number of units and {$$}m{/$$} the number of features. In R the convention has been to store the transpose of these matrices. 

Here is an example from a gene expression dataset:


```r
##can be installed with:
library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables
dim(geneExpression)
```

```
## [1] 8793   24
```

We have RNA expression measurements for 8793 genes from blood taken from 24 individuals (the experimental units). For most statistical analysis we will also need information about the individuals. For example, in this case the data was originally collected to compare gene expression across ethnic groups. However, we have created a subset of this dataset for illustration and separated the data into two groups:



```r
dim(sampleInfo)
```

```
## [1] 24  4
```

```r
head(sampleInfo)
```

```
##     ethnicity       date         filename group
## 107       ASN 2005-06-23 GSM136508.CEL.gz     1
## 122       ASN 2005-06-27 GSM136530.CEL.gz     1
## 113       ASN 2005-06-27 GSM136517.CEL.gz     1
## 163       ASN 2005-10-28 GSM136576.CEL.gz     1
## 153       ASN 2005-10-07 GSM136566.CEL.gz     1
## 161       ASN 2005-10-07 GSM136574.CEL.gz     1
```

```r
sampleInfo$group
```

```
##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0
```

One of the columns, filenames, permits us to connect the rows of this table to the columns of the measurement table.


```r
match(sampleInfo$filename,colnames(geneExpression))
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24
```


Finally, we have a table describing the features:


```r
dim(geneAnnotation)
```

```
## [1] 8793    4
```

```r
head(geneAnnotation)
```

```
##      PROBEID  CHR     CHRLOC SYMBOL
## 1  1007_s_at chr6   30852327   DDR1
## 30   1053_at chr7  -73645832   RFC2
## 31    117_at chr1  161494036  HSPA6
## 32    121_at chr2 -113973574   PAX8
## 33 1255_g_at chr6   42123144 GUCA1A
## 34   1294_at chr3  -49842638   UBA7
```

The table includes an ID that permits us to connect the rows of this table with the rows of the measurement table:

```r
head(match(geneAnnotation$PROBEID,rownames(geneExpression)))
```

```
## [1] 1 2 3 4 5 6
```
The table also includes biological information about the features; namely chromosome location and the gene "name" used by biologists.

#### Examples

Here we list some of the examples of data analysis questions we might be asked to answer with the dataset shown here:

* Inference: For which genes are the population averages different across ethnic groups? 

* Machine learning: Build an algorithm that, given gene expression patterns, predicts ethnic group.

* Clustering: Can we discover subpopulations of individuals from the gene expression patterns? Or can we discover genes pathways based on which cluster together across individuals?

* Exploratory data analysis: Did some experiments failed experiments? Are the assumptions needed to use standard statistical techniques met? 

We will cover all these topics and more in the following sections. 

