---
Title: EDA with PCA exercises
---

A>## Exercises
A>
A>We will use the Biocondcutro pacakge `Biobase` which you can install with `install_bioc` function from `rafalib`:
A>
A>
A>Load the data for this gene expression dataset:
A>
A>
```r
library(Biobase)
library(GSE5859)
data(GSE5859)
```
A>
A>Note that this is the original dataset from which we selected the subset used in `GSE5859Subset`. 
A>
A>We can extract the gene expression data and sample information table using the Bioconductor functions `exprs` and `pData` like this:
A>
A>
```r
geneExpression = exprs(e)
sampleInfo = pData(e)
```
A>
A>1. Familiarize yourself with the `sampleInfo` table. Note that some samples were processed at different times. This is an extraneous variable and should not affect the values in `geneExpression`. However, as we have seen in previous analyses it does appear to have an effect so we will explore this here.
A>
A>    You can extract the year from each date like this:
A>    
A>    
    ```r
    year = format(sampleInfo$date,"%y")
    ```
A>
A>    Note that etnic group and year is almost perfectly co
A>
A>    
    ```r
    length( unique(year) )
    ```
A>  
A>    unique years for which we have data.
A>
A>1. For how many of these years do we have more than one ethnicity represented?
A>
A>
A>2. Repeat the above exercise but now instead of year consider the month as well. Specifically, instead of the `year` variable defined above use:
A>
A>    
    ```r
    month.year = format(sampleInfo$date,"%m%y")
    ```
A>
A>    For what **proportion** of these `month.year` values do we have more than one ethnicity represented?
A>
A>
A>
A>3. Perform a t-test (use `rowttests`) comparing CEU samples processed in 2002 to those processed in 2003. Then use the `qvalue` package to obtain q-values for each gene. 
A>
A>    How many genes have q-values < 0.05 ?
A>
A>
A>4. What is the estimate of `pi0` provided by `qvalue`: 
A>
A>
A>5. Now perform a t-test (use `rowttests`) comparing CEU samples processed in 2003 to those processed in 2004. Then use the `qvalue` package to obtain q-values for each gene. How many genes have q-values less than 0.05?
A>
A>
A>
A>
A>6. Now we are going to compare ethnicities as was done in the original publication in which these data were first presented. Use the `qvalue` function to compare the ASN population to the CEU population. Once again, use the `qvalue` function to obtain q-values.
A>
A>    How many genes have q-values < 0.05 ?
A>
A>
A>
A>
A>7. Note that over 80% of genes are called differentially expressed between ethnic groups. However, due to the confounding with processing date, we need to confirm these differences are actually due to ethnicity. This will not be easy due to the almost perfect confounding. However, above we noted that two groups were represented in 2005. Just like we stratified by majors to remove the "major effect" in our admissions example, here we can stratify by year and perform a t-test comparing ASN and CEU, but only for samples processed in 2005.
A>
A>    How many genes have q-values < 0.05 ?
A>
A>
A>
A>    Note the dramatic drop in the number of gchenes with q-value < 0.05 when we fix the year. However, the sample size is much smaller in this latest analysis which means we have less power:
A>
A>    
    ```r
    table(sampleInfo$ethnicity[index])
    ```
A>
A>
A>8. To provide a more balanced comparison we repeat the analysis but now taking 3 random CEU samples from 2002. Repeat the analysis above but comparing the ASN from 2005 to three random CEU samples from 2002. Set the seed at 3, `set.seed(3)`
A>
A>    How many genes have q-values < 0.05 ?
A>  
A>
A>
