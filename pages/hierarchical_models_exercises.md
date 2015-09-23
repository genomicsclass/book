---
layout: page
title: Hierarchical Models Exercises
---

## Exercises


Load the following data (you can install it from Bioconductor) and extract the data matrix using `exprs`:


```r
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)
```

This dataset comes from an experiment in which RNA was obtained from the same background pool to create six replicate samples. Then RNA from 16 genes were artificially added in different quantities to each sample. These quantities (in picoMolars) and gene IDs are stored here:


```r
pData(rma95)
```


These quantities where the same in the first three arrays and in the last three arrays. So we define two groups like this:


```r
g <- factor(rep(0:1,each=3))
```

and create an index of which rows are associated with the artificially added genes:


```r
spike <- rownames(y) %in% colnames(pData(rma95))
```

1. Only these 16 genes are diferentially expressed since the six samples differ only due to sampling (they all come from the same background pool of RNA). 

    Perform a t-test on each gene using the `rowttest` function. 

    What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?


2. Now compute the within group sample standard deviation for each gene (you can use group 1). Based on the p-value cut-off, split the genes into true positives, false positives, true negatives and false negatives. Create a boxplot comparing the sample SDs for each group. Which of the following best describes the boxplot? 
    - A) The standard deviation is similar across all groups.
    - B) On average, the true negatives have much larger variability.
    - C) The false negatives have larger variability.
    - D) The false positives have smaller standard deviation.




3. In the previous two questions, we observed results consistent with the fact that the random variability associated with the sample standard deviation leads to t-statistics that are large by chance.

    The sample standard deviation we use in the t-test is an estimate and with just a pair of triplicate samples, the variability associated with the denominator in the t-test can be large.

    The following steps perform the basic `limma` analysis. We specify `coef=2` because we are interested in the difference between groups, not the intercept. The `eBayes` step uses a hierarchical model that provides a new estimate of the gene specific standard error.

    
    ```r
    library(limma)
    fit <- lmFit(y, design=model.matrix(~ g))
    colnames(coef(fit))
    fit <- eBayes(fit)
    ```

    Here is a plot of the original, new, hierarchical models based estimate versus the sample based estimate:

    
    ```r
    sampleSD = fit$sigma
    posteriorSD = sqrt(fit$s2.post)
    ```


    Which best describes what the hierarchical model does?
    
    - A) Moves all the estimates of standard deviation closer to 0.12.
    - B) Increases the estimates of standard deviation to increase t.
    - C) Decreases the estimate of standard deviation.
    - D) Decreases the effect size estimates.




4. Use these new estimates of standard deviation in the denominator of the t-test and compute p-values. You can do it like this:

    
    ```r
    library(limma)
    fit = lmFit(y, design=model.matrix(~ g))
    fit = eBayes(fit)
    ##second coefficient relates to diffences between group
    pvals = fit$p.value[,2] 
    ```

    What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?


Compare to the previous volcano plot and notice that we no longer have small p-values for genes with small effect sizes. 
