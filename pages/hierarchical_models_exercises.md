---
Title: Hierarchical Models Exercises
---

A>## Exercises
A>
A>
A>Load the following data (you can isntall it from Biocodcuctor) and extract the data matrix using `exprs`:
A>
A>
```r
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)
```
A>
A>This dataset comes from an experiment in which RNA was obtained from the same background pool to creat six replicate samples. Then RNA from 16 genes were artificially added in different quantities to each sample. These quantities (in picoMolars) and gene IDs are stored here:
A>
A>
```r
pData(rma95)
```
A>
A>
A>Note that these quantities where the same in the first three arrays and in the last three arrays. So we define two groups like this:
A>
A>
```r
g <- factor(rep(0:1,each=3))
```
A>
A>and create an index of which rows are associated with the artificially added genes:
A>
A>
```r
spike <- rownames(y) %in% colnames(pData(rma95))
```
A>
A>1. Note that only these 16 genes are differnetially expressed since these the six samples differ only due sampling (they all come from the same background pool of RNA). 
A>
A>    Perform a t-test on each gene using the `rowttest` function. 
A>
A>    What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?
A>
A>
A>2. Now compute the within group sample  standard deviation for each gene (you can use group 1). Based on the p-value cut-off split the genes into true positives, false positives, true negatives and false negatives. Create a boxplot comparing the sample SDs for each group. Which of the following best described the box-plot? 
A>    - A) The standard deviation is similar across all groups
A>    - B) On average, the true negatives have much larger variability
A>    - C) The false negatives have larger variability
A>    - D) The false positives have smaller standard devition
A>
A>
A>
A>
A>3. In the previous two questions we observed results consistent with the fact that the random variability associated with the sample standard deviation leads to t-statistics that are large by chance.
A>
A>    Note that the sample standard deviation we use in the t-test is an estimate and that with just a pair of triplicate samples the variability associated with the denominator in the t-test can be large.
A>
A>    The following three steps perform the basic `limma` analysis. We specify `coef=2` because we are interested in the difference between groups, not the intercept. The `eBayes` step uses a hierarchichal model that provides a new estimate of the gene specific standard error.
A>
A>    
    ```r
    library(limma)
    fit <- lmFit(y, design=model.matrix(~ g))
    colnames(coef(fit))
    fit <- eBayes(fit)
    ```
A>
A>    Here is a plot of the original new hierarchical models based estimate versus the sample based estimate 
A>
A>    
    ```r
    sampleSD = fit$sigma
    posteriorSD = sqrt(fit$s2.post)
    ```
A>
A>
A>    Which best describes what the hierarchichal model does
A>    - A) Moves all the estimatates of standard devtion closer to 0.12
A>    - B) Increases the estimates of standard deviation to increase t
A>    - C) Decreases the estimate of standard deviation
A>    - D) Decreases the effect size estimates
A>
A>
A>
A>
A>4. Use these new estimates of standard deviation in the denominartor of the t-test and compute p-valus. You can do it like this:
A>
A>    
    ```r
    library(limma)
    fit = lmFit(y, design=model.matrix(~ g))
    fit = eBayes(fit)
    ##second coefficient relates to diffences between group
    pvals = fit$p.value[,2] 
    ```
A>
A>    What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?
A>
A>
A>Compare to the previous volcano plot and notice that we no longer have small p-values for genes with small effect sizes. 
