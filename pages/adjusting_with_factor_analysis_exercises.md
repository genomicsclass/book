---
layout: page
Title: Adjusting with factor analysis exercises
---

## Exercises

In this section we will use the `sva` function in the `sva` package (available from Bioconductor) and apply it to the following data:


```r
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
```

1. In a previous section we estimated factors using PCA, but we noted that the first factor was correlated with our outcome of interest: 

    
    ```r
    s <- svd(geneExpression-rowMeans(geneExpression))
    cor(sampleInfo$group,s$v[,1])
    ```

    The `svafit` function estimates factors, but downweighs the genes that appear to correlate with the outcome of interest. It also tries to estimate the number of factors and returns the estimated factors like this:

    
    ```r
    sex = sampleInfo$group
    mod = model.matrix(~sex)
    svafit = sva(geneExpression,mod)
    head(svafit$sv)
    ```

    The resulting estimated factors are not that different from the PCs.
    
    
    ```r
    for(i in 1:ncol(svafit$sv)){
      print( cor(s$v[,i],svafit$sv[,i]) )
      }
    ```


    Now fit a linear model to each gene that instead of `month` includes these factors in the model. Use the `qvalue` function. 
    
    How many genes have q-value < 0.1?



2. How many of these genes are from chrY or chrX?



