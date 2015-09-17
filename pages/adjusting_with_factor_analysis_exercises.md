---
Title: Adjusting with factor analysis exercises
---

A>## Exercises
A>
A>In this section we will use the `sva` function in the `sva` package (available from Biocondcutor) and apply it to the following data:
A>
A>
```r
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
```
A>
A>1. In a previous section we estimated factors using PCA. But we noted that the first factor was correlated with our outcome of interest: 
A>
A>    
    ```r
    s <- svd(geneExpression-rowMeans(geneExpression))
    cor(sampleInfo$group,s$v[,1])
    ```
A>
A>    The `svafit` function estimates factors, but downweighting the genes that appear to correlate with the outcome of interest. It also tries to estimate the number of factors and returns the estimated factors like this:
A>
A>    
    ```r
    sex = sampleInfo$group
    mod = model.matrix(~sex)
    svafit = sva(geneExpression,mod)
    head(svafit$sv)
    ```
A>
A>    Note that the resulting estimated factors are not that different from  the PCs
A>    
A>    
    ```r
    for(i in 1:ncol(svafit$sv)){
      print( cor(s$v[,i],svafit$sv[,i]) )
      }
    ```
A>
A>
A>    Now fit a linear model to each gene that instead of `month` includes these factors in the model. Use the `qvalue` function. 
A>    
A>    How many genes have q-value < 0.1?
A>
A>
A>
A>2. How many of these genes are from chrY or chrX?
A>
A>
A>
