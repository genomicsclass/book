---
Title: Factor Analsysis Exercises
---

## Exercises

We will continue to use this dataset:

```r
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
```

1. Suppose you want to make an MA plot of the first two samples `y = geneExpression[,1:2]`. Which of the following projections gives us the projection of $$y$$ so that column2 versus column 1 is an MA plot?

    $$
    \mbox{A.} \, \, \,
    y\begin{pmatrix}
    1/\sqrt{2}&1/\sqrt{2} \\ 
    1\sqrt{2}&-1/\sqrt{2}
    \end{pmatrix} \,  \, \,
    \mbox{B.} \, \, \, y\begin{pmatrix} 
    1&1 \\ 
    1&-1
    \end{pmatrix}  \, \, \,
    \mbox{C.} \, \, \, 
    \begin{pmatrix} 
    1&1 \\ 
    1&-1
    \end{pmatrix} y \, \, \,
    \mbox{D.} \, \, \,
    \begin{pmatrix} 
    1&1 \\ 
    1&-1
    \end{pmatrix} y^\top
    $$



2. Say $$Y$$ is $$M \times N$$, in the SVD $$Y=UDV^\top$$ which of the following is not correct?

    - A) $$DV^\top$$ are the new coordinates for the projection $$U^\top Y$$
    - B) $$UD$$ are the new coordinates for the projection $$YV$$
    - C) $$D$$ are the coordinates of the projection $$U^\top Y$$
    - D) $$U^\top Y$$ is a projection from an $$N$$-dimensional to $$M$$-dimensional subspace.  



3. Define

    
    ```r
    y = geneExpression - rowMeans(geneExpression)
    ```

    Compute and plot an image of the correlation for each sample. Make two image plots of these correlations. In the first one, plot the correlation as image. In the second, order the samples by date and then plot the an image of the correlation. The only difference in these plots is the order in which the samples are plotted.

    Based on these plots, which of the following you would say is true:
    - A) The samples appear to be completely independent of each other
    - B) Sex seems to be creating structures as evidenced by the two cluster of highly correlated samples
    - C) There appear to be only two factors completely driven by month
    - D) The fact the in the plot ordered by month we see two groups mainly driven by month and within these, we subgroups driven by date seems to suggest date more than month per se are the hidden factors



4. Based on the correlation plots above, we could argue that there are at least two hidden factors. Using PCA estimate these two factors. Specifically, apply the `svd` to `y` and use the first two PCs as estimates.

    Which command gives us these estimates?

    - A) `pcs = svd(y)$v[1:2,]`
    - B) `pcs = svd(y)$v[,1:2]`
    - C) `pcs = svd(y)$u[,1:2]`
    - D) `pcs = svd(y)$d[1:2]`



5. Plot each of the estimated factor ordered by date. Use color to denote month. The first factor is clearly related to date. 
Which of the following appear to be most different according to this factor?
    - A) June 23 and June 27
    - B) Oct 07 and Oct 28
    - C) June 10 and June 23
    - D) June 15 and June 24



6. Use the `svd` function to obtain the principal components (PCs) for our detrended gene expression data `y`:

    How many PCs explain more than 10% of the variability?




7. Which PC most correlates (negative or positive correlation) with month? 


8. What is this correlation  (in absolute value)?



9. Which PC most correlates (negative or positive correlation) with sex? 


10. What is this correlation  (in absolute value)?


11. Now instead of using month, which we have shown does not quite describe the batch,  add the two estimated factors in exercise 9 to the linear model we used above.

    
    ```r
    X <- model.matrix(~sex+s$v[,1:2])
    ```

    And apply this model to each gene, exact p-values for the sex difference.  How many q-values <0.1 for the sex comparison?
	


12. What proportion of the genes are on chrX and chrY?


