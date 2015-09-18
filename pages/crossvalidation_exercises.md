---
Title:
---

## Exercises

Load the following dataset


```r
library(GSE5859Subset)
data(GSE5859Subset)
```

And define the outcome and predictors. To make the problem more difficult we will only consider autosomal genes:


```r
y = factor(sampleInfo$group)
X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out]
```

1. Use the `createFold` function in the `caret` package, set the seed to 1 `set.seed(1)` and create 10 folds. 

    Question: What is the 2nd entry in the fold 3


2. We are going to use kNN. We are going to consider a smaller set of predictors by using _filtering_ genes using t-tests. Specifically, we will perform a t-test and select the $$m$$ genes with the smallest p-values.

    Let $$m=8$$ and $$k=5$$ and train kNN by leaving out the second fold `idx[[2]]`. How many mistakes do we make on the test set? Remember it is indispensable that you perform the ttest on the training data.



3. Now run through all 5 folds. What is our error rate?



4. Now we are going to select the best values of $$k$$ and $$m$$. Use the expand grid function to try out the following values:

    
    ```r
    ms=2^c(1:11)
    ks=seq(1,9,2)
    params = expand.grid(k=ks,m=ms)
    ```

    Now use apply or a loop to obtain an error rates for each of these pairs of parameters. Which pair of parameters minimizes the error rate?


5. Repeat exercise 4 but now perform the t-test filtering before the cross validation. Note how this biases the entire result and gives us much lower estimated error rates.



6. Repeat exercise 3 but now instead of `sampleInfo$group` use 

    
    ```r
    y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))
    ```

    What is the minimum error rate now?


Note that we achieve much lower error rate when predicting date than when predicting the group. Because group is confounded with date, it is very possible that these predictors have no information about group and that our lower 0.5 error rates are due to the confounding with date. We will learn more about this in the batch effect chapter.








