---
layout: page
title: Cross-validation
---



In this lab, we will explore a method for picking parameters in a prediction / machine learning task, which is called *cross-validation*.

Suppose we have a prediction algorithm which is going to predict the class of some observations using a number of features. For example, we will use the gene expression values to predict the tissue type in our tissues gene expression dataset.

If this algorithm has a parameter which controls the behavior, we might pick the value of this parameter which minimizes the classification error. However, trying to classify the same observations as we use to *train* the model can be misleading.  In lecture, we saw that for K-nearest neighbors, using k=1 will always give 0 classification error in the training set (because we use the single observation to classify itself). Instead, it's better to pick the parameter using the algorithms performance on a set of observations which the algorithm has never seen, a *test* set.

Cross-validation is simply a method which splits the data into a number of *folds*. If we have N folds, then the algorithm typically trains on (N-1) of the folds, and test the algorithms performance on the left-out single fold. This is then repeated N times until each fold has been used as a *test* set.

Let's load in the tissue gene expression dataset:



```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
```

For illustration purposes 
let's drop one of the tissues which doesn't have many samples:


```r
table(tissue)
```

```
## tissue
##  cerebellum       colon endometrium hippocampus      kidney       liver 
##          38          34          15          31          39          26 
##    placenta 
##           6
```

```r
ind <- which(tissue != "placenta")
y <- tissue[ind]
X <- t( e[,ind] )
```

We will use the `createFolds` function from the `caret` package to make 5 folds of the data, which are balanced over the tissues. Don't be confused that the `createFolds` function uses the same letter 'k' as the k in K-nearest neighbors. These 'k' are unrelated.  The caret function `createFolds` is asking for how many folds to create, the 'N' from above. The `knn` function is asking how many closest observations to use to classify the test observations.


```r
# install.packages("caret")
#library(class)
library(caret)
```

```
## Error in library(caret): there is no package called 'caret'
```

```r
set.seed(1)
idx <- createFolds(y, k=10)
```

```
## Error in eval(expr, envir, enclos): could not find function "createFolds"
```

```r
sapply(idx, function(i) table(y[i]))
```

```
## Error in lapply(X = X, FUN = FUN, ...): object 'idx' not found
```


Because tissues have very different gene expression profiles, predicting tissue with all genes will be too easy. For illustration purposes we will try to predict tissue type with just two dimensional data. We will reduce dimension using `cmdscale`


```r
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
Xsmall <- cmdscale(dist(X))
plot(Xsmall,col=as.fumeric(y))
legend("topleft",levels(factor(y)),fill=seq_along(levels(factor(y))))
```

![plot of chunk unnamed-chunk-4](figure/crossvalidation-unnamed-chunk-4-1.png) 

Now we can try out the K-nearest neighbors method on a single fold:


```r
library(class)
i=1
pred <- knn(train =  Xsmall[ -idx[[i]] , ], test = Xsmall[ idx[[i]], ], cl=  y[ -idx[[i]] ], k=5)
```

```
## Error in as.matrix(train): object 'idx' not found
```

```r
table(true=y[ idx[[i]] ], pred)
```

```
## Error in table(true = y[idx[[i]]], pred): object 'idx' not found
```

```r
mean(y[ idx[[i]] ] != pred)
```

```
## Error in mean(y[idx[[i]]] != pred): object 'idx' not found
```

Now we will create a loop, which tries out each value of k from 1 to 12, and runs the K-nearest neighbors algorithm on each fold. We then ask for the proportion of errors for each fold, and report the average from the 5 cross-validation folds:


```r
set.seed(1)
ks <- 1:12
res <- sapply(ks, function(k) {
  # try out each version of k from 1 to 12
  
  res.k <- sapply(seq_along(idx), function(i) {
    # loop over each of the 5 cross-validation folds

    # predict the held-out samples using k nearest neighbors
    pred <- knn(train = Xsmall[ -idx[[i]], ],
                test = Xsmall[ idx[[i]], ],
                cl = y[ -idx[[i]] ], k = k)

    # the ratio of misclassified samples
    mean(y[ idx[[i]] ] != pred)
  })
  
  # average over the 5 folds
  mean(res.k)
})
```

```
## Error in lapply(X = X, FUN = FUN, ...): object 'idx' not found
```

Now we can plot the mean misclassification rate for each value of k:


```r
plot(ks, res, type="o")
```

```
## Error in xy.coords(x, y, xlabel, ylabel, log): object 'res' not found
```


Finally, to show that gene expression can perfectly predict tissue, we use 5 dimensions instead of 2 and note we get perfect prediction


```r
Xsmall <- cmdscale(dist(X),k=5)
set.seed(1)
ks <- 1:12
res <- sapply(ks, function(k) {
  res.k <- sapply(seq_along(idx), function(i) {
    pred <- knn(train = Xsmall[ -idx[[i]], ],
                test = Xsmall[ idx[[i]], ],
                cl = y[ -idx[[i]] ], k = k)
    mean(y[ idx[[i]] ] != pred)
  })
  mean(res.k)
})
```

```
## Error in lapply(X = X, FUN = FUN, ...): object 'idx' not found
```

```r
plot(ks, res, type="o",ylim=c(0,0.20))
```

```
## Error in xy.coords(x, y, xlabel, ylabel, log): object 'res' not found
```

Important note: We applied `cmdscale` to the entire dataset to create a smaller one for illustration purposes. However, in a real machine learning application all transformations of the data must be applied separately on the test and training dataset.
