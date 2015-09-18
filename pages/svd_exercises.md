---
Title: SVD exercises
---

{pagebreak} 

## Exercises

For these exercises we are again going to use:


```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
```

Before we start these exercises it is important to reemphasize that when using the SVD in practice the solution to SVD is not unique. This because $$\mathbf{UDV^\top} = \mathbf{ (-U) D (-V)^\top}$$. In fact we can flip the sign of each column of $$\mathbf{U}$$ and as long as we also flip the respective column in $$\mathbf{V}$$ we will get solution. Here is an example:


```r
s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips
```

Now we switch the sign of each column and check that we get the same answer. We do this using the function `sweep`. If `x` is a matrix and `a` is a vector than `sweep(x,1,y,FUN="*")` applies the fun `FUN` to each row `i` FUN(x[i,],a[i])`, in this case `x[i,]*a[i]`. If instead of 1 we use 2 sweep applies this to columns. To learn about sweep read `?sweep`. 


```r
newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
identical( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))
```

This is important to know because different implementations of the SVD algorithm may give different signs which can result in the same code resulting in different answers when run in different computer systems.

1. Compute the SVD of e 

    
    ```r
    s = svd(e)
    ```

    Now compute the mean of each row:

    
    ```r
    m = rowMeans(e)
    ```

    What is the correlation between the first column of $$\mathbf{U}$$ and `m`?


2. In exercise 1 we saw how the first column relates to the mean of the rows of `e`. Note that if we change these means the distances between columns do not change. Note for example that changing the means does not change the distances:

    
    ```r
    newmeans = rnorm(nrow(e)) ##random values we will add to create new means
    newe = e+newmeans ##we change the means
    sqrt(crossprod(e[,3]-e[,45]))
    sqrt(crossprod(newe[,3]-newe[,45])) 
    ```

    So we might as well make the mean of each row 0 since it does not help us approximate the column distances. We will define `y` as the _detrended_ `e` and recompute the SVD:

    
    ```r
    y = e - rowMeans(e)
    s = svd(y)
    ```

    We showed that $$\mathbf{UDV^\top}$$ is equal to `y` up to numerical error
  
    
    ```r
    resid = y - s$u %*% diag(s$d) %*% t(s$v)
    max(abs(resid))
    ```

    The above can be made more efficient in two ways. First, using the `crossprod` and second not creating a diagonal matrix. Note that in R we can multiply matrices `x` by vector `a`. The result is a matrix with rows `i` equal to `x[i,]*a[i]`. Run this example to see this.

    
    ```r
    x=matrix(rep(c(1,2),each=5),5,2)
    x
    x*c(1:5)
    ```

    which is equivalent to

    
    ```r
    sweep(x,1,1:5,"*")
    ```

    This means that we don't have to convert `s$d` into a matrix. 

    Which of the following gives us the same as `diag(s$d)%*%t(s$v)`
    - A) `s$d %*% t(s$v)`
    - B) `s$d * t(s$v)`
    - C) `t(s$d * s$v)`
    - D) `s$v * s$d`




3. If we define `vd = t(s$d * t(s$v))` then which of the following is not the same  $$UDV^\top$$:
    - A) `tcrossprod(s$u,vd)`
    - B) `s$u %*% s$d * t(s$v)`
    - C) `s$u %*% (s$d * t(s$v) )`
    - D) `tcrossprod( t( s$d*t(s$u)) , s$v)`




4. Let `z = s$d * t(s$v)`. We showed derivation demonstrating that because $$\mathbf{U}$$ is orthogonal the distance between `e[,3]` and `e[,45]` is the same as the distance between `y[,3]` and `y[,45]` which is the same as `vd[,3]` and `vd[,45]`


    
    ```r
    z = s$d * t(s$v)
    ## d was deinfed in question 2.1.5
    sqrt(crossprod(e[,3]-e[,45]))
    sqrt(crossprod(y[,3]-y[,45]))
    sqrt(crossprod(z[,3]-z[,45]))
    ```

    Note that the columns `z` have 189 entries, compared to 22,215 for `e`. 

    What is the difference (in absolute value) between the actual distance `sqrt(crossprod(e[,3]-e[,45]))` and the approximation using only two dimension of `z`





5. How many dimensions do we need to use for the approximation in exercise 4 to be within 10% 




6. Compute distances between sample 3 and all other samples:


7. Recompute this distance using the 2 dimensional approximation. 
    What is the Spearman correlation between this approximate distance and the actual distance?



Note that the last exercise shows how just two dimensions can be useful to get a rough idea about the actual distances.

