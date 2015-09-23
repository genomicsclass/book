---
layout: page
title: Matrix Algebra Examples Exercises
---

{pagebreak} 

## Exercises

1. Suppose we are analyzing a set of 4 samples. The first two samples are from a treatment group A and the second two samples are from a treatment group B. This design can be represented with a model matrix like so:

    
    ```r
    X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
    rownames(X) <- c("a","a","b","b")
    X
    ```
    
    ```
    ##   [,1] [,2]
    ## a    1    0
    ## a    1    0
    ## b    1    1
    ## b    1    1
    ```

    Suppose that the fitted parameters for a linear model give us:

    
    ```r
    beta <- c(5, 2)
    ```

    Use the matrix multiplication operator, `%*%`, in R to answer the following questions:
    
    What is the fitted value for the A samples? (The fitted Y values.)



2. What is the fitted value for the B samples? (The fitted Y values.)



3. Suppose now we are comparing two treatments B and C to a control group A, each with two samples. This design can be represented with a model matrix like so:

    
    ```r
    X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
    rownames(X) <- c("a","a","b","b","c","c")
    X
    ```

    Suppose that the fitted values for the linear model are given by:

    
    ```r
    beta <- c(10,3,-3)
    ```

    What is the fitted value for the B samples?



4. What is the fitted value for the C samples?


