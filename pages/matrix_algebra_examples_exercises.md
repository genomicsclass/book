---
Title: Matrix Algebra Examples Exercises
---

{pagebreak} 

A>## Exercises
A>
A>1. Suppose we are analyzing a set of 4 samples. The first two samples are from a treatment group A and the second two samples are from a treatment group B. This design can be represented with a model matrix like so:
A>
A>    
    ```r
    X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
    rownames(X) <- c("a","a","b","b")
    X
    ```
A>    
A>    ```
A>    ##   [,1] [,2]
A>    ## a    1    0
A>    ## a    1    0
A>    ## b    1    1
A>    ## b    1    1
A>    ```
A>
A>    Suppose that the fitted parameters for a linear model give us:
A>
A>    
    ```r
    beta <- c(5, 2)
    ```
A>
A>    Use the matrix multiplication operator, `%*%`, in R to answer the following questions:
A>    
A>    What is the fitted value for the A samples? (The fitted Y values.)
A>
A>
A>
A>2. What is the fitted value for the B samples? (The fitted Y values.)
A>
A>
A>
A>3. Suppose now we are comparing two treatments B and C to a control group A, each with two samples. This design can be represented with a model matrix like so:
A>
A>    
    ```r
    X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
    rownames(X) <- c("a","a","b","b","c","c")
    ```
A>  
A>    which results in a matrix that looks like:
A>
A>    {$$} 
A>    \,
A>\begin{pmatrix}
A>a&1&0&0\\   
A>a&1&0&0\\
A>b&1&1&0\\
A>b&1&1&0\\
A>c&1&0&1\\
A>c&1&0&1
A>\end{pmatrix}
A>    {/$$}
A>   
A>    Suppose that the fitted values for the linear model are given by:
A>
A>    
    ```r
    beta <- c(10,3,-3)
    ```
A>
A>    What is the fitted value for the B samples?
A>
A>
A>
A>4. What is the fitted value for the C samples?
A>
A>
