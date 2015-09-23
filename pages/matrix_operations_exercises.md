---
layout: page
Title: Matrix Operation Exercises
---


## Exercises

1. Suppose `X` is a matrix in R. Which of the following is **not** equivalent to `X`?
  - A) `t( t(X) )  `
  - B) `X %*% matrix(1,ncol(X) ) `
  - C) `X*1`
  - D) `X%*%diag(ncol(X))`
  


2. Solve the following system of equations using R:

    
    $$\begin{align*}
3a + 4b - 5c + d &= 10\\
2a + 2b + 2c - d &= 5\\
a -b + 5c - 5d &= 7\\
5a + d &= 4
\end{align*}$$
    
    What is the solution for $$c$?


3. Load the following two matrices into R:

    
    ```r
    a <- matrix(1:12, nrow=4)
    b <- matrix(1:15, nrow=3)
    ```

    Note the dimension of `a` and the dimension of `b`.

    In the question below, we will use the matrix multiplication operator in R, `%*%`, to multiply these two matrices.
    
    What is the value in the 3rd row and the 2nd column of the matrix product of `a` and `b`?



4. Multiply the 3rd row of `a` with the 2nd column of `b`, using the element-wise vector multiplication with `*`.

    What is the sum of the elements in the resulting vector?

