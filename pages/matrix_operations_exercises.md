---
Title: Matrix Operation Exercises
---


A>## Exercises
A>
A>1. Suppose `X` is a matrix in R. Which of the following is **not** equivalent to `X`?
A>  - A) `t( t(X) )  `
A>  - B) `X %*% matrix(1,ncol(X) ) `
A>  - C) `X*1`
A>  - D) `X%*%diag(ncol(X))`
A>  
A>
A>
A>2. Solve the following system of equations using R:
A>
A>    
A>    {$$}\begin{align*}
A>3a + 4b - 5c + d &= 10\\
A>2a + 2b + 2c - d &= 5\\
A>a -b + 5c - 5d &= 7\\
A>5a + d &= 4
A>\end{align*}{/$$}
A>    
A>    What is the solution for {$$}c{/$$}?
A>
A>
A>3. Load the following two matrices into R:
A>
A>    
    ```r
    a <- matrix(1:12, nrow=4)
    b <- matrix(1:15, nrow=3)
    ```
A>
A>    Note the dimension of 'a' and the dimension of 'b'.
A>
A>    In the question below, we will use the matrix multiplication operator in R, `%*%`, to multiply these two matrices.
A>    
A>    What is the value in the 3rd row and the 2nd column of the matrix product of 'a' and 'b'?
A>
A>
A>
A>4. Multiply the 3rd row of 'a' with the 2nd column of 'b', using the element-wise vector multiplication with `*`.
A>
A>    What is the sum of the elements in the resulting vector?
A>
