---
Title: Matrix Notation Exercesis
---


A>## Exercises
A>
A>1. In R we have vectors and matrices. You can create your own vectors with the function c.
A>
A>    
    ```r
    c(1,5,3,4)
    ```
A>
A>    They are also the output of many functions such as
A>
A>    
    ```r
    rnorm(10)
    ```
A>
A>    You can turn vectors into matrices using functions such as rbind, cbind or matrix.
A>
A>    Create the matrix from the vector 1:1000 like this:
A>
A>    
    ```r
    X = matrix(1:1000,100,10)
    ```
A>
A>    What is the entry in row 25, column 3 ?
A>
A>
A>
A>2. Using the function cbind, create a 10 x 5 matrix with first column
A>
A>    
    ```r
    x=1:10
    ```
A>    Then columns `2*x`, `3*x`, `4*x` and `5*x` in columns 2 through 5.
A>
A>    What is the sum of the elements of the 7th row?
A>
A>
A>
A>3. Which of the following creates a matrix with multiples of 3 in the third column?
A>  - A) `matrix(1:60,20,3)`
A>  - B) `matrix(1:60,20,3,byrow=TRUE)`
A>  - C) `x=11:20; rbind(x,2*x,3*x)`
A>  - D) `x=1:40; matrix(3*x,20,2)`
A>
A>
A>
A>
