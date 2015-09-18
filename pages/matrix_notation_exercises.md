---
Title: Matrix Notation Exercesis
---


## Exercises

1. In R we have vectors and matrices. You can create your own vectors with the function c.

    
    ```r
    c(1,5,3,4)
    ```

    They are also the output of many functions such as

    
    ```r
    rnorm(10)
    ```

    You can turn vectors into matrices using functions such as rbind, cbind or matrix.

    Create the matrix from the vector 1:1000 like this:

    
    ```r
    X = matrix(1:1000,100,10)
    ```

    What is the entry in row 25, column 3 ?



2. Using the function cbind, create a 10 x 5 matrix with first column `x=1:10`. Then  add `2*x`, `3*x`, `4*x` and `5*x` tp columns 2 through 5. What is the sum of the elements of the 7th row?



3. Which of the following creates a matrix with multiples of 3 in the third column?
  - A) `matrix(1:60,20,3)`
  - B) `matrix(1:60,20,3,byrow=TRUE)`
  - C) `x=11:20; rbind(x,2*x,3*x)`
  - D) `x=1:40; matrix(3*x,20,2)`




