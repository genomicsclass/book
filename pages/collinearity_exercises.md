---
Title: Collinearity  exercises
---

{pagebreak} 

A>## Exercises
A>
A>
A>Consider these design matrices:
A>
A>{$$}
A>A=\begin{pmatrix}
A>1&0& 0& 0\\
A>1& 0& 0& 0\\
A>1& 1& 1& 0\\
A>1& 1& 0& 1
A>\end{pmatrix}\,\,
A>B=\begin{pmatrix}
A>1& 0& 0& 1\\
A>1& 0& 1& 1\\
A>1& 1& 0& 0\\
A>1& 1& 1& 0
A>\end{pmatrix}\,\,
A>C=\begin{pmatrix}
A>1& 0& 0\\
A>1& 1& 2\\
A>1& 2& 4\\
A>1& 3& 6
A>\end{pmatrix} 
A>{/$$}
A>
A>{$$}
A>D=\begin{pmatrix}
A>1& 0& 0& 0& 0\\
A>1& 0& 0& 0& 1\\
A>1& 1& 0& 1& 0\\
A>1& 1& 0& 1& 1\\
A>1& 0& 1& 1& 0\\
A>1& 0& 1& 1& 1
A>\end{pmatrix}\,\,
A>E=\begin{pmatrix}
A>1& 0& 0& 0\\
A>1& 0& 1& 0\\
A>1& 1& 0& 0\\
A>1& 1& 1& 1
A>\end{pmatrix} \,\,
A>F=\begin{pmatrix}
A>1& 0& 0& 1\\
A>1& 0& 0& 1\\
A>1& 0& 1& 1\\
A>1& 1& 0& 0\\
A>1& 1& 0& 0\\
A>1& 1& 1& 0
A>\end{pmatrix}
A>{/$$}
A>
A>
A>1. Which of the above design matrices does NOT have the problem of collinearity?
A>
A>
A>2. Let's use the example from the lecture to visualize how there is not a single best beta-hat, when the design matrix has collinearity of columns. An example can be made with:
A>
A>    
    ```r
    sex <- factor(rep(c("female","male"),each=4))
    trt <- factor(c("A","A","B","B","C","C","D","D"))
    ```
A>
A>    The model matrix can then be formed with:
A>
A>    
    ```r
    X <- model.matrix( ~ sex + trt)
    ```
A>
A>    And we can see that the number of independent columns is less than the number of columns of X:
A>
A>    
    ```r
    qr(X)$rank
    ```
A>    Suppose we observe some outcome, Y. For simplicity we will use synthetic data:
A>
A>    
    ```r
    Y <- 1:8
    ```
A>
A>    Now, we will fix the value for two beta's and optimize the remaining betas. We will fix {$$}\beta_{male}{/$$} and {$$}\beta_D{/$$}. And then we will find the optimal value for the remaining betas, in terms of minimizing {$$}\sum((Y - X \beta)^2){/$$}. The optimal value for the other entries in {$$}\beta{/$$}, we find those that minimizes:
A>
A>    {$$}
A>\sum_{i=1}^8  ( (Y_i - X_{i,male} \beta_{male} - X_{i,D} \beta_{i,D}) - Z_i \gamma )^2 )
A>    {/$$}
A>
A>    Where {$$}X_{male}{/$$} is the male column of the design matrix, {$$}X_D{/$$} is the D column, Z has the remaining columns, and {$$}\gamma{/$$} the remaining parameters.
A>
A>    So all we need to do is redefine {$$}Y{/$$} as {$$}Y^* = Y - X_{male} \beta_{male} - X_{D} \beta_D{/$$} and fit a linear model. The following line of code creates this  variable {$$}Y^*{/$$}, after fixing {$$}\beta_{male}{/$$} to a value `a`, and beta_D to a value, `b`:
A>
A>    
    ```r
    makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b
    ```
A>
A>    Now we'll construct a function which, for a given value a and b, gives us back the the sum of squared residuals after fitting the other terms.
A>
A>    
    ```r
    fitTheRest <- function(a,b) {
      Ystar <- makeYstar(a,b)
       Xrest <- X[,-c(2,5)]
      betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
      residuals <- Ystar - Xrest %*% betarest
      sum(residuals^2)
    }
    ```
A>
A>    What is the sum of squared residuals when the male coefficient is 1 and the D coefficient is 2, and the other coefficients are fit using the linear model solution?
A>
A>
A>
A>
A>3. We can apply our function `fitTheRest` to a grid of values for {$$}\beta_{male}{/$$} and {$$}\beta_D{/$$}, using the `outer` function in R. `outer` takes three arguments, a grid of values for the first argument, a grid of values for the second argument, and finally a function which takes two arguments.
A>
A>    Try it out: 
A>
A>    
    ```r
    outer(1:3,1:3,`*`)
    ```
A>    
A>    We can run fitTheRest on a grid of values, using the following code (the Vectorize() is necessary as outer() requires only vectorized functions):
A>
A>    
    ```r
    outer(-2:8,-2:8,Vectorize(fitTheRest))
    ```
A>    
A>    In the grid of values, what is the smallest sum of squared residuals?
A>    
A>
A>
A>
