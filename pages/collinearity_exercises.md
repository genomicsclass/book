---
layout: page
Title: Collinearity  exercises
---

{pagebreak} 

## Exercises


Consider these design matrices:

$$
A=\begin{pmatrix}
1&0& 0& 0\\
1& 0& 0& 0\\
1& 1& 1& 0\\
1& 1& 0& 1
\end{pmatrix}\,\,
B=\begin{pmatrix}
1& 0& 0& 1\\
1& 0& 1& 1\\
1& 1& 0& 0\\
1& 1& 1& 0
\end{pmatrix}\,\,
C=\begin{pmatrix}
1& 0& 0\\
1& 1& 2\\
1& 2& 4\\
1& 3& 6
\end{pmatrix} 
$$

$$
D=\begin{pmatrix}
1& 0& 0& 0& 0\\
1& 0& 0& 0& 1\\
1& 1& 0& 1& 0\\
1& 1& 0& 1& 1\\
1& 0& 1& 1& 0\\
1& 0& 1& 1& 1
\end{pmatrix}\,\,
E=\begin{pmatrix}
1& 0& 0& 0\\
1& 0& 1& 0\\
1& 1& 0& 0\\
1& 1& 1& 1
\end{pmatrix} \,\,
F=\begin{pmatrix}
1& 0& 0& 1\\
1& 0& 0& 1\\
1& 0& 1& 1\\
1& 1& 0& 0\\
1& 1& 0& 0\\
1& 1& 1& 0
\end{pmatrix}
$$


1. Which of the above design matrices does NOT have the problem of collinearity?


2. The following exercises are advanced. Let's use the example from the lecture to visualize how there is not a single best $$\hat{\beta}$$, when the design matrix has collinearity of columns. An example can be made with:

    
    ```r
    sex <- factor(rep(c("female","male"),each=4))
    trt <- factor(c("A","A","B","B","C","C","D","D"))
    ```

    The model matrix can then be formed with:

    
    ```r
    X <- model.matrix( ~ sex + trt)
    ```

    And we can see that the number of independent columns is less than the number of columns of X:

    
    ```r
    qr(X)$rank
    ```
    Suppose we observe some outcome Y. For simplicity, we will use synthetic data:

    
    ```r
    Y <- 1:8
    ```

    Now we will fix the value for two coefficients and optimize the remaining ones. We will fix $$\beta_{male}$$ and $$\beta_D$$. Then we will find the optimal value for the remaining betas, in terms of minimizing the residual sum of squares. We find the value that minimize:

    $$
\sum_{i=1}^8  \{ (Y_i - X_{i,male} \beta_{male} - X_{i,D} \beta_{i,D}) - \mathbf{Z}_i \boldsymbol{\gamma} )^2 \}
    $$

    where $$X_{male}$$ is the male column of the design matrix, $$X_D$$ is the D column, $$\mathbf{Z}_i$$ is a 1 by 3 matrix with the remaining column entries for unit $$i$$, and $$\boldsymbol{\gamma}$$ is a 3 x 1 matrix with the remaining parameters.

    So all we need to do is redefine $$Y$$ as $$Y^* = Y - X_{male} \beta_{male} - X_{D} \beta_D$$ and fit a linear model. The following line of code creates this  variable $$Y^*$$, after fixing $$\beta_{male}$$ to a value `a`, and $$\beta_D$$ to a value, `b`:

    
    ```r
    makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b
    ```

    Now we'll construct a function which, for a given value a and b, gives us back the sum of squared residuals after fitting the other terms.

    
    ```r
    fitTheRest <- function(a,b) {
      Ystar <- makeYstar(a,b)
      Xrest <- X[,-c(2,5)]
      betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
      residuals <- Ystar - Xrest %*% betarest
      sum(residuals^2)
    }
    ```

    What is the sum of squared residuals when the male coefficient is 1 and the D coefficient is 2, and the other coefficients are fit using the linear model solution?




3. We can apply our function `fitTheRest` to a grid of values for $$\beta_{male}$$ and $$\beta_D$$, using the `outer` function in R. `outer` takes three arguments: a grid of values for the first argument, a grid of values for the second argument, and finally a function which takes two arguments.

    Try it out: 

    
    ```r
    outer(1:3,1:3,`*`)
    ```
    
    We can run `fitTheRest` on a grid of values, using the following code (the `Vectorize` is necessary as `outer` requires only vectorized functions):

    
    ```r
    outer(-2:8,-2:8,Vectorize(fitTheRest))
    ```
    
    In the grid of values, what is the smallest sum of squared residuals?
    



