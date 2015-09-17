---
Title: Standard Errors Exercises
---



A>## Exercises
A>
A>In the previous assessment, we used a Monte Carlo technique to see that the linear model coefficients are random variables when the data is a random sample. Now we will use the matrix algebra from the previous video to try to estimate the standard error of the linear model coefficients. Again, take a random sample of the father.son heights data:
A>
A>
```r
library(UsingR)
N <- 50
set.seed(1)
index <- sample(n,N)
sampledat <- father.son[index,]
x <- sampledat$fheight
y <- sampledat$sheight
betahat <- lm(y~x)$coef
```
A>
A>The formula for the standard error is:
A>
A>{$$}
A>\mbox{SE}(\hat{\beta}) = \sqrt{\mbox{var}(\hat{\beta})}
A>{/$$}
A>
A>with
A>
A>{$$}
A>\mbox{var}(\hat{\beta}) = \sigma^2 (X^\top X)^{-1}
A>{/$$}
A>
A>We will estimate or calculate each part of this equation and then combine them.
A>
A>First, we want to estimate {$$}\sigma^2{/$$}, the variance of {$$}Y{/$$}. As we have seen in the previous unit, the random part of {$$}Y{/$$} is only coming from epsilon, because we assume {$$}X\beta{/$$} is fixed. So we can try to estimate the variance of the epsilons from the residuals, the {$$}Y_i{/$$} minus the fitted values from the linear model.
A>
A>
A>1. Note that the fitted values (Y-hat) from a linear model can be obtained with:
A>
A>    
    ```r
    fit <- lm(y ~ x)
    fit$fitted.values
    ```
A>    
A>    What is the sum of the squared residuals (where residuals are given by {$$}r_i = Y_i - \hat{Y}_i{/$$} ?
A>
A>
A>
A>2. Our estimate of {$$}\sigma^2{/$$} will be the sum of squared residuals divided by {$$}N - p{/$$}, the sample size minus the number of terms in the model. Since we have a sample of 50 and 2 terms in the model (an intercept and a slope), our estimate of {$$}\sigma^2{/$$} will be the sum of squared residuals divided by 48. Use the answer to exercise 1 to provide an estimate of {$$}\sigma^2{/$$}:
A>
A>
A>  
A>3. Form the design matrix X (note: use a capital X!). This can be done by combining a column of 1's with a 
A>column of 'x' the father's heights.
A>
A>    
    ```r
    N <- 50
    X <- cbind(rep(1,N), x)
    ```
A>
A>    Now calculate {$$}(X^\top X)^{-1}{/$$}. Use the `solve` function for the inverse and `t` for the transpose. What is the element in the first row, first column?
A>
A>
A>4. Now we are one step away from the standard error of beta-hat. Take the diagonals from the {$$}(X^\top X)^{-1}{/$$} matrix above, using the `diag` function. Now multiply our estimate of {$$}\sigma^2{/$$} and the diagonals of this matrix. This is the estimated variance of {$$}\hat{beta}{/$$}, so take the square root of this. You should end up with two numbers, the standard error for the intercept and the standard error for the slope.
A>
A>    What is the standard error for the slope?
A>
A>
A>
A>Compare your answer to the last question to the value you estimated using Monte Carlo in the previous set of exercises. It will not be the same, because we are only estimating the standard error given a particular sample of 50 (which we obtained with set.seed(1)).
A>
A>Note that the standard error estimate is also printed in the second column of:
A>
A>
```r
summary(fit)
```
