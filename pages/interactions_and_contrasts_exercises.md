---
layout: page
Title: Interactions and Contrasts
---

## Exercises

Suppose we have an experiment with two species A and B, and two conditions, control and treated.


```r
species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
```

We will use the formula of `~ species + condition` to create the model matrix:


```r
model.matrix(~ species + condition)
```

1. Suppose we want to build a contrast of coefficients for the above experimental design.

    You can either figure this question out by looking at the design matrix, or by using the `contrast` function from the contrast library with random numbers for `y`. The contrast vector will be returned as `contrast(...)$X`.

    What should the contrast vector be, to obtain the difference between the species B control group and the species A treatment group (species B control - species A treatment)? Assume that the coefficients (columns of design matrix) are: Intercept, speciesB, conditiontreated.
    
    - A) 0 0 1  
    - B) 0 -1 0 
    - C) 0 1 1  
    - D) 0 1 -1
    - E) 0 -1 1
    - F) 1 0 1
    

   
   
2. Use the Rmd script to load the spider dataset. Suppose we build a model using two `variables: ~ type + leg`.

    What is the t-statistic for the contrast of leg pair L4 vs. leg pair L2?



3. The t-statistic for the contrast of leg pair L4 vs. leg pair L2 is constructed by taking the difference of the estimated coefficients legL4 and legL2, and then dividing by the standard error of the difference.

    For a contrast vector $$\mathbf{C}$$, the standard error of the contrast $$\mathbf{C}\hat{\boldsymbol{\beta}}$$ is:
 
    $$
    \sqrt{ \mathbf{C} \hat{\boldsymbol{\Sigma}} \mathbf{C}^\top }
    $$
    
    with $$\hat{\boldsymbol{\Sigma}}$$ the estimated covariance matrix of the coefficient estimates $$\hat{\boldsymbol{\beta}}$$. The covariance matrix contains elements which give the variance or covariance of elements in \hat{\beta}. The elements on the diagonal of the $$\boldsymbol{\hat{\Sigma}}$$ matrix give the estimated variance of each element in $$\hat{\boldsymbol{\beta}}$$. The square root of these is the standard error of the elements in $$\hat{\boldsymbol{\beta}}$$. The off-diagonal elements of $$\hat{\Sigma}$$ give the estimated covariance of two different elements of $$\hat{\boldsymbol{\beta}}$$. So $$\hat{\Sigma}_{1,2}$$ gives the covariance of the first and second element of $$\hat{\boldsymbol{\beta}}$$. The $$\hat{\boldsymbol{\beta}}$$ matrix is symmetric, which means $$\hat{\Sigma}_{i,j}=\hat{\Sigma}_{j,i}$$.

    For the difference, we have that:

    $$
    \mbox{var}(\hat{\beta }_{L4} - \hat{\beta }_{L2}) = \mbox{var}(\hat{\beta }_{L4}) + \mbox{var}(\hat{\beta }_{L2}) - 2 \mbox{cov}(\hat{\beta }_{L4}, \hat{\beta }_{L2})
    $$
    
    We showed how to estimate $$\hat{\Sigma}$$ using:

    
    ```r
    X <- model.matrix(~ type + leg, data=spider)
    Sigma.hat <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X)
    ```

    Using the estimate of $$\Sigma$$, what is your estimate of  $$\mbox{cov}(\hat{\beta}_{L4}, \hat{\beta}_{L2})$$ ?



    Our contrast matrix for the desired comparison is:
    
    ```r
    C <- matrix(c(0,0,-1,0,1),1,5)
    ```

4. Suppose that we notice that the within-group variances for the groups with smaller frictional coefficients are generally smaller, and so we try to apply a transformation to the frictional coefficients to make the within-group variances more constant.

    Add a new variable `log2friction` to the spider dataframe:

    
    ```r
    spider$log2friction <- log2(spider$friction)
    ```

    The `Y` values now look like:

    
    ```r
    boxplot(log2friction ~ type*leg, data=spider)
    ```

    Run a linear model of log2friction with type, leg and interactions between type and leg.
    
    What is the t-statistic for the interaction of type push and leg L4? If this t-statistic is sufficiently large, we would reject the null hypothesis that the push vs. pull effect on `log2(friction)` is the same in L4 as in L1.



5. Using the same analysis of log2 transformed data, What is the F-value for all of the `type:leg` interaction terms in an ANOVA? If this value is sufficiently large, we would reject the null hypothesis that the push vs. pull effect on `log2(friction)` is the same for all leg pairs.




6. What is the L2 vs. L1 estimate in `log2(friction)` for the pull samples?



7. What is the L2 vs. L1 estimate in `log2(friction)` for the push samples? Remember, because of the interaction terms, this is not the same as the L2 vs L1 difference for the pull samples. If you're not sure use the `contrast` function. Another hint: consider the arrows plot for the model with interactions.


   
