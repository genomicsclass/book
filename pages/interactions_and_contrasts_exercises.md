---
Title: Interactions and Contrasts
---

## Exercises

1. Suppose we have an experiment with two species A and B, and two conditions: control and treated.


```r
species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
```

And we will use a formula of '~ species + condition'.

The model matrix is then:


```r
model.matrix(~ species + condition)
```

1. Suppose we want to build a contrast of coefficients for the above experimental design.

    You can either figure this question out through logic, by looking at the design matrix, or using the `contrast` function from the contrast library. The contrast vector is returned as `contrast(...)$X`.

    What should the contrast vector be, for the contrast of (species=B and condition=control) vs (species=A and condition=treatment)? Assume that the beta vector from the model fit by R is: Intercept, speciesB, conditiontreated.
    - A) 0 0 1  
    - B) 0 -1 0 
    - C) 0 1 1  
    - D) 0 1 -1
    - E. 0 -1 1
    - F. 1 0 1
    

   
   
2. Use the Rmd script of the spider dataset. Suppose we build a model using two `variables: ~ type + leg`.

    What is the t-statistic for the contrast of leg pair L4 vs leg pair L2?



3. The t-statistic for the contrast of leg pair L4 vs leg pair L2 is constructed by taking the difference of the coefficients legL4 and legL2, and then dividing by the standard error of the difference. In the last question we will explore how the standard error of the difference is calculated here.

    For a contrast vector $$\mathbf{C}$$, the standard error of the contrast $$\mathbf{C}\hat{\boldsymbol{\beta}}$$ is:
 
    $$
    \sqrt{ \mathbf{C} \boldsymbol{\Sigma} \mathbf{C}^\top }
    $$
    
    $$\Sigma$$ is the covariance matrix of the coeffcient estimates $$\hat{\boldsymbol{\beta}}$$. The covariance matrix contains elements which give the variance or covariance of elements in beta-hat. The elements on the diagonal of the $$\boldsymbol{\Sigma}$$ matrix give the variance of each element in beta-hat. The square root of these is the standard error of the elements in $$\hat{\boldsymbol{\beta}}$$. The off-diagonal elements of Sigma give the covariance of two different elements of the $$\hat{\boldsymbol{\beta}}$$ matrix. So $$\boldsymbol{\Sigma}[1,2]$$ gives the covariance of the first and second element of $$\hat{\boldsymbol{\beta}}$$. The $$\hat{\boldsymbol{\beta}}$$ matrix is symmetric, which means $$\boldsymbol{\Sigma}[i,j]=\boldsymbol{\Sigma}[j,i]$$.
    

    $$
    \mbox{var}(\hat{\beta }_{L4} - \hat{\beta }_{L2}) = \mbox{var}(\hat{\beta }_{L4}) + \mbox{var}(\hat{\beta }_{L2}) - 2 \mbox{cov}(\hat{\beta }_{L4}, \hat{\beta }_{L2})
    $$
    
    In the book page, we computed Sigma using:

    
    ```r
    X <- model.matrix(~ type + leg, data=spider)
    Sigma <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X)
    ```

    Our contrast matrix is:
    
    ```r
    C <- matrix(c(0,0,-1,0,1),1,5)
    ```
    
    Using $$\boldsymbol{\Sigma}$$, what is $$\mbox{cov}(\hat{\beta}_{L4}, \hat{\beta}_{L2})$$ ?



    To see a further application of the matrix algebra, we can confirm that 

    $$
\mbox{var}(\hat{\beta}_{L4} - \hat{\beta}_{L2}) = \mbox{var}(\hat{\beta}_{L4}) + \mbox{var}(\hat{\beta}_{L2} - 2 \mbox{cov}(\hat{\beta}_{L4}, \hat{\beta}_{L2})
    $$

    is equal to
    
    ```r
    C %*% Sigma %*% t(C))
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
    
    What is the t-statistic for the interaction of type push and leg L4? If this t-statistic is sufficiently large, we would reject the null hypothesis that the push vs pull effect on `log2(friction)` is the same in L4 as in L1.



5. What is the F-value for all of the type:leg interaction terms, in an analysis of variance? If this value is sufficiently large, we would reject the null hypothesis that the push vs pull effect on `log2(friction)` is the same for all leg pairs.




6. What is the L2 vs L1 estimate in log2friction for the pull samples?



7. What is the L2 vs L1 estimate in log2friction for the push samples? Remember, because of the interaction terms, this is not the same as the L2 vs L1 difference for the pull samples. If you're not sure use the `contrast` function. Another hint: consider the arrows plot for the model with interactions.


   
    Note that taking the log2 of a Y value and then performing a linear model has a meaningful effect on the coefficients. If we have, $$\log_2(X) = \beta_0$$ and $$log2(Y) = \beta_0 + \beta_1$$, then $$Y/X = 2^(\beta_0 + \beta_1) / 2^(\beta_0)= 2^\beta_1$$, so $$\beta_1$$ represents a log2 fold change of $$Y$$ over $$X$$. If $$\beta_1 = 1$$, then $$Y$$ is 2 times $$X$$. If $$\beta_1 = -1$$, then $$Y$$ is half of $$X$$, etc. 

8. Snalysis of variance (ANOVA), performed in R using the `anova`, allows us to test whether a number of coefficients are equal to zero, by comparing a linear model including these terms to a linear model where these terms are set to 0. In this last question, we will use Monte Carlo techniques to observe the distribution of the ANOVA's _F-value_ under the null hypothesis, that there are no differences between groups.

    Suppose we have 4 groups, and 10 samples per group, so 40 samples overall:

    
    ```r
    N <- 40
    p <- 4
    group <- factor(rep(1:p,each=N/p))
    X <- model.matrix(~ group)
    ```

    We will show here how to calculate the "F-value", and then we will use random number to observe the distribution of the F-value under the null hypothesis.

    The F-value is the mean sum of squares explained by the terms of interest (in our case, the 'group' terms) divided by the mean sum of squares of the residuals of a model including the terms of interest. So it is the explanatory power of the terms divided by the leftover variance.

    Intuitively, if this number is large, it means that the group variable explains a lot of the variance in the data, compared to the amount of variance left in the data after using group information. We will calculate these values exactly here:

    First generate some random, null data, where the mean is the same for all groups:

    
    ```r
    Y <- rnorm(N,mean=42,7)
    ```
    
    The base model we wil compare against is simply $$\hat{Y}$$ = `mean(Y)`, which we will call $$\mu_0$$, and the initial sum of squares is the Y values minus $$\mu_0$$:

    
    ```r
    mu0 <- mean(Y)
    initial.ss <- sum((Y - mu0)^2)
    ```
    
    We then need to calculate the fitted values for each group, which is simply the mean of each group, and the residuals from this model, which we will call "after.group.ss" for the sum of squares after using the group information:

    
    ```r
    s <- split(Y, group)
    after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
    ```
    
    Then the explanatory power of the group variable is the initial sum of squares minus the residual sum of squares:

    
    ```r
    (group.ss <- initial.ss - after.group.ss)
    ```

    We calculate the mean of these values, but we divide by terms which remove the number of fitted parameters. For the group sum of squares, this is number of parameters used to fit the groups (3, because the intercept is in the initial model). For the after group sum of squares, this is the number of samples minus the number of parameters total (So $$N - 4$$, including the intercept).
    
    ```r
    group.ms <- group.ss / (p - 1)
    after.group.ms <- after.group.ss / (N - p)
    ```

    The F-value is simply the ratio of these mean sum of squares.

    
    ```r
    f.value <- group.ms / after.group.ms
    ```
    What Is the point of all these calculations? The point is that, after following these steps, the exact distribution of the F-value has a nice mathematical formula under the null hypothesis. 
    
    Set the seed at 1 and calculate the F-value for 1000 random versions of Y. What is the mean of these F-values?



    If you save the values from the simulation into `Fs` we can plot the distribution of the 1000 F-values and overlay the theoretical F-distribution, with parameters `df1=p - 1`, `df2=N - p`. Note the similarity:

    
    ```r
    hist(Fs, col="grey", border="white", breaks=50, freq=FALSE)
    xs <- seq(from=0,to=6,length=100)
    lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")
    ```

    This is the distribution which is used to calculate the p-values for the ANOVA table produced by `anova`. 
