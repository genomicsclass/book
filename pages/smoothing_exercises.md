---
Title: Smoothing exercises
---


## Exercises


1. Generate the following data:
  
    
    ```r
    n = 10000
    set.seed(1)
    men = rnorm(n,176,7) #height in centimeters
    women = rnorm(n,162,7) #height in centimeters
    y = c(rep(0,n),rep(1,n))
    x = round(c(men,women))
    ##mix it up
    ind = sample(seq(along=y))
    y = y[ind]
    x = x[ind]
    ```

    Set the seed at 5, `set.seed(5)` and take a random sample of 250 from:

    
    ```r
    set.seed(5)
    N = 250
    ind = sample(length(y),N)
    Y = y[ind]
    X = x[ind]
    ```

    Use loess to estimate $$f(x)=E(Y|X=x)$$ using the default parameters. What is the predicted $$f(168)$?



2. The loess estimate above is a random variable. We can compute standard errors for it. Here we use Monte Carlo to demonstrate that it is a random variable. Use Monte Carlo simulation to estimate the standard error of your estimate of $$f(168)$$. 

    Set the seed to 5, `set.seed(5)` and perform 10000 simulations and report the SE of the loess based estimate.




