---
Title: Conditional Expectations Exercises
---

## Exercises

Throughout this exercises it will be useful to remember that when our data are 0s and 1s, probabilities and expectations are the same thing. We can do the math, but here is some R code:


```r
n = 1000
y = rbinom(n,1,0.25)
##proportion of ones Pr(Y)
sum(y==1)/length(y)
##expectaion of Y
mean(y)
```


1. Generate some random data to imitate heights for men (0) and women (1):

    
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

1. Using the data generated above, what is the $$E(Y|X=176)$?


2. Now make a plot of $$E(Y|X=x)$$ for `x=seq(160,178)` using the data generated in exercise 1.

    If you want your chance of predicting female or male based on height, and you want your probability of success to be larger than 0.5, at what is the largest height where you predict female ?


