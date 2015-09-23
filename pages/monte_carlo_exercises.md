## Exercises

We have used Monte Carlo simulation throughout this chapter to demonstrate statistical concepts; namely, sampling from the population. We mostly applied this to demonstrate the statistical properties related to inference on differences in averages. Here, we will consider examples of how Monte Carlo simulations are used in practice. 

1. Imagine you are [William_Sealy_Gosset](https://en.wikipedia.org/wiki/William_Sealy_Gosset) and have just mathematically derived the distribution of the t-statistic when the sample comes from a normal distribution. Unlike Gosset you have access to computers and can use them to check the results. 

    Let's start by creating an outcome.
Set the seed at 1, use `rnorm` to generate a random sample of size 5, $$X_1, \dots, X_5$$ from a standard normal distribution, then compute the t-statistic $$t = \sqrt{5} \, \bar{X} / s$$ with $$s$$ the sample standard deviation. What value do you observe?



2. You have just performed a Monte Carlo simulation using `rnorm` , a random number generator for normally distributed data. Gosset's mathematical calculation tells us that this random variable follows a t-distribution with $$N-1$$ degrees of freedom. Monte Carlo simulations can be used to check the theory: we generate many outcomes and compare them to the theoretical result. Set the seed to 1, generate $$B=1000$$ t-statistics as done in exercise 1. What percent are larger than 2?


3. The answer to exercise 2 is very similar to the theoretical prediction: `1-pt(2,df=4)`. We can check several such quantiles using the `qqplot` function. 

    To obtain quantiles for the t-distribution we can generate percentiles from just above 0 to just below 1: `B=100; ps = seq(1/(B+1), 1-1/(B+1),len=B)` and compute the quantiles with `qt(ps,df=4)`. Now we can use `qqplot` to compare these theoretical quantiles to those obtained in the Monte Carlo simulation. Use Monte Carlo simulation developed for exercise 2 to corroborate that the t-statistic $$t = \sqrt{N} \, \bar{X} / s$$ follows a t-distribution for several values of $$N$$. 

    For which sample sizes does the approximation best work?
    
    - A) Larger sample sizes.
    - B) Smaller sample sizes.
    - C) The approximations are spot on for all sample sizes.
    - D) None. We should use CLT instead.




4. Use Monte Carlo simulation to corroborate that the t-statistic comparing two means and obtained with normally distributed (mean 0 and sd) data follows a t-distribution. In this case we will use the `t.test` function with `var.equal=TRUE`. With this argument the degrees of freedom will be `df=2*N-2` with `N` the sample size.  For which sample sizes does the approximation best work?
    - A) Larger sample sizes.
    - B) Smaller sample sizes.
    - C) The approximations are spot on for all sample sizes.
    - D) None. We should use CLT instead.



5. Is the following statement true or false? If instead of generating the sample with `X=rnorm(15)`, we generate it with binary data `X=rbinom(n=15,size=1,prob=0.5)` then the t-statistic

    
    ```r
    tstat <- sqrt(15)*mean(X) / sd(X)
    ```

    is approximated by a t-distribution with 14 degrees of freedom. 


6. Is the following statement true or false? If instead of generating the sample with `X=rnorm(N)` with `N=500`, we generate the data with binary data `X=rbinom(n=500,size=1,prob=0.5)`, then the t-statistic `sqrt(N)*mean(X)/sd(X)` is approximated by a t-distribution with 499 degrees of freedom. 



7. We can derive approximation of the distribution of the sample average or the t-statistic theoretically. However, suppose we are interested in the distribution of a statistic for which a theoretical approximation is not immediately obvious. 

    Consider the sample median as an example. Use a Monte Carlo to determine which of the following best approximates the median of a sample taken from normally distributed population with mean 0 and standard deviation 1.
    
    - A) Just like for the average, the sample median is approximately normal with mean 0 and SD $$1/\sqrt{N}$$.
    - B) The sample median is not approximately normal.
    - C) The sample median is t-distributed for small samples and normally distributed for large ones.
    - D) The sample median is approximately normal with mean 0 and SD larger than $$1 / \sqrt{N}$$.


