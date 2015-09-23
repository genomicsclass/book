---
layout: page
title:  Power Calculations Exercises
---

## Exercises

For these exercises we will load the babies dataset from `babies.txt`. We will use this data to review the concepts behind the p-values and then test confidence interval concepts.


```r
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
```

This is a large dataset (1,236 cases), and we will pretend that it contains the entire population in which we are interested. We will study the differences in birth weight between babies born to smoking and non-smoking mothers.

First, let's split this into two birth weight datasets: one of birth weights to non-smoking mothers and the other of birth weights to smoking mothers.


```r
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
```

Now, we can look for the true population difference in means between smoking and non-smoking birth weights.


```r
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)
```

The population difference of mean birth weights is about 8.9 ounces. The standard deviations of the nonsmoking and smoking groups are about 17.4 and 18.1 ounces, respectively.

As we did with the mouse weight data, this assessment interactively reviews inference concepts using simulations in R. We will treat the babies dataset as the full population and draw samples from it to simulate individual experiments. We will then ask whether somebody who only received the random samples would be able to draw correct conclusions about the population. 

We are interested in testing whether the birth weights of babies born to non-smoking mothers are significantly different from the birth weights of babies born to smoking mothers.

1. Set the seed at 1 and obtain two samples, each of size $$N=25$$, from non-smoking mothers (`dat.ns`) and smoking mothers (`dat.s`). Compute the t-statistic (call it `tval`).



2. Recall that we summarize our data using a t-statistics because we know that in situations where the null hypothesis is true (what we mean when we say "under the null") and the sample size is relatively large, this t-value will have an approximate standard normal distribution. Because we know the distribution of the t-value under the null, we can quantitatively determine how unusual the observed t-value would be if the null hypothesis were true. 

    The standard procedure is to examine the probability a t-statistic that actually does follow the null hypothesis would have larger absolute value than the absolute value of the t-value we just observed -- this is called a two-sided test.

    We have computed these by taking one minus the area under the standard normal curve between `-abs(tval)` and `abs(tval)`. In R, we can do this by using the `pnorm` function, which computes the area under a normal curve from negative infinity up to the value given as its first argument:



3. Because of the symmetry of the standard normal distribution, there is a simpler way to calculate the probability that a t-value under the null could have a larger absolute value than `tval`. Choose the simplified calculation from the following:
    - A) `1-2*pnorm(abs(tval))`
    - B) `1-2*pnorm(-abs(tval))`
    - C) `1-pnorm(-abs(tval))`
    - D) `2*pnorm(-abs(tval))`


4. By reporting only p-values, many scientific publications provide an incomplete story of their findings. As we have mentioned, with very large sample sizes, scientifically insignificant differences between two groups can lead to small p-values. Confidence intervals are more informative as they include the estimate itself. Our estimate of the difference between babies of smoker and non-smokers: `mean(dat.s) - mean( dat.ns)`. If we use the CLT, what quantity would we add and subtract to this estimate to obtain a 99% confidence interval?


5. If instead of CLT, we use the t-distribution approximation, what do we add and subtract (use `2*N-2` degrees of freedom)?



6. Why are the values from 4 and 5 so similar?
    - A) Coincidence.
    - B) They are both related to 99% confidence intervals.
    - C) `N` and thus the degrees of freedom is large enough to make the normal and t-distributions very similar.
    - D) They are actually quite different, differing by more than 1 ounce.
  
7. No matter which way you compute it, the p-value `pval` is the probability that the null hypothesis could have generated a t-statistic more extreme than than what we observed: `tval`. If the p-value is very small, this means that observing a value more extreme than `tval` would be very rare if the null hypothesis were true, and would give strong evidence that we should **reject** the null hypothesis. We determine how small the p-value needs to be to reject the null by deciding how often we would be willing to mistakenly reject the null hypothesis.

    The standard decision rule is the following: choose some small value $$\alpha$$ (in most disciplines the conventional choice is  $$\alpha = 0.05$$) and reject the null hypothesis if the p-value is less than $$\alpha$$. We call $$\alpha$$ the _significance level_ of the test.

    It turns out that if we follow this decision rule, the probability that we will reject the null hypothesis by mistake is equal to $$\alpha$$. (This fact is not immediately obvious and requires some probability theory to show.) We call the _event_ of rejecting the null hypothesis, when it is in fact true, a _Type I error_, we call the _probability_ of making a Type I error, the _Type I error rate_, and we say that rejecting the null hypothesis when the p-value is less than $$\alpha$$, _controls_ the Type I error rate so that it is equal to $$\alpha$$. We will see a number of decision rules that we use in order to control the probabilities of other types of errors. Often, we will guarantee that the probability of an error is less than some level, but, in this case, we can guarantee that the probability of a Type I error is _exactly equal_ to $$\alpha$$.

    Which of the following sentences about a Type I error is **not** true?
    
    - A) The following is another way to describe a Type I error: you decided to reject the null hypothesis on the basis of data that was actually generated by the null hypothesis.
    - B) The following is the another way to describe a Type I error: due to random fluctuations, even though the data you observed were actually generated by the null hypothesis, the p-value calculated from the observed data was small, so you rejected it.
    - C) From the original data alone, you can tell whether you have made a Type I error.
    - D) In scientific practice, a Type I error constitutes reporting a "significant" result when there is actually no result.
    


8. In the simulation we have set up here, we know the null hypothesis is false -- the true value of difference in means is actually around  $$8.9$$. Thus, we are concerned with how often the decision rule outlined in the last section allows us to conclude that the null hypothesis is actually false. In other words, we would like to quantify the _Type II error rate_ of the test, or the probability that we fail to reject the null hypothesis when the alternative hypothesis is true.

    Unlike the Type I error rate, which we can characterize by assuming that the null hypothesis of "no difference" is true, the Type II error rate cannot be computed by assuming the alternative hypothesis alone because the alternative hypothesis alone does not specify a particular value for the difference. It thus does not nail down a specific distribution for the t-value under the alternative.

    For this reason, when we study the Type II error rate of a hypothesis testing procedure, we need to assume a particular _effect size_, or hypothetical size of the difference between population means, that we wish to target. We ask questions such as "what is the smallest difference I could reliably distinguish from 0 given my sample size $$N$?" or, more commonly, "How big does $$N$$ have to be in order to detect that the absolute value of the difference is greater than zero?" Type II error control plays a major role in designing data collection procedures **before** you actually see the data, so that you know the test you will run has enough sensitivity or _power_. Power is one minus the Type II error rate, or the probability that you will reject the null hypothesis when the alternative hypothesis is true.

    There are several aspects of a hypothesis test that affect its power for a particular effect size. Intuitively, setting a lower $$\alpha$$ decreases the power of the test for a given effect size because the null hypothesis will be more difficult to reject. This means that for an experiment with fixed parameters (i.e., with a predetermined sample size, recording mechanism, etc), the power of the hypothesis test trades off with its Type I error rate, no matter what effect size you target.

    We can explore the trade off of power and Type I error concretely using the babies data. Since we have the full population, we know what the true effect size is (about 8.93) and we can compute the power of the test for true difference between populations.

    Set the seed at 1 and take a random sample of $$N=5$$ measurements from each of the smoking and nonsmoking datasets. What is the p-value (use the `t-test` function)?


9. The p-value is larger than 0.05 so using the typical cut-off, we would not reject. This is a type II error. Which of the following is *not* a way to decrease this type of error?
    - A) Increase our chance of a type I error.
    - B) Take a larger sample size.
    - C) Find a population for which the null is not true.
    - D) Use a higher $$\alpha$$ level.



10. Set the seed at 1, then use the `replicate` function to repeat the code used in exercise 9 10,000 times. What proportion of the time do we reject at the 0.05 level?



11. Note that, not surprisingly, the power is lower than 10%. Repeat the exercise above for samples sizes of 30, 60, 90 and 120. Which of those four gives you power of about 80%?



11. Repeat problem 11, but now require an $$\alpha$$ level of 0.01.  Which of those four gives you power of about 80%?


