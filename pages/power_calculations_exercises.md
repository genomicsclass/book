
A>## Exercises
A>
A>For these exercises we will load the babies dataset from `babies.txt`. We will use these data to review the concepts behind the p-values and then test confidence interval concepts.
A>
A>
```r
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
```
A>
A>This is a large dataset (1,236 cases), and we will pretend that it contains the entire population that we are interested in. We will study the differences in birthweight between babies born to smoking and non-smoking mothers.
A>
A>First, let's split this into two birthweight datasets, one of birthweights to non-smoking mothers, and the other of birthweights to smoking mothers.
A>
A>
```r
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
```
A>
A>Now, we can look for the true population difference in means between smoking and non-smoking birthweights.
A>
A>
```r
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)
```
A>
A>The population difference of mean birthweights is about 8.9 ounces. The standard deviations of the nonsmoking and smoking groups are about 17.4 and 18.1 ounces, respectively.
A>
A>As we did with the mouse weight data, this assessment interactively reviews inference concepts using simulations in R. We will treat the babies dataset as the full population, and draw samples from it to simulate individual experiments. We will then ask whether somebody who only received the random samples would be able to draw correct conclusions about the population. 
A>
A>We are interested in testing whether the birth weights of babies born to non-smoking mothers are significantly different from the birth weights of babies born to smoking mothers.
A>
A>1. Set the seed at 1 and obtain two samples, each of size {$$}N=25{/$$}, from non-smoking mothers (`dat.ns`) and smoking mothers (`dat.s`). Following lecture, compute the t-statistic (call it `tval`).
A>
A>
A>
A>2. Recall that we summarize our data using a t-statistics because we know that in situations where the null hypothesis is true (what we mean when we say "under the null") and the sample size is relatively large, this t-value will have an approximate standard normal distribution. Because we know the distribution of the t-value under the null, we can quantitatively determine how unusual the observed t-value would be if the null hypothesis were true. 
A>
A>    The standard procedure is to examine the probability a t-statistic that actually does follow the null hypothesis would have larger absolute value than the absolute value of thet-value we just observed -- this is called a two-sided test.
A>
A>    We have computed these by computing the probability is to take one minus the area under the standard normal curve between `-abs(tval)` and `abs(tval)`. In R, we can do this using the `pnorm` function, which computes the area under a normal curve from negative infinity up to the value given as its first argument:
A>
A>
A>
A>3. Because of the symmetry of the standard normal distribution, there is a simpler way to calculate the probability that a t-value under the null could have a larger absolute value than `tval`. Choose the simplified calculation from the following:
A>    - A) `1-2*pnorm(abs(tval))`
A>    - B) `1-2*pnorm(-abs(tval))`
A>    - C) `1-pnorm(-abs(tval))`
A>    - D) `2*pnorm(-abs(tval))`
A>
A>
A>4. By reporting only p-values, many scientific publications provide an incomplete story of their findings. As we have mentioned, with vey large sample sizes scientifically insignificant differences between two groups can lead to small p-values. Confidence intervals are more informative as they include the estimate itself. Our estimate of the difference between babises of smoker and non-smokers: `mean(dat.s) - mean( dat.ns)`. If we use the CLT, what quantity would we add and substract to this estimate to obtain a 99% confidence interval
A>
A>
A>5. If instead of CLT, we use the t-distribution approximation, what do we add and subtract (use `2*N-2` degrees of freedom)
A>
A>
A>
A>6. Why are the values from 4 and 5 so similar?
A>    - A) Coincidence
A>    - B) They are both related to 99% confidence intervals
A>    - C) `N` and, thus the degrees of freedom, is large enough to make the normal and t-distributions very similar
A>    - D) They are actually quite different, differing by more than 1 ounce.
A>  
A>7. No matter which way you compute it, the p-value `pval` is the probability that the null hypothesis could have generated a t-statistic more extreme than than what we observed: `tval`. If the p-value is very small, this means that observing a value more extreme than `tval` would be very rare if the null hypothesis were true, and would give strong evidence that we should **reject** the null hypothesis. We determine how small the p-value needs to be to reject the null by deciding how often we would be willing to mistakenly reject the null hypothesis.
A>
A>    The standard decision rule is the following: Choose some small value {$$}\alpha{/$$} (in most disciplines the conventional choice is  {$$}\alpha = 0.05{/$$}), and reject the null hypothesis if the p-value is less than {$$}\alpha{/$$}. We call {$$}\alpha{/$$} the _significance level_ of the test.
A>
A>    It turns out if we follow this decision rule, the probability that we will reject the null hypothesis by mistake is equal to {$$}\alpha{/$$}. (This fact is not immediately obvious and requires some probability theory to show.) We call the _event_ of rejecting the null hypothesis when it is in fact true a _Type I error_, we call the _probability_ of making a Type I error the _Type I error rate_, and we say that rejecting the null hypothesis when the p-value is less than {$$}\alpha{/$$} _controls_ the Type I error rate so that it is equal to {$$}\alpha{/$$}. (Over the course of this class, we will see a number of decision rules that we use in order to control the proabilities of other types of errors. Often we will guarantee that the probability of an error is less than some level, but in this case, we can guarnatee that the probability of a Type I error is _exactly equal_ to {$$}\alpha{/$$}.)
A>
A>    Which of the following sentences about a Type I error is **not** true?
A>    
A>    - A) The following is another way to describe Type I error: You decided to reject the null hypothesis on the basis of data that was actually generated by the null hypothesis.
A>    - B) The following is the another way to describe Type I error: Due to random fluctuatons, even though the data you observed were actually generated by the null hypothesis, the p-value calculated from the observed data was small, so you rejected it.
A>    - C) From the original data alone, you can tell whether you have made a Type I error.
A>    - D) In scientific practice, a Type I error constitutes reporting a "significant" result when there is actually no result.
A>    
A>
A>
A>8. In the simulation we have set up here, we know the null hypothesis is false -- the true value of difference in means is actually around  {$$}8.9{/$$}. Thus, we are concerned with how often the decision rule outlined in the last section allows us to conclude the that the null hypothesis is actually false. Put another way, we would like to quantify the _Type II error rate_ of the test, or the probability that we fail to reject the null hypothesis when the alternative hypothesis is true.
A>
A>    Unlike the Type I error rate, which we can characterize by assuming that the null hypothesis of "no difference" is true, the Type II error rate cannot be computed by assuming the alternative hypothesis alone because the alternative hypothesis alone does not specify a particular value for the difference, and thus does not nail down a specific distrbution for the t-value under the alternative.
A>
A>    For this reason, when we study the Type II error rate of a hypothesis testing procedure, we need to assume a particular _effect size_, or hypothetical size of the difference between population means, that we wish to target. We ask questions like "What is the smallest difference I could reliably distinguish from 0 given my sample size {$$}N{/$$}?", or more commonly, "How big does {$$}N{/$$} have to be in order to detect that the absolute value of the difference is greater than zero?" Type II error control plays a major role in designing data collection procedures **before** you actually see the data so that you know the test you will run has enough sensitivity or _power_. Power is one minus the Type II error rate, or the probability that you will reject the null hypothesis when the alternative hypothesis is true.
A>
A>    There are several aspects of a hypothesis test that affect its power for a particular effect size. Intuitively, setting a lower {$$}\alpha{/$$} decreases the power of the test for a given effect size because the null hypothesis will be more difficult to reject. This means that for an experiment with fixed parameters (i.e., with a predetermined sample size, recording mechanism, etc), the power o the hypothesis test trades off with its Type I error rate, no matter what effect size you target.
A>
A>    We can explore the tradeoff of power and Type I error concretely using the babies data. Because we have the full population, we know what the true effect size is (about 8.93), and we can compute the power of the test for true difference between populations.
A>
A>    Set the seed at 1 and take a random sample of {$$}N=5{/$$} measurements from each of the smoking and nonsmoking datasets. What is the p-value (use the `t-test` function)?
A>
A>
A>9. Note that the p-value is larger than 0.05 so using the typical cut-off, we would no reject. This is a type II error. Which of the following is *not* a way to decrease  this type of error?
A>    - A) Increase our chance of a type I error 
A>    - B) Take a larger sample size 
A>    - C) Find a population for which the null is not true
A>    - D) Use a higher {$$}\alpha{/$$} level
A>
A>
A>
A>10. Set the seed at 1, then use the `replicate` function to repeat the code us used in exercise 9 10,000 times. What proportion of the time do we reject at the 0.05 level?
A>
A>
A>
A>11. Note that, not surprisingly, the power is lower than 10%. Repeat the exercise above for samples sizes of 30, 60, 90 and 120. Which of those four gives you power of about 80%?
A>
A>
A>
A>11. Repear problem 11, but now require an {$$}\alpha{/$$} level of 0.01.  Which of those four gives you power of about 80%?
A>
A>
