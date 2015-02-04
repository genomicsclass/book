---
title: "Confidence Intervals"
layout: page
---







# Introduction
We have described how to compute p-values, which are ubiquitous in the life sciences. However, we do not recommend reporting p-values as the only statistical summary of your results. The reason is simple: statistical significance does not guarantee scientific significance. With large enough sample sizes, one might detect a statistically significance difference in weight of, say, 1 microgram. But is this an important finding? Would we say a diet results in higher weight if the increase is less than a fraction of a percent? The problem with reporting only p-values is that you will not provide a very important piece of information: the effect sizes.

An much more attractive alternative is to report confidence intervals. A confidence interval includes relays information about your estimated effect size and the uncertainty associated with this estimate. Here we use the mice data to illustrate the concept behind confidence intervals.

# Confidence interval for population mean

We show how to construct a confidence interval for the population mean of control female mice. 
We start by reading in the data and selecting the appropriate rows:


```r
library(downloader)
```

```
## Error in library(downloader): there is no package called 'downloader'
```

```r
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- tempfile()
download(url,destfile=filename)
```

```
## Error in eval(expr, envir, enclos): could not find function "download"
```

```r
dat <- read.csv(filename)
```

```
## Warning in file(file, "rt"): cannot open file
## '/var/folders/6d/d_8pbllx7318htlp5wv_rm580000gn/T//RtmpbVIFv7/file2ac5c9696c3':
## No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
chowPopulation <- dat[dat$Sex=="F" & dat$Diet=="chow",3]
```

```
## Error in eval(expr, envir, enclos): object 'dat' not found
```

The population average $\mu_X$ is our parameter of interest here:


```r
mu_chow <- mean(chowPopulation)
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
print(mu_chow)
```

```
## Error in print(mu_chow): object 'mu_chow' not found
```

We are interested in estimating this parameter. In practice we do not get to see the entire population so, as we did for p-values, we demonstrate how we can use samples to do this. Let's start with a sample of size 30:


```r
N <- 30
hf <- sample(chowPopulation,N)
```

```
## Error in sample(chowPopulation, N): object 'chowPopulation' not found
```

We know this is a random variable, so the sample average will not be a perfect estimate. In fact, because in this illustrative example and we know the value of the parameter, we can see that they are not exactly the same. A confidence interval is a statistical way of reporting our finding, the sample average, in a way that explicitly summarizes the variability of our random variable.

With a sample size of 30, we will use the CLT. The CLT tells us that $\bar{X}$ or `mean(hf)` follows a normal distribution with mean $\mu_X$ or `mean(chowPopulation)` and standard error approximately  $s_X/\sqrt{N}$ or


```r
se <- sd(hf)/sqrt(N)
```

```
## Error in is.data.frame(x): object 'hf' not found
```

```r
print(se)
```

```
## Error in print(se): object 'se' not found
```

<a name="interval"></a>

## Defining the interval

A 95% confidence interval (we can use percentages other than 95%) is an random intervals with probability 95% of falling on the parameter we are estimating. To construct it we note that CLT tells us that $\sqrt{N} (\bar{X}-\mu_X)/s_X$ follows a normal distribution with mean 0 and SD 1. This implies that the probability that the probability of this event:

$$-2 \leq \sqrt{N} (\bar{X}-\mu_X)/s_X \leq 2$$  


```r
pnorm(2)-pnorm(-2)
```

```
## [1] 0.9544997
```

is about 95% (to get closer use `qnorm(1-0.05/2)` instead of 2). Now do some basic algebra to clear out everything and leave $\mu_X$  alone in the middle and you get that the following event

$$\bar{X}-2 s_X/\sqrt{N} \leq \mu_X \leq \bar{X}+2s_X/\sqrt{N}$$  

has probability 95%. 

Note: it is important to note that it is the edges of the interval $\bar{X} \pm 2 s_X/\sqrt{N}$, not $\mu_X$, that are random. 


What does this mean? We can construct this interval with R relatively easily:

```r
Q <- qnorm(1- 0.05/2)
interval <- c(mean(hf)-Q*se, mean(hf)+Q*se )
```

```
## Error in mean(hf): object 'hf' not found
```

```r
interval
```

```
## Error in eval(expr, envir, enclos): object 'interval' not found
```

Which covers $\mu_X$ or `mean(chowPopulation)`. However, we can take another sample and we might not be as lucky. In fact the theory tells us that we will cover $\mu_X$ 95% of the time. Because we have access to the population data we can confirm this by taking several new samples:


```r
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
B <- 250
mypar2(1,1)
plot(mean(chowPopulation)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
abline(v=mean(chowPopulation))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
for(i in 1:B){
  hf <- sample(chowPopulation,N)
  se=sd(hf)/sqrt(N)
  interval <- c(mean(hf)-Q*se, mean(hf)+Q*se )
  covered<-mean(chowPopulation)<= interval[2] & mean(chowPopulation)>=interval[1]
  color <- ifelse(covered,1,2)
  lines( interval, c(i,i),col=color)
}
```

```
## Error in sample(chowPopulation, N): object 'chowPopulation' not found
```

You can run this over and over again to see what happens. You will see that about 5% we fail to cover $\mu_X$.

<a name="smallsample"></a>

## Small sample size and the CLT

For $N=30$ the CLT works very well. However what if $N=5$, do these confidence interval work as well? We used the CLT to create our intervals, and with $N=5$ it may not be a useful approximation. We can confirm this with a simulation



```r
plot(mean(chowPopulation)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
abline(v=mean(chowPopulation))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
Q <- qnorm(1- 0.05/2)
N<-5
for(i in 1:B){
  hf <- sample(chowPopulation,N)
  se=sd(hf)/sqrt(N)
  interval <- c(mean(hf)-Q*se, mean(hf)+Q*se )
  covered<-mean(chowPopulation)<= interval[2] & mean(chowPopulation)>=interval[1]
  color <- ifelse(covered,1,2)
  lines( interval, c(i,i),col=color)
}
```

```
## Error in sample(chowPopulation, N): object 'chowPopulation' not found
```

Note that, despite the intervals being larger (we are dividing by $\sqrt{5}$ instead of $\sqrt{30}$) we see many more intervals not covering $\mu_X$. This is because the CLT is incorrectly telling us that the distribution of the `mean(hf)` is approximately normal when in fact it has fatter tail. This mistake affects us in the the calculation of `Q` which uses assumes a normal distribution and uses `qnorm`. The t-distribution might be more appropriate. All we have to do is re-run the above but change how we calculate `Q`: use `qt` instead of `qnorm`



```r
plot(mean(chowPopulation)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
abline(v=mean(chowPopulation))
```

```
## Error in mean(chowPopulation): object 'chowPopulation' not found
```

```r
##Q <- qnorm(1- 0.05/2) ##no longer normal so use:
Q <- qt(1- 0.05/2, df=4)
N<-5
for(i in 1:B){
  hf <- sample(chowPopulation,N)
  se=sd(hf)/sqrt(N)
  interval <- c(mean(hf)-Q*se, mean(hf)+Q*se )
  covered<-mean(chowPopulation)<= interval[2] & mean(chowPopulation)>=interval[1]
  color <- ifelse(covered,1,2)
  lines( interval, c(i,i),col=color)
}
```

```
## Error in sample(chowPopulation, N): object 'chowPopulation' not found
```

Note that now the intervals are made bigger. This is because the t-distribution has fatter tails and thus

```r
qt(1- 0.05/2, df=4) ##is bigger than 
```

```
## [1] 2.776445
```

```r
qnorm(1- 0.05/2)
```

```
## [1] 1.959964
```

which makes the intervals larger and thus cover $\mu_X$ more frequently. In fact, about 95% of the time.


# Connection between confidence intervals and p-values

We recommend that in practice confidence intervals be reported instead of p-values. If for some reason you are required to provide p-values, or required that you results are significant at the 0.05 of 0.01 leves, confidence intervals to provide this information. Note that if you use CLT and report $\bar{X} \pm 2 s_X/\sqrt{N}$ as a 95% confidence interval and this interval does not include 0. Because the interval does not include 0, this implies that either 
$\bar{X} - 2 s_X/\sqrt{N}  > 0$ or $\bar{X} + 2 s_X/\sqrt{N}  < 0$ which in turn implies that either $\sqrt{N}\bar{X}/s_X > 2$ or $\sqrt{N}\bar{X}/s_X < 2$ which implies the t-statistics is more extreme than 2 which implies the p-value must be smaller than 0.05. The same calculation can be made if we use the t-distribution instead of CLT. In summary, if a 95% or 99% confidence interval does not include 0 then the p-value must be smaller than 0.05 or 0.01 respectively.  

Note that the confidence interval is provided by the `t.test` function:




```r
library(downloader)
```

```
## Error in library(downloader): there is no package called 'downloader'
```

```r
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- tempfile()
download(url,destfile=filename)
```

```
## Error in eval(expr, envir, enclos): could not find function "download"
```

```r
dat <- read.csv(filename)
```

```
## Warning in file(file, "rt"): cannot open file
## '/var/folders/6d/d_8pbllx7318htlp5wv_rm580000gn/T//RtmpbVIFv7/file2ac22840bb2':
## No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
controlIndex <- which(dat$Diet=="chow")
```

```
## Error in which(dat$Diet == "chow"): object 'dat' not found
```

```r
treatmentIndex <- which(dat$Diet=="hf")
```

```
## Error in which(dat$Diet == "hf"): object 'dat' not found
```

```r
control <- dat[controlIndex,2]
```

```
## Error in eval(expr, envir, enclos): object 'dat' not found
```

```r
treatment <- dat[treatmentIndex,2]
```

```
## Error in eval(expr, envir, enclos): object 'dat' not found
```

```r
t.test(treatment,control)
```

```
## Error in t.test(treatment, control): object 'treatment' not found
```


In this case the 95% confidence interval does include 0 and we note that the p-value is larger than 0.05 as predicted. Note what happens if we change this to a 90% confidence interval:


```r
t.test(treatment,control,conf.level=0.9)
```

```
## Error in t.test(treatment, control, conf.level = 0.9): object 'treatment' not found
```

0 is no longer in the confidence interval and the p-value is smaller than 0.10.

