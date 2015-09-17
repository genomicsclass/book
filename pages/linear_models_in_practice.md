---
title: Linear models in practice
layout: page
---



## Linear Models In Practice

The R markdown document for this section is available [here](https://github.com/genomicsclass/labs/tree/master/linear/linear_models_in_practice.Rmd).

We will demonstrate how to analyze the high fat diet data using linear models instead of directly applying a t-test. We will demonstrate how, ultimately, these two approaches are equivalent. 

We start by reading in the data and creating a quick stripchart:




```r
set.seed(1) #same jitter in stripchart
dat <- read.csv("femaleMiceWeights.csv") ##previously downloaded
stripchart(dat$Bodyweight ~ dat$Diet, vertical=TRUE, method="jitter",
           main="Bodyweight over Diet")
```

![Mice bodyweights stratified by diet.](images/R/linear_models_in_practice-tmp-bodyweight_by_diet_stripchart-1.png) 

We can see that the high fat diet group appears to have higher weights on average, although there is overlap between the two samples.

#### A linear model with one variable

For demonstration purposes, we will build the design matrix {$$}\mathbf{X}{/$$} using the formula `~ Diet`. The group with the 1's in the second column is determined by the level of `Diet` which comes second; that is, the non-reference level. 


```r
levels(dat$Diet)
```

```
## [1] "chow" "hf"
```

```r
X <- model.matrix(~ Diet, data=dat)
X
```

```
##    (Intercept) Diethf
## 1            1      0
## 2            1      0
## 3            1      0
## 4            1      0
## 5            1      0
## 6            1      0
## 7            1      0
## 8            1      0
## 9            1      0
## 10           1      0
## 11           1      0
## 12           1      0
## 13           1      1
## 14           1      1
## 15           1      1
## 16           1      1
## 17           1      1
## 18           1      1
## 19           1      1
## 20           1      1
## 21           1      1
## 22           1      1
## 23           1      1
## 24           1      1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$Diet
## [1] "contr.treatment"
```

```r
colnames(X)
```

```
## [1] "(Intercept)" "Diethf"
```

```r
dat$Diet <- relevel(dat$Diet, ref="hf")
model.matrix(~ Diet, data=dat)
```

```
##    (Intercept) Dietchow
## 1            1        1
## 2            1        1
## 3            1        1
## 4            1        1
## 5            1        1
## 6            1        1
## 7            1        1
## 8            1        1
## 9            1        1
## 10           1        1
## 11           1        1
## 12           1        1
## 13           1        0
## 14           1        0
## 15           1        0
## 16           1        0
## 17           1        0
## 18           1        0
## 19           1        0
## 20           1        0
## 21           1        0
## 22           1        0
## 23           1        0
## 24           1        0
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$Diet
## [1] "contr.treatment"
```

After trying out the `relevel` function, we finally reset `chow` as the reference level because we want the comparison to be {$$}hf - chow{/$$}:


```r
dat$Diet <- relevel(dat$Diet, ref="chow")
```

## The Mathematics Behind lm()

The R markdown document for this section is available [here](https://github.com/genomicsclass/labs/tree/master/linear/linear_models_in_practice.Rmd).

Before we use our shortcut for running linear models, `lm`, we want to review what will happen internally. Inside of `lm`, we will form the design matrix {$$}\mathbf{X}{/$$}, and calculate the {$$}\boldsymbol{\beta}{/$$} which minimizes the sum of squares, as described in a previous lecture. The formula for this solution is:

{$$} \hat{\boldsymbol{\beta}} = (\mathbf{X}^t \mathbf{X})^{-1} \mathbf{X}^t \mathbf{Y} {/$$}

We can calculate this in R using our matrix multiplication operator `%*%`, the inverse function `solve` and the transpose function `t`.



```r
Y <- dat$Bodyweight
X <- model.matrix(~ Diet, data=dat)
solve(t(X) %*% X) %*% t(X) %*% Y
```

```
##                  [,1]
## (Intercept) 23.813333
## Diethf       3.020833
```

These coefficients are the average of the control group and the difference of the averages:



```r
s <- split(dat$Bodyweight, dat$Diet)
mean(s[["chow"]])
```

```
## [1] 23.81333
```

```r
mean(s[["hf"]]) - mean(s[["chow"]])
```

```
## [1] 3.020833
```

Finally, we use our shortcut, `lm`, to run the linear model:


```r
fit <- lm(Bodyweight ~ Diet, data=dat)
summary(fit)
```

```
## 
## Call:
## lm(formula = Bodyweight ~ Diet, data = dat)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -6.1042 -2.4358 -0.4138  2.8335  7.1858 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   23.813      1.039  22.912   <2e-16 ***
## Diethf         3.021      1.470   2.055   0.0519 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 3.6 on 22 degrees of freedom
## Multiple R-squared:  0.1611,	Adjusted R-squared:  0.1229 
## F-statistic: 4.224 on 1 and 22 DF,  p-value: 0.05192
```

```r
(coefs <- coef(fit))
```

```
## (Intercept)      Diethf 
##   23.813333    3.020833
```

#### Examining the coefficients

The following large and clunky piece of code allows us to visualize the meaning of the coefficients with colored arrows:


```r
stripchart(dat$Bodyweight ~ dat$Diet, vertical=TRUE, method="jitter",
           main="Bodyweight over Diet", ylim=c(0,40), xlim=c(0,3))
a <- -0.25
lgth <- .1
library(RColorBrewer)
cols <- brewer.pal(3,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
abline(h=coefs[1]+coefs[2],col=cols[2])
legend("right",names(coefs),fill=cols,cex=.75,bg="white")
```

![Estimated linear model coefficients for bodyweight data illustrated with arrows.](images/R/linear_models_in_practice-tmp-parameter_estimate_illustration-1.png) 

## Comparing Simple Two Group lm to a t-test

The R markdown document for this section is available [here](https://github.com/genomicsclass/labs/tree/master/linear/linear_models_in_practice.Rmd).

To make a connection with material presented earlier, this simple linear model is actually giving us the same result (the t-statistic and p-value) for the difference as a specific kind of t-test. This is the t-test between two groups with the assumption that both groups have the same variance. This was encoded into our linear model when we assumed that the errors {$$}\boldsymbol{\varepsilon}{/$$} were all equally distributed.

Though, in this case, the linear model is equivalent to a t-test, we will soon explore more complicated designs, where the linear model is a useful extension.

Our `lm` coefficients were:


```r
summary(fit)$coefficients
```

```
##              Estimate Std. Error   t value     Pr(>|t|)
## (Intercept) 23.813333   1.039353 22.911684 7.642256e-17
## Diethf       3.020833   1.469867  2.055174 5.192480e-02
```

And the t-statistic of the t-test is the same, with a flipped sign:


```r
(ttest <- t.test(s[["chow"]], s[["hf"]], var.equal=TRUE))
```

```
## 
## 	Two Sample t-test
## 
## data:  s[["chow"]] and s[["hf"]]
## t = -2.0552, df = 22, p-value = 0.05192
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -6.06915183  0.02748516
## sample estimates:
## mean of x mean of y 
##  23.81333  26.83417
```

```r
summary(fit)$coefficients[2,3]
```

```
## [1] 2.055174
```

```r
ttest$statistic
```

```
##         t 
## -2.055174
```

If we put the high fat group first, we get the same sign as the linear model:


```r
t.test(s[["hf"]], s[["chow"]], var.equal=TRUE)$statistic
```

```
##        t 
## 2.055174
```

