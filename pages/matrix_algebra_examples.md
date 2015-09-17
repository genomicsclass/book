---
title: "Linear Algebra Examples "
author: "Rafa"
date: "February 18, 2015"
output: html_document
layout: page
---



## Examples 

The R markdown document for this section is available [here](https://github.com/genomicsclass/labs/tree/master/matrixalg/matrix_algebra_examples.Rmd).

Now we are ready to see how matrix algebra can be useful when analyzing data. We start with some simple examples and eventually arrive at the main one: how to write linear models with matrix algebra notation and solve the least squares problem.


#### The average

To compute the sample average and variance of our data, we use these formulas {$$}\bar{Y}=\frac{1}{N} Y_i{/$$} and {$$}\mbox{var}(Y)=\frac{1}{N} \sum_{i=1}^N (Y_i - \bar{Y})^2{/$$}. We can represent these with matrix multiplication. First, define this {$$}N \times 1{/$$} matrix made just of 1s:

{$$}
A=\begin{pmatrix}
1\\
1\\
\vdots\\
1
\end{pmatrix}
{/$$}

This implies that:

{$$}
\frac{1}{N}
\mathbf{A}^\top Y = \frac{1}{N}
\begin{pmatrix}1&1&,\dots&1\end{pmatrix}
\begin{pmatrix}
Y_1\\
Y_2\\
\vdots\\
Y_N
\end{pmatrix}=
\frac{1}{N} \sum_{i=1}^N Y_i
= \bar{Y}
{/$$}

Note that we are multiplying by the scalar {$$}1/N{/$$}. In R, we multiply matrix using `%*%`:


```r
library(UsingR)
y <- father.son$sheight
print(mean(y))
```

```
## [1] 68.68407
```

```r
N <- length(y)
Y<- matrix(y,N,1)
A <- matrix(1,N,1)
barY=t(A)%*%Y / N

print(barY)
```

```
##          [,1]
## [1,] 68.68407
```

#### The variance

As we will see later, multiplying the transpose of a matrix with another is very common in statistics. In fact, it is so common that there is a function in R:


```r
barY=crossprod(A,Y) / N
print(barY)
```

```
##          [,1]
## [1,] 68.68407
```

For the variance we note that if:

{$$}
\mathbf{r}\equiv \begin{pmatrix}
Y_1 - \bar{Y}\\
\vdots\\
Y_N - \bar{Y}
\end{pmatrix}, \,\,
\frac{1}{N} \mathbf{r}^\top\mathbf{r} = 
\frac{1}{N}\sum_{i=1}^N (Y_i - \bar{Y})^2
{/$$}

And in R if you only send one matrix into `crossprod`, it computes: {$$}r^\top r{/$$} so we can simply type:


```r
r <- y - barY
crossprod(r)/N
```

```
##          [,1]
## [1,] 7.915196
```

Which is almost equivalent to:

```r
var(y) 
```

```
## [1] 7.922545
```
The difference is due to the fact that `var` is for the sample estimate which divides by {$$}N-1{/$$}, so this:


```r
var(y) * (N-1) / N
```

```
## [1] 7.915196
```
gives us the same answer as our matrix multiplication example.

#### Linear models

Now we are ready to put all this to use. Let's start with Galton's example. If we define these matrices:
 
{$$}
\mathbf{Y} = \begin{pmatrix}
Y_1\\
Y_2\\
\vdots\\
Y_N
\end{pmatrix}
,
\mathbf{X} = \begin{pmatrix}
1&x_1\\
1&x_2\\
\vdots\\
1&x_N
\end{pmatrix}
,
\mathbf{\beta} = \begin{pmatrix}
\beta_0\\
\beta_1
\end{pmatrix} \mbox{ and }
\mathbf{\varepsilon} = \begin{pmatrix}
\varepsilon_1\\
\varepsilon_2\\
\vdots\\
\varepsilon_N
\end{pmatrix}
{/$$}



Then we can write the model:

{$$} 
Y_i = \beta_0 + \beta_1 x_i + \varepsilon, i=1,\dots,N 
{/$$}

as: 


{$$}
\,
\begin{pmatrix}
Y_1\\
Y_2\\
\vdots\\
Y_N
\end{pmatrix} = 
\begin{pmatrix}
1&x_1\\
1&x_2\\
\vdots\\
1&x_N
\end{pmatrix}
\begin{pmatrix}
\beta_0\\
\beta_1
\end{pmatrix} +
\begin{pmatrix}
\varepsilon_1\\
\varepsilon_2\\
\vdots\\
\varepsilon_N
\end{pmatrix}
{/$$}

or simply: 

{$$}
\mathbf{Y}=\mathbf{X}\boldsymbol{\beta}+\boldsymbol{\varepsilon}
{/$$}

which is a much simpler way to write it. 


The least squares equation becomes simpler as well since it is the following cross-product:

{$$}
(\mathbf{Y}-\mathbf{X}\boldsymbol{\beta})^\top
(\mathbf{Y}-\mathbf{X}\boldsymbol{\beta})
{/$$}

So now we are ready to determine which values of {$$}\beta{/$$} minimize the above. There are a series of rules that permit us to compute partial derivatives equations in matrix notation. By equating the derivative to 0 and solving for the {$$}\beta{/$$}, we will have our solution. The only one we need here tells us that the derivative of the above equation is:

{$$}
2 \mathbf{X}^\top (\mathbf{Y} - \mathbf{X} \boldsymbol{\hat{\beta}})=0
{/$$}

{$$}
\mathbf{X}^\top \mathbf{X} \boldsymbol{\hat{\beta}} = \mathbf{X}^\top \mathbf{Y}   
{/$$}


{$$}
\boldsymbol{\hat{\beta}} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y}   
{/$$}

and we have our solution. We usually put a hat on the {$$}\beta{/$$} that solves this, {$$}\hat{\beta}{/$$} as it is an estimate of the "real" {$$}\beta{/$$} that generated the data.

Remember that the least squares are like a square (multiply something by itself) and that this formula is similar to the derivative of {$$}f(x)^2{/$$} being {$$}2f(x)f\prime (x){/$$}. 

Let's see how it works in R:


```r
library(UsingR)
x=father.son$fheight
y=father.son$sheight
X <- cbind(1,x)
betahat <- solve( t(X) %*% X ) %*% t(X) %*% y
###or
betahat <- solve( crossprod(X) ) %*% crossprod( X, y )
```


Now we can see the results of this by computing the estimated {$$}\hat{\beta}_0+\hat{\beta}_1 x{/$$} for any value of {$$}x{/$$}:


```r
newx <- seq(min(x),max(x),len=100)
X <- cbind(1,newx)
fitted <- X%*%betahat
plot(x,y,xlab="Father's height",ylab="Son's height")
lines(newx,fitted,col=2)
```

![Galton's data with fitted regression line.](images/R/matrix_algebra_examples-tmp-galton_regression_line-1.png) 

This {$$}\hat{\boldsymbol{\beta}}=(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y}{/$$} is one of the most widely used results in data analysis. One of the advantages of this approach is that we can use it in many different situations.  For example, in our falling object problem: 
 

```r
set.seed(1)
g <- 9.8 #meters per second
n <- 25
tt <- seq(0,3.4,len=n) #time in secs, t is a base function
d <- 56.67  - 0.5*g*tt^2 + rnorm(n,sd=1)
```

Note that we are using almost the same exact code:



```r
X <- cbind(1,tt,tt^2)
y <- d
betahat <- solve(crossprod(X))%*%crossprod(X,y)
newtt <- seq(min(tt),max(tt),len=100)
X <- cbind(1,newtt,newtt^2)
fitted <- X%*%betahat
plot(tt,y,xlab="Time",ylab="Height")
lines(newtt,fitted,col=2)
```

![Fitted parabola to simulated data for distance travelled versus time of falling object measured with error.](images/R/matrix_algebra_examples-tmp-gravity_with_fitted_parabola-1.png) 

And the resulting estimates are what we expect:


```r
betahat
```

```
##          [,1]
##    56.5317368
## tt  0.5013565
##    -5.0386455
```

The Tower of Pisa is about 56 meters high, there is no initial velocity and half the constant of gravity is 9.8/2=4.9.

#### The `lm` Function
R has a very convenient function that fits these models. We will learn more about this function later, but here is a preview:


```r
X <- cbind(tt,tt^2)
fit=lm(y~X)
summary(fit)
```

```
## 
## Call:
## lm(formula = y ~ X)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.5295 -0.4882  0.2537  0.6560  1.5455 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  56.5317     0.5451 103.701   <2e-16 ***
## Xtt           0.5014     0.7426   0.675    0.507    
## X            -5.0386     0.2110 -23.884   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9822 on 22 degrees of freedom
## Multiple R-squared:  0.9973,	Adjusted R-squared:  0.997 
## F-statistic:  4025 on 2 and 22 DF,  p-value: < 2.2e-16
```

We obtain the same values as above.


### Summary

We have shown how to write linear models using linear algebra. We are going to do this for several examples, many of which are related to designed experiments. We showed how to obtain least squares estimates. Keep in mind, however, that because {$$}y{/$$} is a random variable, these estimates are random as well. In a later section, we will learn how to compute standard error for these estimates and use this to perform inference.




