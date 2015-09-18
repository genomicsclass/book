---
layout: page
title: Smoothing
---



<a name="smoothing"></a>

## Smoothing 

Smoothing is a very powerful technique used all across data analysis. The general idea is to group data points that are expected to have similar expectations and [RAFA]

The following data are from measurements from replicated RNA. We consider the data used in an MA-plot ( $$Y$$ = log ratios and $$A$$ = averages) and take down-sample in a way that balances the number of points for different strata of $$A$$:


```r
##Following three packages are available from Bioconductor
library(Biobase)
library(SpikeIn)
```

```
## Error in library(SpikeIn): there is no package called 'SpikeIn'
```

```r
library(hgu95acdf)
```

```
## Error in library(hgu95acdf): there is no package called 'hgu95acdf'
```

```r
data(SpikeIn95)
```

```
## Warning in data(SpikeIn95): data set 'SpikeIn95' not found
```

```r
##Example with two columns
i=10;j=9

##remove the spiked in genes and take random sample
siNames<-colnames(pData(SpikeIn95))
```

```
## Error in colnames(pData(SpikeIn95)): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error in pData(SpikeIn95) : 
##   error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'SpikeIn95' not found
```

```r
ind <- which(!probeNames(SpikeIn95)%in%siNames)
```

```
## Error in match(x, table, nomatch = 0L): could not find function "probeNames"
```

```r
pms <- pm(SpikeIn95)[ ind ,c(i,j)]
```

```
## Error in eval(expr, envir, enclos): could not find function "pm"
```

```r
##pick a representative sample for A and order A
Y=log2(pms[,1])-log2(pms[,2])
```

```
## Error in eval(expr, envir, enclos): object 'pms' not found
```

```r
X=(log2(pms[,1])+log2(pms[,2]))/2
```

```
## Error in eval(expr, envir, enclos): object 'pms' not found
```

```r
set.seed(4)
ind <- tapply(seq(along=X),round(X*5),function(i)
  if(length(i)>20) return(sample(i,20)) else return(NULL))
```

```
## Error in tapply(seq(along = X), round(X * 5), function(i) if (length(i) > : error in evaluating the argument 'X' in selecting a method for function 'tapply': Error in seq(along = X) : object 'X' not found
```

```r
ind <- unlist(ind)
```

```
## Error in unlist(ind): error in evaluating the argument 'x' in selecting a method for function 'unlist': Error: object 'ind' not found
```

```r
X <- X[ind]
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```r
Y <- Y[ind]
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

```r
o <-order(X)
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```r
X <- X[o]
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```r
Y <- Y[o]
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

In the MA plot we see that $$Y$$ depends on $$X$$. This dependence must be a bias because these are based on replicates which means $$Y$$ should be 0 on average regardless of $$X$$. We want to predict $$f(x)=\mbox{E}(Y \mid X=x)$$ so that we can remove this bias.


```r
library(rafalib)
mypar()
plot(X,Y)
```

```
## Error in plot(X, Y): object 'X' not found
```

Linear regression does not capture the apparent curvature in $$f(x)$$:


```r
mypar()
plot(X,Y)
```

```
## Error in plot(X, Y): object 'X' not found
```

```r
fit <- lm(Y~X)
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

```r
points(X,Y,pch=21,bg=ifelse(Y>fit$fitted,1,3))
```

```
## Error in points(X, Y, pch = 21, bg = ifelse(Y > fit$fitted, 1, 3)): object 'X' not found
```

```r
abline(fit,col=2,lwd=4,lty=2)
```

```
## Error in abline(fit, col = 2, lwd = 4, lty = 2): object 'fit' not found
```

Note that the points above the fitted line (green) and those below (purple) are not evenly distributed.

## Bin Smoothing

Instead of fitting a line, let's go back to the idea of stratifying and computing the mean. This is referred to as _bin smoothing_. Now, if we stratify by $$x$$ , the general idea is that the underlying curve is "smooth" enough that in small bins it is approximately constant, which implies that all the $$Y$$ in that bin have the same expected value. For example, in the plot below we highlight points in a bin centered at 8.6 as well as the points of a bin centered at 12.1 if we use bins of size 1. We also show and the fitted mean values for the $$Y$$ in those bin (dashed lines):


```r
mypar()
centers <- seq(min(X),max(X),0.1)
```

```
## Error in seq(min(X), max(X), 0.1): object 'X' not found
```

```r
plot(X,Y,col="grey",pch=16)
```

```
## Error in plot(X, Y, col = "grey", pch = 16): object 'X' not found
```

```r
windowSize <- .5
i <- 25
center<-centers[i]
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
ind=which(X>center-windowSize & X<center+windowSize)
```

```
## Error in which(X > center - windowSize & X < center + windowSize): object 'X' not found
```

```r
fit<-mean(Y)
```

```
## Error in mean(Y): object 'Y' not found
```

```r
points(X[ind],Y[ind],bg=3,pch=21)
```

```
## Error in points(X[ind], Y[ind], bg = 3, pch = 21): object 'X' not found
```

```r
lines(c(min(X[ind]),max(X[ind])),c(fit,fit),col=2,lty=2,lwd=4)
```

```
## Error in lines(c(min(X[ind]), max(X[ind])), c(fit, fit), col = 2, lty = 2, : object 'X' not found
```

```r
i <- 60
center<-centers[i]
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
ind=which(X>center-windowSize & X<center+windowSize)
```

```
## Error in which(X > center - windowSize & X < center + windowSize): object 'X' not found
```

```r
fit<-mean(Y[ind])
```

```
## Error in mean(Y[ind]): object 'Y' not found
```

```r
points(X[ind],Y[ind],bg=3,pch=21)
```

```
## Error in points(X[ind], Y[ind], bg = 3, pch = 21): object 'X' not found
```

```r
lines(c(min(X[ind]),max(X[ind])),c(fit,fit),col=2,lty=2,lwd=4)
```

```
## Error in lines(c(min(X[ind]), max(X[ind])), c(fit, fit), col = 2, lty = 2, : object 'X' not found
```

By computing this mean for bins around every point we form an estimate of the underlying curve $$f(x)$$ :


```r
windowSize<-0.5
smooth<-rep(NA,length(centers))
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
mypar (4,3)
for(i in seq(along=centers)){
  center<-centers[i]
  ind=which(X>center-windowSize & X<center+windowSize)
  smooth[i]<-mean(Y[ind])
  if(i%%round(length(centers)/12)==1){ ##we show 12
    plot(X,Y,col="grey",pch=16)
    points(X[ind],Y[ind],bg=3,pch=21)
    lines(c(min(X[ind]),max(X[ind])),c(smooth[i],smooth[i]),col=2,lwd=2)
    lines(centers[1:i],smooth[1:i],col="black")
    points(centers[i],smooth[i],col="black",pch=16,cex=1.5)
  }
}
```

```
## Error in seq(along = centers): object 'centers' not found
```

The final result looks like this:


```r
mypar (1,1)
plot(X,Y,col="darkgrey",pch=16)
```

```
## Error in plot(X, Y, col = "darkgrey", pch = 16): object 'X' not found
```

```r
lines(centers,smooth,col="black",lwd=3)
```

```
## Error in lines(centers, smooth, col = "black", lwd = 3): object 'centers' not found
```


## Loess
 
Local weighted regression (loess) is similar to bin smoothing. The difference is that we approximate the local behavior with a line or a parabola. This permits us to expand the bin sizes as seen below:


```r
centers <- seq(min(X),max(X),0.1)
```

```
## Error in seq(min(X), max(X), 0.1): object 'X' not found
```

```r
mypar (1,1)
plot(X,Y,col="darkgrey",pch=16)
```

```
## Error in plot(X, Y, col = "darkgrey", pch = 16): object 'X' not found
```

```r
windowSize <- 1.25

i <- 25
center<-centers[i]
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
ind=which(X>center-windowSize & X<center+windowSize)
```

```
## Error in which(X > center - windowSize & X < center + windowSize): object 'X' not found
```

```r
fit<-lm(Y~X,subset=ind)
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

```r
points(X[ind],Y[ind],bg=3,pch=21)
```

```
## Error in points(X[ind], Y[ind], bg = 3, pch = 21): object 'X' not found
```

```r
a <- min(X[ind]);b <- max(X[ind])
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```r
lines(c(a,b),fit$coef[1]+fit$coef[2]*c(a,b),col=2,lty=2,lwd=3)
```

```
## Error in lines(c(a, b), fit$coef[1] + fit$coef[2] * c(a, b), col = 2, : object 'a' not found
```

```r
i <- 60
center<-centers[i]
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
ind=which(X>center-windowSize & X<center+windowSize)
```

```
## Error in which(X > center - windowSize & X < center + windowSize): object 'X' not found
```

```r
fit<-lm(Y~X,subset=ind)
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

```r
points(X[ind],Y[ind],bg=3,pch=21)
```

```
## Error in points(X[ind], Y[ind], bg = 3, pch = 21): object 'X' not found
```

```r
a <- min(X[ind]);b <- max(X[ind])
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```
## Error in eval(expr, envir, enclos): object 'X' not found
```

```r
lines(c(a,b),fit$coef[1]+fit$coef[2]*c(a,b),col=2,lty=2,lwd=3)
```

```
## Error in lines(c(a, b), fit$coef[1] + fit$coef[2] * c(a, b), col = 2, : object 'a' not found
```

Here are 12 steps of the process:

```r
mypar (4,3)
windowSize<-1.25
smooth<-rep(NA,length(centers))
```

```
## Error in eval(expr, envir, enclos): object 'centers' not found
```

```r
for(i in seq(along=centers)){
  center<-centers[i]
  ind=which(X>center-windowSize & X<center+windowSize)
  fit<-lm(Y~X,subset=ind)
  smooth[i]<-fit$coef[1]+fit$coef[2]*center

  if(i%%round(length(centers)/12)==1){ ##we show 12
    plot(X,Y,col="grey",pch=16)
    points(X[ind],Y[ind],bg=3,pch=21)
    a <- min(X[ind]);b <- max(X[ind])
    lines(c(a,b),fit$coef[1]+fit$coef[2]*c(a,b),col=2,lwd=2)
  
    lines(centers[1:i],smooth[1:i],col="black")
    points(centers[i],smooth[i],col="black",pch=16,cex=1.5)
  }
}
```

```
## Error in seq(along = centers): object 'centers' not found
```

This results in a smoother fit since we use larger sample sizes to estimate our local parameters:


```r
mypar (1,1)
plot(X,Y,col="darkgrey",pch=16)
```

```
## Error in plot(X, Y, col = "darkgrey", pch = 16): object 'X' not found
```

```r
lines(centers,smooth,col="black",lwd=3)
```

```
## Error in lines(centers, smooth, col = "black", lwd = 3): object 'centers' not found
```

The function `loess` performs this analysis for us:


```r
fit <- loess(Y~X, degree=1, span=1/3)
```

```
## Error in eval(expr, envir, enclos): object 'Y' not found
```

```r
newx <- seq(min(X),max(X),len=100) 
```

```
## Error in seq(min(X), max(X), len = 100): object 'X' not found
```

```r
smooth <- predict(fit,newdata=data.frame(X=newx))
```

```
## Error in predict(fit, newdata = data.frame(X = newx)): object 'fit' not found
```

```r
mypar ()
plot(X,Y,col="darkgrey",pch=16)
```

```
## Error in plot(X, Y, col = "darkgrey", pch = 16): object 'X' not found
```

```r
lines(newx,smooth,col="black",lwd=3)
```

```
## Error in lines(newx, smooth, col = "black", lwd = 3): object 'newx' not found
```

