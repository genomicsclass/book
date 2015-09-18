---
layout: page
title: Modeling Batch Effects with Factor Analysis
---



##  Modeling Batch Effects with Factor Analysis

To illustrate how we can adjust for batch effects using statistical methods, we continue to use our subset example with sex as the outcome of interest, and with months purposely somewhat confounded with sex. 


```r
library(GSE5859Subset)
```

```
## Error in library(GSE5859Subset): there is no package called 'GSE5859Subset'
```

```r
data(GSE5859Subset)
```

```
## Warning in data(GSE5859Subset): data set 'GSE5859Subset' not found
```

Below is an image we showed earlier with a subset of genes showing both the sex effect and the time effects, along with sample to sample correlations (computed on all genes) showing the complex structure of the data:



```r
library(rafalib)
library(RColorBrewer)
library(genefilter)


sex <- sampleInfo$group
```

```
## Error in eval(expr, envir, enclos): object 'sampleInfo' not found
```

```r
batch <- factor(format(sampleInfo$date,"%m"))
```

```
## Error in format(sampleInfo$date, "%m"): object 'sampleInfo' not found
```

```r
chr <- geneAnnotation$CHR
```

```
## Error in eval(expr, envir, enclos): object 'geneAnnotation' not found
```

```r
tt<-rowttests(geneExpression,batch)
```

```
## Error in rowttests(geneExpression, batch): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'geneExpression' not found
```

```r
ind1 <- which(chr=="chrY") #real differences
```

```
## Error in which(chr == "chrY"): object 'chr' not found
```

```r
ind2 <- setdiff(c(order(tt$dm)[1:25],order(-tt$dm)[1:25]),ind1)
```

```
## Error in order(tt$dm): object 'tt' not found
```

```r
set.seed(1)
ind0 <- setdiff(sample(seq(along=tt$dm),50),c(ind2,ind1))
```

```
## Error in seq(along = tt$dm): object 'tt' not found
```

```r
geneindex<-c(ind2,ind0,ind1)
```

```
## Error in eval(expr, envir, enclos): object 'ind2' not found
```

```r
mat<-geneExpression[geneindex,]
```

```
## Error in eval(expr, envir, enclos): object 'geneExpression' not found
```

```r
mat <- mat -rowMeans(mat)
```

```
## Error in eval(expr, envir, enclos): object 'mat' not found
```

```r
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

mypar(1,2)
image(t(mat),xaxt="n",yaxt="n",col=icolors)
```

```
## Error in t(mat): object 'mat' not found
```

```r
y <- geneExpression - rowMeans(geneExpression)
```

```
## Error in eval(expr, envir, enclos): object 'geneExpression' not found
```

```r
image(1:ncol(y),1:ncol(y),cor(y),col=icolors,zlim=c(-1,1),
       xaxt="n",xlab="",yaxt="n",ylab="")
```

```
## Error in ncol(y): object 'y' not found
```

```r
axis(2,1:ncol(y),sex,las=2)
```

```
## Error in ncol(y): object 'y' not found
```

```r
axis(1,1:ncol(y),sex,las=2)
```

```
## Error in ncol(y): object 'y' not found
```



We have seen how approaches that assume month explains the batch and use linear models perform relatively well. However, there was still room for improvement. This is most likely due to the fact that month is only a surrogate for some other variable that actually induces structure or between sample correlation.

#### What is a batch?

Here is a plot of dates for each sample, with color representing month:


```r
times <-sampleInfo$date 
```

```
## Error in eval(expr, envir, enclos): object 'sampleInfo' not found
```

```r
mypar(1,1)
o=order(times)
```

```
## Error in order(times): object 'times' not found
```

```r
plot(times[o],pch=21,bg=as.numeric(batch)[o],ylab="date")
```

```
## Error in plot(times[o], pch = 21, bg = as.numeric(batch)[o], ylab = "date"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'times' not found
```

```r
o=order(times)
```

```
## Error in order(times): object 'times' not found
```

```r
plot(times[o],pch=21,bg=as.numeric(batch)[o],ylab="date")
```

```
## Error in plot(times[o], pch = 21, bg = as.numeric(batch)[o], ylab = "date"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'times' not found
```
There is more than one day per month. Could day have an effect as well?


#### PCA

Here is a plot of the first principal component ordered by date:

```r
s <- svd(y)
```

```
## Error in as.matrix(x): object 'y' not found
```

```r
mypar(1,1)
o<-order(times)
```

```
## Error in order(times): object 'times' not found
```

```r
cols <- as.numeric( batch)
```

```
## Error in eval(expr, envir, enclos): object 'batch' not found
```

```r
plot(s$v[o,1],pch=21,cex=1.25,bg=cols[o],ylab="First PC",xaxt="n",xlab="")
```

```
## Error in plot(s$v[o, 1], pch = 21, cex = 1.25, bg = cols[o], ylab = "First PC", : error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 's' not found
```

```r
legend("topleft",c("Month 1","Month 2"),col=1:2,pch=16,box.lwd=0)
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```

Day seems to be highly correlated and explained a high percentage of the variability:


```r
mypar(1,1)
plot(s$d^2/sum(s$d^2),ylab="% variance explained",xlab="Principal component")
```

```
## Error in plot(s$d^2/sum(s$d^2), ylab = "% variance explained", xlab = "Principal component"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 's' not found
```

In fact, the first six or so PC seem to be at least partially driven by date:

```r
mypar(3,4)
for(i in 1:12){
  days <- gsub("2005-","",times)  
  boxplot(split(s$v[,i],gsub("2005-","",days)))
}
```

```
## Error in gsub("2005-", "", times): object 'times' not found
```


What happens if we simply remove the top six PC from the data and then perform a t-test? 


```r
D <- s$d; D[1:4]<-0 #take out first 2
```

```
## Error in eval(expr, envir, enclos): object 's' not found
```

```
## Error in D[1:4] <- 0: object of type 'closure' is not subsettable
```

```r
cleandat <- sweep(s$u,2,D,"*")%*%t(s$v)
```

```
## Error in sweep(s$u, 2, D, "*"): object 's' not found
```

```r
res <-rowttests(cleandat,factor(sex))
```

```
## Error in rowttests(cleandat, factor(sex)): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'cleandat' not found
```

This does remove the batch effect, but it seems we have also removed much of the biological effect we are interested in. In fact, no genes have q-value <0.1 anymore.



```r
mypar(1,2)
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))
```

```
## Error in hist(res$p.value[which(!chr %in% c("chrX", "chrY"))], main = "", : object 'res' not found
```

```r
plot(res$dm,-log10(res$p.value))
```

```
## Error in plot(res$dm, -log10(res$p.value)): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'res' not found
```

```r
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
```

```
## Error in points(res$dm[which(chr == "chrX")], -log10(res$p.value[which(chr == : object 'res' not found
```

```r
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
```

```
## Error in points(res$dm[which(chr == "chrY")], -log10(res$p.value[which(chr == : object 'res' not found
```

```r
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```

<a name="sva"></a>
#### Surrogate Variable Analysis

An alternative is to fit models with both the covariate of interest, as well as those believed to be batches. An example of an approach that does this is Surrogate Variable Analysis (SVA).

The basic idea of SVA is to first estimate the factors, but taking care not to include the outcome of interest. To do this, an interactive approach is used in which each row is given a weight. These weights are then used in the SVD calculation with higher weights given to rows not associated with the outcome of interest and associated with batches. Below is a demonstration of two iterations. The three images are the data (for a subset of genes), the weights, and the estimated first factor.



```r
library(sva)
```

```
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-7. For overview type 'help("mgcv-package")'.
```

```r
library(limma)
mod <- model.matrix(~sex)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
cind <- order( as.Date(sampleInfo$date) )
```

```
## Error in as.Date(sampleInfo$date): object 'sampleInfo' not found
```

```r
dates <- gsub("2005-","",sampleInfo$date)
```

```
## Error in gsub("2005-", "", sampleInfo$date): object 'sampleInfo' not found
```

```r
weights=rep(1,nrow(y))
```

```
## Error in nrow(y): object 'y' not found
```

```r
par(mar = c(4.1, 2.1, 3.5, 2.1), 
    mgp = c(1.5, 0.5, 0))
layout(matrix(c(1:6),nrow=2,byrow=TRUE),widths=c(5,1.5,5))
for(b in 1:2){
  image(1:ncol(mat),1:nrow(mat),t(mat[,cind]*weights[geneindex]),xaxt="n",yaxt="n",col=icolors,xlab="",ylab="")
  axis(side=1,seq(along=dates),dates[cind],las=2)
  abline(v=12.5)
  
  svafit <- sva(y,mod,B=b,n.sv=5)
  weights = svafit$pprob.gam*(1-svafit$pprob.b)
  
  surrogate <- svd( y*weights)$v[,1]#Weighted SVD
  
  image(matrix(weights[geneindex],nrow=1),xaxt="n",yaxt="n",col=brewer.pal(9,"Blues"))
  plot(surrogate[cind],bg=sex[cind]+1,pch=21,xlab="",xaxt="n",ylab="Surrogate variable",ylim=c(-.5,.5),cex=1.5)
  axis(side=1,seq(along=dates),dates[cind],las=2)
  abline(v=12.5)
  text(1,0.5,"June")
  text(13.5,0.5,"Oct")
  legend("bottomright",c("0","1"),col=c(1,2),pch=16)
}
```

```
## Error in ncol(mat): object 'mat' not found
```


The above is an illustration of the algorithm. To actually run SVA, we follow the code. In this case, SVA picks the number of surrogate values or factors for us.



```r
library(limma)
svafit <- sva(geneExpression,mod)
```

```
## Error in ncol(dat): object 'geneExpression' not found
```

```r
svaX<-model.matrix(~sex+svafit$sv)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
lmfit <- lmFit(geneExpression,svaX)
```

```
## Error in is(object, "list"): object 'geneExpression' not found
```

```r
tt<- lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)
```

```
## Error in eval(expr, envir, enclos): object 'lmfit' not found
```

There is an observable improvement over all other approaches:


```r
res <- data.frame(dm= -lmfit$coef[,2],
                  p.value=2*(1-pt(abs(tt),lmfit$df.residual[1]) ) )
```

```
## Error in data.frame(dm = -lmfit$coef[, 2], p.value = 2 * (1 - pt(abs(tt), : object 'lmfit' not found
```

```r
mypar(1,2)
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))
```

```
## Error in hist(res$p.value[which(!chr %in% c("chrX", "chrY"))], main = "", : object 'res' not found
```

```r
plot(res$dm,-log10(res$p.value))
```

```
## Error in plot(res$dm, -log10(res$p.value)): error in evaluating the argument 'x' in selecting a method for function 'plot': Error: object 'res' not found
```

```r
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
```

```
## Error in points(res$dm[which(chr == "chrX")], -log10(res$p.value[which(chr == : object 'res' not found
```

```r
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
```

```
## Error in points(res$dm[which(chr == "chrY")], -log10(res$p.value[which(chr == : object 'res' not found
```

```r
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```


And here is a decompose of the data into sex effects, surrogate variables, and independent noise:


```r
Batch<- lmfit$coef[geneindex,3:7]%*%t(svaX[,3:7])
```

```
## Error in eval(expr, envir, enclos): object 'lmfit' not found
```

```r
Signal<-lmfit$coef[geneindex,1:2]%*%t(svaX[,1:2])
```

```
## Error in eval(expr, envir, enclos): object 'lmfit' not found
```

```r
error <- geneExpression[geneindex,]-Signal-Batch
```

```
## Error in eval(expr, envir, enclos): object 'geneExpression' not found
```

```r
##demean for plot
Signal <-Signal-rowMeans(Signal)
```

```
## Error in eval(expr, envir, enclos): object 'Signal' not found
```

```r
mat <- geneExpression[geneindex,]-rowMeans(geneExpression[geneindex,])
```

```
## Error in eval(expr, envir, enclos): object 'geneExpression' not found
```

```r
mypar(1,4,mar = c(2.75, 4.5, 2.6, 1.1))
image(t(mat),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
```

```
## Error in t(mat): object 'mat' not found
```

```r
image(t(Signal),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
```

```
## Error in t(Signal): object 'Signal' not found
```

```r
image(t(Batch),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
```

```
## Error in t(Batch): object 'Batch' not found
```

```r
image(t(error),col=icolors,zlim=c(-5,5),xaxt="n",yaxt="n")
```

```
## Error in t(error): object 'error' not found
```


