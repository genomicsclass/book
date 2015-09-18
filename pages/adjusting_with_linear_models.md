---
layout: page
title: Adjusting for Batch Effects with Linear Models
---




## Data Example

To illustrate how we can adjust for batch effects using statistical methods, we will create a data example in which the outcome of interest is confounded with batch, but not completely. We will also select an outcome for which we have an expectation of what genes should be diferentially expressed. Namely, we make sex the outcome of interest and expect genes on the Y chromosome to be diferentially expressed. We may also see genes from the X chromosome as differentially expressed since some escape X inactivation. The example dataset is below.


```r
##available from course github repository
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

To illustrate the confounding, we will pick some genes to show in a heatmap plot. We pick all Y chromosome genes, some genes that we see correlate with batch, and then some randomly selected genes.


```r
library(rafalib)
library(RColorBrewer)
library(genefilter)

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
ind1 <- which(chr=="chrY") ##real differences
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

The object `mat` contains the subset of the data we will show in the plots. By looking at an image we see the Y chromosome genes as well as those most affected by month:


```r
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
mypar(1,1)
image(t(mat),xaxt="n",yaxt="n",col=icolors)
```

```
## Error in t(mat): object 'mat' not found
```

In what follows, we will imitate the typical analysis we would do in practice. We will act as if we don't know which genes are supposed to be differentially expressed between males and females. 

#### Exploratory Data Analysis for Evaluation

Another reason we are using the dataset described above for illustrating different approaches is because we actually have a reasonable idea of what to expect. Autosomal (not on chrX or chrY) genes on the list are likely false positives and chrY are likely true positives. ChrX genes could go either way. This gives us the opportunity to compare different procedures. Since in practice we rarely know the "truth", these evaluations are not possible. Simulations are therefore commonly used for evaluation purposes: we know the truth because we construct the data. But simulations are at risk of not capturing all the nuances of real experimental data. This dataset is an experimental dataset. 

In the next sections we will use the histogram p-values to evaluate the specificity (low false positive rates) of the batch adjustment procedures presented here. Because the autosomal genes are not expected to be differentially expressed, we should see a a flat p-value histogram. To evaluate sensitivity (low false negative rates), we will report the number of the reported genes on chrX and chrY. Below are the results for when we don't adjust and report genes with q-values smaller than 0.1. We also include a volcano plot with a horizontal dashed line separating the genes called significant from those that are not, and colors used to highlight chrX and chrY genes.


```r
library(qvalue)
```

```
## Error in library(qvalue): there is no package called 'qvalue'
```

```r
res <- rowttests(geneExpression,as.factor( sampleInfo$group ))
```

```
## Error in rowttests(geneExpression, as.factor(sampleInfo$group)): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'geneExpression' not found
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

```r
qvals <- qvalue(res$p.value)$qvalue
```

```
## Error in eval(expr, envir, enclos): could not find function "qvalue"
```

```r
index <- which(qvals<0.1)
```

```
## Error in which(qvals < 0.1): object 'qvals' not found
```

```r
abline(h=-log10(max(res$p.value[index])))
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): object 'res' not found
```

```r
cat("Total genes with q-value < 0.1:",length(index))
```

```
## Error in cat("Total genes with q-value < 0.1:", length(index)): object 'index' not found
```

```r
cat("Number of selected genes on chrY:", sum(chr[index]=="chrY",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrY:", sum(chr[index] == "chrY", : object 'chr' not found
```

```r
cat("Number of selected genes on chrX:", sum(chr[index]=="chrX",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrX:", sum(chr[index] == "chrX", : object 'chr' not found
```

The histogram is not flat. Instead, low p-values are over-represented. More than half of the genes on the final list are autosomal.

## Adjusting for Batch Effects with Linear Models

We have already observed that processing date has an effect on gene expression.  We will therefore try to _adjust_ for this by including it in a model.  When we perform a t-test comparing the two groups, it is equivalent to fitting the following linear model:

$$Y_{ij} = \alpha_j + x_i \beta_{j} + \varepsilon_{ij}$$

to each gene $$j$$ with $$x_i=1$$ if subject $$i$$ is female and 0 otherwise. Note that $$\beta_{j}$$ represents the estimated difference for gene $$j$$ and $$\varepsilon_{ij}$$ represents the within group variation. So what is the problem?

The theory we described in the linear models chapter assumes that the error terms are independent. We know that this is not the case for all genes because we know the error terms from October will be more alike to each other than the June error terms. We can _adjust_ for this by including a term that models this effect:


$$Y_{ij} = \alpha_j + x_i \beta_{j} + z_i \gamma_j+\varepsilon_{ij}.$$

Here $$z_i=1$$ if sample $$i$$ was processed in October and 0 otherwise and $$\gamma_j$$ is the month effect for gene $$j$$. This an example of how linear models give us much more flexibility than procedures such as the t-test.

We construct a model matrix that includes batch.

```r
sex <- sampleInfo$group
```

```
## Error in eval(expr, envir, enclos): object 'sampleInfo' not found
```

```r
X <- model.matrix(~sex+batch)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

Now we can fit a model for each gene. For example, note the difference between the original model and one that has been adjusted for batch:


```r
j <- 7635
y <- geneExpression[j,]
```

```
## Error in eval(expr, envir, enclos): object 'geneExpression' not found
```

```r
X0 <- model.matrix(~sex) 
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
fit <- lm(y~X0-1)
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
summary(fit)$coef
```

```
## Error in summary(fit): object 'fit' not found
```

```r
X <- model.matrix(~sex+batch)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
fit <- lm(y~X)
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
summary(fit)$coef
```

```
## Error in summary(fit): object 'fit' not found
```

We then fit this new model for each gene. For instance, we can use `sapply` to recover the estimated coefficient and p-value in the following way:


```r
res <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
} ) )
```

```
## Error in nrow(geneExpression): object 'geneExpression' not found
```

```r
##turn into data.frame so we can use the same code for plots as above
res <- data.frame(res)
```

```
## Error in data.frame(res): object 'res' not found
```

```r
names(res) <- c("dm","p.value")
```

```
## Error in names(res) <- c("dm", "p.value"): object 'res' not found
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

```r
qvals <- qvalue(res$p.value)$qvalue
```

```
## Error in eval(expr, envir, enclos): could not find function "qvalue"
```

```r
index <- which(qvals<0.1)
```

```
## Error in which(qvals < 0.1): object 'qvals' not found
```

```r
abline(h=-log10(max(res$p.value[index])))
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): object 'res' not found
```

```r
cat("Total genes with q-value < 0.1:",length(index))
```

```
## Error in cat("Total genes with q-value < 0.1:", length(index)): object 'index' not found
```

```r
cat("Number of selected genes on chrY:", sum(chr[index]=="chrY",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrY:", sum(chr[index] == "chrY", : object 'chr' not found
```

```r
cat("Number of selected genes on chrX:", sum(chr[index]=="chrX",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrX:", sum(chr[index] == "chrX", : object 'chr' not found
```

There is a great improvement in specificity (less false positives) without much loss in sensitivity (we still find many chrY genes). However, we still see some bias in the histogram. In the following sections we will see that month does not perfectly account for the batch effect and that better estimates are possible.


## A Note on Computing Efficiency

In the code above, the design matrix does not change within the iterations we are computing $$(X^\top X)^{-1}$$ repeatedly and applying to each gene. Instead we can perform this calculation in one matrix algebra calculation by computing it once and then obtaining all the betas by multiplying $$(X^\top X)^{-1}X^\top Y$$ with the columns of $$Y$$ representing genes in this case. The `limma` package has an implementation of this idea (using the QR decomposition). Notice how much faster this is:


```r
library(limma)
X <- model.matrix(~sex+batch)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
fit <- lmFit(geneExpression,X)
```

```
## Error in is(object, "list"): object 'geneExpression' not found
```

The estimated regression coefficients for each gene are obtained like this:

```r
dim( fit$coef)
```

```
## Error in eval(expr, envir, enclos): object 'fit' not found
```
We have one estimate for each gene. To obtain p-values for one of these, we have to construct the ratios:


```r
k <- 2 ##second coef
ses <- fit$stdev.unscaled[,k]*fit$sigma
```

```
## Error in eval(expr, envir, enclos): object 'fit' not found
```

```r
ttest <- fit$coef[,k]/ses
```

```
## Error in eval(expr, envir, enclos): object 'fit' not found
```

```r
pvals <- 2*pt(-abs(ttest),fit$df)
```

```
## Error in abs(ttest): non-numeric argument to mathematical function
```

#### Combat

 Combat [NEED CITATION] is a popular method and is based on using linear models to adjust for batch effects. It fits a hierarchical model (we will learn about these in the next section) to estimate and remove row specific batch effects. Combat uses a modular approach. In a first step, what is considered to be a batch effect is removed:


```r
library(sva) #available from Bioconductor
```

```
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-7. For overview type 'help("mgcv-package")'.
```

```r
mod <- model.matrix(~sex)
```

```
## Error in eval(expr, envir, enclos): object 'sex' not found
```

```r
cleandat <- ComBat(geneExpression,batch,mod)
```

```
## Error in ComBat(geneExpression, batch, mod): object 'batch' not found
```

Then the results can be used to fit a model with our variable of interest:



```r
res<-genefilter::rowttests(cleandat,factor(sex))
```

```
## Error in genefilter::rowttests(cleandat, factor(sex)): error in evaluating the argument 'x' in selecting a method for function 'rowttests': Error: object 'cleandat' not found
```

In this case, the results are less specific than what we obtain by fitting the simple linear model:


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

```r
qvals <- qvalue(res$p.value)$qvalue
```

```
## Error in eval(expr, envir, enclos): could not find function "qvalue"
```

```r
index <- which(qvals<0.1)
```

```
## Error in which(qvals < 0.1): object 'qvals' not found
```

```r
abline(h=-log10(max(res$p.value[index])))
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): object 'res' not found
```

```r
cat("Total genes with q-value < 0.1:",length(index))
```

```
## Error in cat("Total genes with q-value < 0.1:", length(index)): object 'index' not found
```

```r
cat("Number of selected genes on chrY:", sum(chr[index]=="chrY",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrY:", sum(chr[index] == "chrY", : object 'chr' not found
```

```r
cat("Number of selected genes on chrX:", sum(chr[index]=="chrX",na.rm=TRUE))
```

```
## Error in cat("Number of selected genes on chrX:", sum(chr[index] == "chrX", : object 'chr' not found
```

