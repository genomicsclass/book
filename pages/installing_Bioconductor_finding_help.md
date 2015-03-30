---
layout: page
title: Installing Bioconductor and finding help
---

## Installing Bioconductor

In order to install Bioconductor, copy the following two lines into your R console.


```r
source("http://bioconductor.org/biocLite.R")
biocLite()
```

This will install the core Bioconductor packages. Further packages can be installed using the `biocLite` function and specifying a character vector of which packages to install. For example, to install the "affy" and "genefilter" libraries you would type:


```r
biocLite(c("genefilter","geneplotter"))
```

Remember, you still need to load the library, e.g., `library(genefilter)`, if you want to use the package.

More information on installing and updating Bioconductor packages can be found at:

http://bioconductor.org/install/

## Finding help 

There are many ways to find help directly from within R. Typically, every function will have its own manual page which is accessible by typing a question mark ?, followed by the function name and hitting return.


```r
?mean
?mad
example(mad)
example(boxplot)
```

Simply typing the name of the function, without parentheses, and hitting return will show the source code of the function.

The manual page contains a **description**, example **usage**, explanation of all **arguments**, further **details**, explanation of the returned **value**, **references**, **see also** linking to other functions, and **examples**.

## Package help

To read over all the help files for all functions of a given package, we can launch a clickable index in your default web browser.


```r
help(package="genefilter", help_type="html")
```

If you are using RStudio, you can launch the same index within RStudio's help pane with:


```r
help(package="genefilter")
```

Note that you need to have the package installed in order to access this documentation.

Another quick way to look up functions in a given package is to write out the package name, two ":" symbols and then trying tab-completion to get a list of functions which are exported by that package.


```r
library(geneplotter)
```


```r
geneplotter::
```

## Object help

If you have a question about a particular object in R, you might want to look up the help for the "class" of that object, which will tell you how to construct such an object and what methods are available for manipulating such objects. For example, we can find the name of the class of an object and look up help:


```r
class(6)
?numeric
?"numeric-class"
```

Sometimes, the constructor function and the class name will point to the same help page, although this is not necessarily true for all packages.


```r
library(Biobase)
```


```r
?ExpressionSet
?"ExpressionSet-class"
```

A quick way to find out what methods are available for a given class:


```r
methods(class="ExpressionSet")
methods(class="lm")
```

In `rafalib` we have also worked on a function (with help from Gabe Becker) which will provide a list of methods which are available for a given object:


```r
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
data(sample.ExpressionSet)
whatMethods(sample.ExpressionSet)
```

```
## [1] "coerce"          "esApply"         "exprs"           "exprs<-"        
## [5] "initialize"      "makeDataPackage" "rowMedians"      "rowQ"           
## [9] "write.exprs"
```

This is unfortunately only a partial list. If we look up `?ExpressionSet` we see that this class of object inherits from another class `eSet`. Knowing this, we can see more methods which are available:


```r
whatMethods("eSet")
```

```
##  [1] "["                "[["               "[[<-"            
##  [4] "$"                "$<-"              "abstract"        
##  [7] "annotation"       "annotation<-"     "assayData"       
## [10] "assayData<-"      "coerce"           "combine"         
## [13] "description"      "description<-"    "dim"             
## [16] "dimnames"         "dimnames<-"       "dims"            
## [19] "experimentData"   "experimentData<-" "fData"           
## [22] "fData<-"          "featureData"      "featureData<-"   
## [25] "featureNames"     "featureNames<-"   "fvarLabels"      
## [28] "fvarLabels<-"     "fvarMetadata"     "fvarMetadata<-"  
## [31] "initialize"       "notes"            "notes<-"         
## [34] "pData"            "pData<-"          "phenoData"       
## [37] "phenoData<-"      "preproc"          "preproc<-"       
## [40] "protocolData"     "protocolData<-"   "pubMedIds"       
## [43] "pubMedIds<-"      "sampleNames"      "sampleNames<-"   
## [46] "show"             "storageMode"      "storageMode<-"   
## [49] "updateObject"     "updateObjectTo"   "varLabels"       
## [52] "varLabels<-"      "varMetadata"      "varMetadata<-"
```

## Source code

You can find the source code for many functions by typing out the name of the function without () and pressing enter.


```r
mad
```

```
## standardGeneric for "mad" defined from package "stats"
## 
## function (x, center = median(x), constant = 1.4826, na.rm = FALSE, 
##     low = FALSE, high = FALSE) 
## standardGeneric("mad")
## <environment: 0x107373838>
## Methods may be defined for arguments: x, center, constant, na.rm, low, high
## Use  showMethods("mad")  for currently available ones.
```

You might have to specify a particular class if you want source code for a method:


```r
plotMA
```

```
## nonstandardGenericFunction for "plotMA" defined from package "BiocGenerics"
## 
## function (object, ...) 
## {
##     standardGeneric("plotMA")
## }
## <environment: 0x102ed30e0>
## Methods may be defined for arguments: object
## Use  showMethods("plotMA")  for currently available ones.
```

```r
showMethods("plotMA")
```

```
## Function: plotMA (package BiocGenerics)
## object="ANY"
## object="data.frame"
```

```r
getMethod("plotMA","data.frame")
```

```
## Method Definition:
## 
## function (object, ...) 
## {
##     .local <- function (object, ylim = NULL, colNonSig = "gray32", 
##         colSig = "red3", colLine = "#ff000080", log = "x", cex = 0.45, 
##         xlab = "mean expression", ylab = "log fold change", ...) 
##     {
##         if (!(ncol(object) == 3 & inherits(object[[1]], "numeric") & 
##             inherits(object[[2]], "numeric") & inherits(object[[3]], 
##             "logical"))) {
##             stop("When called with a data.frame, plotMA expects the data frame to have 3 columns, two numeric ones for mean and log fold change, and a logical one for significance.")
##         }
##         colnames(object) <- c("mean", "lfc", "sig")
##         object = subset(object, mean != 0)
##         py = object$lfc
##         if (is.null(ylim)) 
##             ylim = c(-1, 1) * quantile(abs(py[is.finite(py)]), 
##                 probs = 0.99) * 1.1
##         plot(object$mean, pmax(ylim[1], pmin(ylim[2], py)), log = log, 
##             pch = ifelse(py < ylim[1], 6, ifelse(py > ylim[2], 
##                 2, 16)), cex = cex, col = ifelse(object$sig, 
##                 colSig, colNonSig), xlab = xlab, ylab = ylab, 
##             ylim = ylim, ...)
##         abline(h = 0, lwd = 4, col = colLine)
##     }
##     .local(object, ...)
## }
## <environment: namespace:geneplotter>
## 
## Signatures:
##         object      
## target  "data.frame"
## defined "data.frame"
```

## Vignettes

"Vignettes" are documents which accompany R packages and are required for every Bioconductor package. They typically show an example workflow of the functions of the package using "chunks" of code with descriptive text, exactly as the document you are currently reading. 

You can find Bioconductor vignettes in PDF or R script form on the Bioconductor website, but they are even easier to access directly through R. Furthermore, accessing vignettes through R guarantees that the vignette is for the correct version of the package that you are using. The following code will list the names of the vignettes for a given package:


```r
vignette(package="Biobase")
```

A further call to `vignette` with the name of the vignette will launch a PDF viewer:


```r
vignette("ExpressionSetIntroduction")
```

In addition, an HTML browser can be launched with links to the various vignettes of a package:



```r
browseVignettes(package="Biobase")
```

