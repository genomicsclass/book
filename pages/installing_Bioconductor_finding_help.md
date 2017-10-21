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

There are many ways to find help directly from within R. 
A very comprehensive help page is generated with the command `help.start()` from the
command line.  This gives links to manuals and manual pages for all documentation in the environment.
This documentation set is dependent on the set of packages that are installed.  You may
have seen a help page for a function in a given R session, that cannot be found in
another R session of the same version.  Only if the package to which the function belongs
has been installed will the help page be found.

The [Bioconductor support site](https://support.bioconductor.org/) is modeled on [biostar](https://github.com/ialbert/biostar-central) and is a very active and congenial forum for detailed questions about
Bioconductor tools and analyses.  The [guide for posting](http://www.bioconductor.org/help/support/posting-guide/)
is worth a visit.

## Help for functions

Typically, every function will have its own manual page which is accessible by typing a question mark ?, followed by the function name and hitting return.


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

R has good capabilities for self-description.  Classes can be formally linked to methods
that operate usefully on their instances.  The methods available can be
listed using the `methods` function.


```r
data(sample.ExpressionSet)
methods(class=class(sample.ExpressionSet))
```

```
##  [1] [                [[               [[<-             $$               
##  [5] $$<-              abstract         annotation       annotation<-    
##  [9] as.data.frame    assayData        assayData<-      classVersion    
## [13] classVersion<-   coerce           combine          description     
## [17] description<-    dim              dimnames         dimnames<-      
## [21] dims             esApply          experimentData   experimentData<-
## [25] exprs            exprs<-          fData            fData<-         
## [29] featureData      featureData<-    featureNames     featureNames<-  
## [33] fvarLabels       fvarLabels<-     fvarMetadata     fvarMetadata<-  
## [37] initialize       isCurrent        isVersioned      KEGG2heatmap    
## [41] KEGGmnplot       makeDataPackage  Makesense        notes           
## [45] notes<-          pData            pData<-          phenoData       
## [49] phenoData<-      preproc          preproc<-        protocolData    
## [53] protocolData<-   pubMedIds        pubMedIds<-      rowMedians      
## [57] rowQ             sampleNames      sampleNames<-    show            
## [61] storageMode      storageMode<-    updateObject     updateObjectTo  
## [65] varLabels        varLabels<-      varMetadata      varMetadata<-   
## [69] write.exprs     
## see '?methods' for accessing help and source code
```

## Source code

You can find the source code for many functions by typing out the name of the function without () and pressing enter.


```r
read.csv
```

```
## function (file, header = TRUE, sep = ",", quote = "\\"", dec = ".", 
##     fill = TRUE, comment.char = "", ...) 
## read.table(file = file, header = header, sep = sep, quote = quote, 
##     dec = dec, fill = fill, comment.char = comment.char, ...)
## <bytecode: 0x7f8e35b19200>
## <environment: namespace:utils>
```

Note that this just 'wraps' a call to `read.table`, which is much more involved, but can
be printed to the terminal by mentioning it to the interpreter.

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
## <environment: 0x7f8e338381b8>
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

## Analysis help

The [workflows repository](http://www.bioconductor.org/help/workflows/) 
contains extensive _computable documents_ describing how
certain common analyses can be performed, typically involving published
data from realistic studies.

## Summary

- All R functions from Bioconductor packages will have _manual pages_ and
most will have _running examples_;
- All Bioconductor software packages will have _vignettes_ that describe
how functions work together to achieve the objectives of the package;
- Many _workflow documents_ are available that describe how functions from
different packages can be used together to complete analyses;
- R has significant _built-in utilities_ for finding documents about basic infrastructure
and statistical procedures;
- R's [mailing list for elementary questions](https://stat.ethz.ch/mailman/listinfo/r-help)
and Bioconductor's [support site](https://support.bioconductor.org/) are useful fora
with extensive archives.
