---
title: "Going further with linear models"
author: "Mike"
date: "02/19/2015"
output: html_document
layout: page
---

Linear models can be extended in many directions. Here are some examples of extensions, which you might come across in analyzing data in the life sciences:

## Robust linear models

In calculating the solution and its estimated error in the standard linear model, we minimize the squared errors. This involves a sum of squares from all the data points, which means that a few *outlier* data points can have a large influence on the solution. In addition, the errors are assumed to be have constant variance (called *homoskedasticity*), which might not always hold true (when this is not true, it is called *heteroskedasticity*). Methods have been developed therefore to generate more *robust* solutions, which behave well in the presence of outliers, or when the distributional assumptions are not met. A number of these are mentioned on the [robust statistics](http://cran.r-project.org/web/views/Robust.html) page on the CRAN website. For more background, there is also a [Wikipedia article](http://en.wikipedia.org/wiki/Robust_regression) with references.

## Generalized linear models

In the standard linear model, we did not make any assumptions about the distribution of $\mathbf{Y}$, though in some cases we can gain better estimates if we know that $\mathbf{Y}$ is, for example restricted to non-negative integers $0,1,2,\dots$, or restricted to the interval $[0,1]$. A framework for analyzing such cases is referred to as *generalized linear models*, commonly abbreviated as GLMs. The two key components of the GLM are the *link function* and a probability distribution. The link function $g$ connects our familiar matrix product $\mathbf{X} \boldsymbol{\beta}$ to the $\mathbf{Y}$ values through:

$$ \textrm{E}(\mathbf{Y}) = g^{-1}( \mathbf{X} \boldsymbol{\beta} ) $$

There is a function in base R for fitting GLMs, which is `glm` and uses a familiar form as `lm`, with additional arguments including `family` which specifies the distributional assumption of $\mathbf{Y}$. Some examples of the use of GLMs are shown at the [Quick R](http://www.statmethods.net/advstats/glm.html) website. There are a number of references for GLMs on the [Wikipedia page](http://en.wikipedia.org/wiki/Generalized_linear_model). 

## Mixed effects linear models

In the standard linear model, we assumed that the matrix $\mathbf{X}$ was *fixed* and not random. For example, we measured the frictional coefficients for each leg pair, and in the push and pull direction. The fact that an observation had a $1$ for a given column in $\mathbf{X}$ was not random, but dictated by the experimental design. However, in the father and son heights example, we did not fix the values of the father's heights, but observed these (and likely these were measured with some error). A framework for studying the effect of the randomness for various columns in $X$ is referred to as *mixed effects* models, which implies that some effects are *fixed* and some effects are *random*. One of the most popular packages in R for fitting linear mixed effects models is [lme4](http://lme4.r-forge.r-project.org/) which has an accompanying paper on [Fitting Linear Mixed-Effects Models using lme4](http://arxiv.org/abs/1406.5823). There is also a [Wikipedia page](http://en.wikipedia.org/wiki/Mixed_model) with more references.

## Bayesian linear models

The approach presented here focused the estimation of $\boldsymbol{beta}$ which minimized the squared error and then calculating its standard error. An alternative approach to the statistical inference about $\boldsymbol{\beta}$ is using *Bayesian* methods. Bayesian methods first involve the formulation of *prior* distributional beliefs about model parameters, and then calculating either directly or computing using computational methods the *posterior* distribution of $\boldsymbol{\beta}$, from which we can examine the *posterior mode* (the most likely value), and [credible intervals](http://en.wikipedia.org/wiki/Credible_interval). In addition, many models can be connected together in what is referred to as a *hierarchical model*. A good reference for Bayesian hierarchical models is [Bayesian Data Analysis](http://www.stat.columbia.edu/~gelman/book/), and some software for computing Bayesian linear models can be found on the [Bayes](http://cran.r-project.org/web/views/Bayesian.html) page on CRAN. Some well known software for computing Bayesian models are [stan](http://mc-stan.org/) and [BUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/).

## Penalized linear models

Penalized linear models introduce an penalty term to the minimization in addition to the squares of the residuals. These are typically of the form, $\lambda \sum_i \|\beta_i\|^k$, for different numbers $k$. The motivation for these is to introduce some stability of the $\boldsymbol{\beta}$.

https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/lm.ridge.html

http://statweb.stanford.edu/~tibs/ElemStatLearn/

http://cran.r-project.org/web/packages/lars/index.html

http://cran.r-project.org/web/packages/glmnet/index.html

## Many simultaneous linear models

$\mathbf{Y}$ is not M x 1, but M x N

limma

Statistica Sinica 12(2002), 31-46. REPLICATED MICROARRAY DATA. Ingrid Lönnstedt and Terry Speed

http://www.statsci.org/smyth/pubs/

http://www.ncbi.nlm.nih.gov/pubmed/16646809

Stat Appl Genet Mol Biol. 2004;3:Article3. Epub 2004 Feb 12.
Linear models and empirical bayes methods for assessing differential expression in microarray experiments.
Smyth GK1.

http://www.bioconductor.org/packages/release/bioc/html/limma.html

UserGuide

