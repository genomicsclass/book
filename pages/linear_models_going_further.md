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

stan, BUGS

## Penalized linear models

Lasso, ridge, elastic neet

## Many simultaneous linear models

(Y is not M x 1, but M x N)

limma
