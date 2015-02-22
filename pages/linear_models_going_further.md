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

In the standard linear model, 

glm

[Quick R](http://www.statmethods.net/advstats/glm.html)

## Mixed effects linear models

Linear models with fixed and random effects

lme4

## Bayesian linear models

stan, BUGS

## Penalized linear models

Lasso, ridge, elastic neet

## Many simultaneous linear models

(Y is not M x 1, but M x N)

limma
