---
title: Population, Samples, and Estimates Exercises
layout: page
---

A>## Exercises
A>
A>For these exercises, we will be using the following dataset:
A>
A>
```r
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 
```
A>
A>We will remove the lines that contain missing values:
A>
```r
dat <- na.omit( dat )
```
A>
A>1. Use `dplyr` to create a vector `x` with the bodyweight of all males on the control (`chow`) diet. What is this population's average?
A>
A>
A>2. Now use the `rafalib` package and use the `popsd` function to compute the population standard deviation.
A>
A>
A>3. Set the seed at 1. Take a random sample {$$}X{/$$} of size 25 from `x`. What is the sample average?
A>
A>
A>
A>4. Use `dplyr` to create a vector `y` with the bodyweight of all males on the high fat (`hf`) diet. What is this population's average?
A>
A>
A>5. Now use the `rafalib` package and use the `popsd` function to compute the population standard deviation.
A>
A>
A>6. Set the seed at 1. Take a random sample {$$}Y{/$$} of size 25 from `y`. What is the sample average?
A>
A>
A>7. What is the difference in absolute value between {$$}\bar{y} - \bar{x}{/$$} and {$$}\bar{X}-\bar{Y}{/$$}?
A>
A>
A>8. Repeat the above for females. Make sure to set the seed to 1 before each `sample` call. What is the difference in absolute value between {$$}\bar{y} - \bar{x}{/$$} and {$$}\bar{X}-\bar{Y}{/$$}?
A>
A>
A>9. For the females, our sample estimates we were closer to the population difference than with males. What is a possible explanation for this?
A>  - A) The population variance of the females is smaller than that of the males; thus the sample variable has less variability.
A>  - B) Statistical estimates are more precise for females.
A>  - C) The sample size was larger for females.
A>  - D) The sample size was smaller for females.
A>  
A>
A>
