---
title: Population, Samples, and Estimates Exercises
layout: page
---

## Exercises

For these exercises, we will be using the following dataset:


```r
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 
```

We will remove the lines that contain missing values:

```r
dat <- na.omit( dat )
```

1. Use `dplyr` to create a vector `x` with the bodyweight of all males on the control (`chow`) diet. What is this population's average?


2. Now use the `rafalib` package and use the `popsd` function to compute the population standard deviation.


3. Set the seed at 1. Take a random sample $$X$$ of size 25 from `x`. What is the sample average?



4. Use `dplyr` to create a vector `y` with the bodyweight of all males on the high fat (`hf`) diet. What is this population's average?


5. Now use the `rafalib` package and use the `popsd` function to compute the population standard deviation.


6. Set the seed at 1. Take a random sample $$Y$$ of size 25 from `y`. What is the sample average?


7. What is the difference in absolute value between $$\bar{y} - \bar{x}$$ and $$\bar{X}-\bar{Y}$?


8. Repeat the above for females. Make sure to set the seed to 1 before each `sample` call. What is the difference in absolute value between $$\bar{y} - \bar{x}$$ and $$\bar{X}-\bar{Y}$?


9. For the females, our sample estimates we were closer to the population difference than with males. What is a possible explanation for this?
  - A) The population variance of the females is smaller than that of the males; thus the sample variable has less variability.
  - B) Statistical estimates are more precise for females.
  - C) The sample size was larger for females.
  - D) The sample size was smaller for females.
  


