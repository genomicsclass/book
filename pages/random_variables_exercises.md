---
layout: page
title:  Random Variables Exercises
---

A>## Exercises
A>
A>For these exercises, we will be using the following dataset:
A>
A>
```r
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )
```
A>
A>Here `x` represents the weights for the entire population.
A>
A>1. What is the average of these weights?
A>
A>
A>
A>2. After setting the seed at 1, `set.seed(1)` take a random sample of size 5. What is the absolute value (use `abs`) of the difference between the average of the sample and the average of all the values?
A>
A>
A>
A>
A>3. After setting the seed at 5, `set.seed(5)` take a random sample of size 5. What is the absolute value of the difference between the average of the sample and the average of all the values?
A>
A>
A>
A>
A>4. Why are the answers from 2 and 3 different?
A>  - A) Because we made a coding mistake.
A>  - B) Because the average of the `x` is random.
A>  - C) Because the average of the samples is a random variable.
A>  - D) All of the above.
A>
A>
A>
A>  
A>5. Set the seed at 1, then using a for-loop take a random sample of 5 mice 1,000 times. Save these averages. What percent of these 1,000 averages are more than 1 ounce away from the average of `x` ?
A>
A>
A>
A>6. We are now going to increase the number of times we redo the sample from 1,000 to 10,000. Set the seed at 1, then using a for-loop take a random sample of 5 mice 10,000 times. Save these averages. What percent of these 10,000 averages are more than 1 ounce away from the average of `x` ?
A>
A>
A>
A>7. Note that the answers of 4 and 5 barely changed. This is expected. The way we think about the random value distributions is as the distribution of the list of values obtained if we repeated the experiment an infinite number of times. On a computer, we can't perform an infinite number of iterations so instead, for our examples, we consider 1,000 to be large enough, thus 10,000 is as well. Now, if instead we change the sample size, then we change the random variable and thus its distribution. 
A>
A>    Set the seed at 1, then using a for-loop take a random sample of 50 mice 1,000 times. Save these averages. What percent of these 1,000 averages are more than 1 ounce away from the average of `x` ?
A>
A>
A>
A>8. Use a histogram to "look" at the distribution of averages we get with a sample size of 5 and a sample size of 50. How would you say they differ?
A>  - A) They are actually the same.
A>  - B) They both look roughly normal, but with a sample size of 50 the spread is smaller.
A>  - C) They both look roughly normal, but with a sample size of 50 the spread is larger.
A>  - D) The second distribution does not look normal at all.
A>
A>
A>
A>9. For the last set of averages, the ones obtained from a sample size of 50, what percent are between 23 and 25?
A>
A>
A>
A>10. Now ask the same question of a normal distribution with average 23.9 and standard deviation 0.43
A>
A>
A>The answer to 9 and 10 were very similar. This is because we can approximate the distribution of the sample average with a normal distribution. We will learn more about why this is next. 
A>
A>
A>
