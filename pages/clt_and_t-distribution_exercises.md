---
title: "Central Limit Theorem and t-distribution Exercise"
layout: page
---


## Exercises

For these exercises, we will be using the following dataset:


```r
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- na.omit( read.csv(filename) )
```


1. If a list of numbers has a distribution that is well approximated by the normal distribution, what proportion of these numbers are within one standard deviation away from the list's average?


2. What proportion of these numbers are within two standard deviations away from the list's average?


3. What proportion of these numbers are within three standard deviations away from the list's average?


4. Define `y` to be the weights of males on the control diet. What proportion of the mice are within one standard deviation away from the average weight (remember to use `popsd` for the population `sd`)? 


5. What proportion of these numbers are within two standard deviations away from the list's average?


6. What proportion of these numbers are within three standard deviations away from the list's average?


7. Note that the numbers for the normal distribution and our weights are relatively close. Also, notice that we are indirectly comparing quantiles of the normal distribution to quantiles of the mouse weight distribution. We can actually compare all quantiles using a qqplot. Which of the following best describes the qq-plot comparing mouse weights to the normal distribution?
    - A) The points on the qq-plot fall exactly on the identity line.
    - B) The average of the mouse weights is not 0 and thus it can't follow a normal distribution.
    - C) The mouse weights are well approximated by the normal distribution, although the larger values (right tail) are larger than predicted by the normal. This is consistent with the differences seen between question 3 and 6. 
    - D) These are not random variables and thus they can't follow a normal distribution.
  


8. Create the above qq-plot for the four populations: male/females on each of the two diets. What is the most likely explanation for the mouse weights being well approximated? What is the best explanation for all these being well approximated by the normal distribution?
    - A) The CLT tells us that sample averages are approximately normal.
    - B) This just happens to be how nature behaves. Perhaps the result of many biological factors averaging out.
    - C) Everything measured in nature follows a normal distribution.
    - D) Measurement error is normally distributed.
  


9. Here we are going to use the function `replicate` to learn about the distribution of random variables. All the above exercises relate to the normal distribution as an approximation of the distribution of a fixed list of numbers or a population. We have not yet discussed probability in these exercises. If the distribution of a list of numbers is approximately normal, then if we pick a number at random from this distribution, it will follow a normal distribution. However, it is important to remember that stating that some quantity has a distribution does not necessarily imply this quantity is random. Also, keep in mind that this is not related to the central limit theorem. The central limit applies to averages of random variables. Let's explore this concept. 

    We will now take a sample of size 25 from the population of males on the chow diet. The average of this sample is our random variable. We will use the `replicate` to observe 10,000 realizations of this random variable. Set the seed at 1, generate these 10,000 averages. Make a histogram and qq-plot of these 10,000 numbers against the normal distribution. 
    
    We can see that, as predicted by the CLT, the distribution of the random variable is very well approximated by the normal distribution.

    
    ```r
    y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
    avgs <- replicate(10000, mean( sample(y, 25)))
    mypar(1,2)
    hist(avgs)
    qqnorm(avgs)qqline(avgs)
    ```

    What is the average of the distribution of the sample average?


10. What is the standard deviation of the distribution of sample averages?


11. According to the CLT, the answer to exercise 9 should be the same as `mean(y)`. You should be able to confirm that these two numbers are very close. Which of the following does the CLT tell us should be close to your answer to exercise 10?
    - A) `popsd(y)`
    - B) `popsd(avgs)/sqrt(25)`
    - C) `sqrt(25) / popsd(y)`
    - D) `popsd(y)/sqrt(25)`
  


12. In practice we do not know $$\sigma$$ (`popsd(y)`) which is why we can't use the CLT directly. This is because we see a sample and not the entire distribution. We also can't use `popsd(avgs)` because to construct averages, we have to take 10,000 samples and this is never practical. We usually just get one sample. Instead we have to estimate `popsd(y)`. As described, what we use is the sample standard deviation. Set the seed at 1, using the `replicate` function, create 10,000 samples of 25 and now, instead of the sample average, keep the standard deviation. Look at the distribution of the sample standard deviations. It is a random variable. The real population SD is about 4.5. What proportion of the sample SDs are below 3.5?


13. What the answer to question 12 reveals is that the denominator of the t-test is a random variable. By decreasing the sample size, you can see how this variability can increase. It therefore adds variability. The smaller the sample size, the more variability is added. The normal distribution stops providing a useful approximation. When the distribution of the population values is approximately normal, as it is for the weights, the t-distribution provides a better approximation. We will see this later on. Here we will look at the difference between the t-distribution and normal. Use the function `qt` and `qnorm` to get the quantiles of `x=seq(0.0001,0.9999,len=300)`. Do this for degrees of freedom 3, 10, 30, and 100. Which of the following is true?
    - A) The t-distribution and normal distribution are always the same.
    - B) The t-distribution has a higher average than the normal distribution.
    - C) The t-distribution has larger tails up until 30 degrees of freedom, at which point it is practically the same as the normal distribution.
    - D) The variance of the t-distribution grows as the degrees of freedom grow.






