---
layout: page
title: Inference for high throughput exercises
---

## Exercises

These exercises will help clarify that p-values are random variables and some of the properties of these p-values. Note that just like the sample average is a random variable because it is based on a random sample, the p-values are based on random variables (sample mean and sample standard deviation for example) and thus it is also a random variable.

To see this, let's see how p-values change when we take different samples.

```r
set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)
```

1. What proportion of the p-values is below 0.05?



2. What proportion of the p-values is below 0.01?



3. Assume you are testing the effectiveness of 20 diets on mice weight. For each of the 20 diets, you run an experiment with 10 control mice and 10 treated mice. Assume the null hypothesis, that the diet has no effect, is true for all 20 diets and that mice weights follow a normal distribution, with mean 30 grams and a standard deviation of 2 grams. Run a Monte Carlo simulation for one of these studies:

    
    ```r
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)
    ```

    Question: Now run a Monte Carlo simulation imitating the results for the experiment for all 20 diets. If you set the seed at 100, `set.seed(100)`, how many of p-values are below 0.05?




4. Now create a simulation to learn about the distribution of the number of p-values that are less than 0.05. In question 3, we ran the 20 diet experiment once. Now we will run the experiment 1,000 times and each time save the number of p-values that are less than 0.05. 

    Set the seed at 100, `set.seed(100)`, run the code from Question 3 1,000 times, and save the number of times the p-value is less than 0.05 for each of the 1,000 instances. What is the average of these numbers? This is the expected number of tests (out of the 20 we run) that we will reject when the null is true.



5. What this says is that on average, we expect some p-value to be 0.05 even when the null is true for all diets. Use the same simulation data and report for what percent of the 1,000 replications did we make more than 0 false positives?


