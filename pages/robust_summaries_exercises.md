---
layout: page
title: Robust summaries exercises
---

## Exercises

We are going to explore the properties of robust statistics. We will use one of the datasets included in R, which contains weight of chicks in grams as they grow from day 0 to day 21. This dataset also splits up the chicks by different protein diets, which are coded from 1 to 4. We use this dataset to also show an important operation in R (not related to robust summaries): `reshape`.

This dataset is built into R and can be loaded with:


```r
data(ChickWeight)
```

To begin, take a look at the weights of all observations over time and color the points to represent the Diet:


```r
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
```


First, notice that the rows here represent time points rather than individuals. To facilitate the comparison of weights at different time points and across the different chicks, we will reshape the data so that each row is a chick. In R we can do this with the `reshape` function:


```r
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")
```

The meaning of this line is: reshape the data from _long_ to _wide_, where the columns Chick and Diet are the ID's and the column Time indicates different observations for each ID. Now examine the head of this dataset:


```r
head(chick)
```
We also want to remove any chicks that have missing observations at any time points (NA for "not available"). The following line of code identifies these rows and then removes them:


```r
chick = na.omit(chick)
```
1. Focus on the chick weights on day 4 (check the column names of 'chick' and note the numbers). How much does the average of chick weights at day 4 increase if we add an outlier measurement of 3000 grams? Specifically, what is the average weight of the day 4 chicks, including the outlier chick, divided by the average of the weight of the day 4 chicks without the outlier. Hint: use `c`to add a number to a vector.





2. In exercise 1, we saw how sensitive the mean is to outliers. Now let's see what happens when we use the median instead of the mean.  Compute the same ratio, but now using median instead of mean. Specifically, what is the median weight of the day 4 chicks, including the outlier chick, divided by the median of the weight of the day 4 chicks without the outlier.



3. Now try the same thing with the sample standard deviation (the `sd` function in R). Add a chick with weight 3000 grams to the chick weights from day 4. How much does the standard deviation change? What's the standard deviation with the outlier chick divided by the standard deviation without the outlier chick?



4. Compare the result above to the median absolute deviation in R, which is calculated with the `mad` function. Note that the mad is unaffected by the addition of a single outlier. The `mad` function in R includes the scaling factor 1.4826, such that `mad` and `sd` are very similar for a sample from a normal distribution. What's the MAD with the outlier chick divided by the MAD without the outlier chick?



5. Our last question relates to how the Pearson correlation is affected by an outlier as compared to the Spearman correlation. The Pearson correlation between x and y is given in R by `cor(x,y)`. The Spearman correlation is given by `cor(x,y,method="spearman")`. 

    Plot the weights of chicks from day 4 and day 21. We can see that there is some general trend, with the lower weight chicks on day 4 having low weight again on day 21, and likewise for the high weight chicks.
    
    Calculate the Pearson correlation of the weights of chicks from day 4 and day 21. Now calculate how much the Pearson correlation changes if we add a chick that weighs 3000 on day4 and 3000 on day 21. Again, divide the Pearson correlation with the outlier chick over the Pearson correlation computed without the outliers.




6. Save the weights of the chicks on day 4 from diet 1 as a vector `x`. Save the weights of the chicks on day 4 from diet 4 as a vector `y`. Perform a t-test comparing `x` and `y` (in R the function `t.test(x,y)` will perform the test). Then perform a Wilcoxon test of `x` and `y` (in R the function `wilcox.test(x,y)` will perform the test). A warning will appear that an exact p-value cannot be calculated with ties, so an approximation is used, which is fine for our purposes.

    Perform a t-test of x and y, after adding a single chick of weight 200 grams to `x` (the diet 1 chicks). What is the p-value from this test? The p-value of a test is available with the following code: `t.test(x,y)$p.value`



7. Do the same for the Wilcoxon test. The Wilcoxon test is robust to the outlier.  In addition, it has less assumptions that the t-test on the distribution of the underlying data.



8. We will now investigate a possible downside to the Wilcoxon-Mann-Whitney test statistic. Using the following code to make three boxplots, showing the true Diet 1 vs 4 weights, and then two altered versions: one with an additional difference of 10 grams and one with an additional difference of 100 grams. Use the x and y as defined above, NOT the ones with the added outlier.

    
    ```r
    library(rafalib)
    mypar(1,3)
    boxplot(x,y)
    boxplot(x,y+10)
    boxplot(x,y+100)
    ```
  
    What is the difference in t-test statistic (obtained by `t.test(x,y)$statistic`) between adding 10 and adding 100 to all the values in the group 'y'? Take the the t-test statistic with x and y+10 and subtract the t-test statistic with x and y+100. The value should be positive.



9. Examine the Wilcoxon test statistic for x and y+10 and for x and y+100. Because the Wilcoxon works on ranks, once the two groups show complete separation, that is all points from group 'y' are above all points from group 'x', the statistic will not change, regardless of how large the difference grows. Likewise, the p-value has a minimum value, regardless of how far apart the groups are. This means that the Wilcoxon test can be considered less powerful than the t-test in certain contexts. In fact, for small sample sizes, the p-value can't be very small, even when the difference is very large. What is the p-value if we compare `c(1,2,3)` to `c(4,5,6)` using a Wilcoxon test?


10. What is the p-value if we compare `c(1,2,3)` to `c(400,500,600)` using a Wilcoxon test?


