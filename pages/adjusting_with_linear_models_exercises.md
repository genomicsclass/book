---
Title: Adjusting with linear models exercises
---

A>## Exercises
A>
A>For the dataset we have been working with models do not help due to the almost perfect confounding. This is one reason we created the subset dataset:
A>
A>
```r
library(GSE5859Subset)
data(GSE5859Subset)
```
A>
A>Here we purposely confounded month and group (sex) but not completely:
A>
A>
```r
sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)
```
A>
A>
A>1. Using the functions `rowttests` and `qvalue` compare the two groups. Because this is a smaller dataset which decreases our power, we will use more lenient FDR cut-off of 10%.
A>
A>    How many gene have q-values less than 0.1? 
A>
A>
A>
A>2. Note that `sampleInfo$group` here presents males and females. Thus we expect differences to be in on chrY and, for genes that escape inactivation, chrX. Note that we do not expect many autosomal genes to be difference between males and females. This gives us an opportunity to evaluate false and true positives with experimental data. For example, we evaluate results using the proportion genes of the list that are on chrX or chrY.
A>
A>    For the list calculated above, what proportion of this list is on chrX or chrY?
A>
A>
A>
A>3. We can also check how many of the chrX and chrY we detected as different. How many are on Y?
A>
A>
A>4. Now for the autosomal genes (not on chrX and chrY) for which q-value < 0.1 perform a t-test comparing
A>samples processed in June to those processed in October. 
A>
A>    What proportion of these have p-values <0.05 ?
A>
A>
A>
A>5. The above result shows that the great majority of the autosomal genes show differences due to processing data. This provides further evidence that confounding is resulting in false positives. So we are going to try to model the month effect to better estimate the sex effect. We are going to use a linear model:
A>
A>    Which of the following creates the appropriate design matrix?
A>    - A) `X = model.matrix(~sex+ethnicity)`
A>    - B) `X = cbind(sex,as.numeric(month))`  
A>    - C) It can't be done with one line
A>    - D) `X = model.matrix(~sex+month)`
A>   
A>
A>
A>6. Now use the `X` defined above to fit a regression model using `lm` for each gene. Note that you can obtain p-values for estimated parameters using `summary`. Here is an example
A>
A>    
    ```r
    X = model.matrix(~sex+month)
    i = 234
    y = geneExpression[i,]
    fit = lm(y~X)
    summary(fit)$coef
    ```
A>
A>
A>    How many of the q-values for the group comparison are <0.1 now?
A>
A>
A>    Note the big drop from what we obtained without the correction. 
A>
A>7. With this new list, what proportion of these are chrX and chrY?
A>
A>
A>    Note the big improvement.
A>
A>8. How many on Y or X?
A>
A>
A>9. Now, from the linear model above, extract the p-values related to the coefficient representing the October versus June differences using the same linear model.
A>
A>    How many of the q-values for the month comparison are <0.1 now?
A>
A>
A>
A>    Note that this approach is basically the approach implemented by Combat.
A>
