---
layout: page
title: Getting Started Exercises
---

## Exercises

Here we will test some of the basics of R data manipulation which you should know or should have learned by following the tutorials above. You will need to have the file `femaleMiceWeights.csv` in your working directory. As we showed above, one way to do this is using the `downloader` package:


```r
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv" 
download(url, destfile=filename)
```

1. Read in the file `femaleMiceWeights.csv` and report the body weight of the mouse in the exact name of the column containing the weights



2. The `[` and `]` symbols can be used to extract specific rows and specific columns of the table.  What is the entry in the 12th row and second column?



3. You should have learned how the `$` character is used to extract a column from a table and return it as a vector. Use `$` to extract the weight column and report the weight of the mouse in the 11th row.



4. The `length` function returns the number of elements in a vector. How many mice are included in our dataset?



5. To create a vector with the numbers 3 to 7, we can use `seq(3,7)` or, because they are consecutive, `3:7`. View the data and determine what rows are associated with the high fat or `hf` diet. Then use the `mean` function to compute the average weight of these mice.




6. One of the functions we will be using often is `sample`. Read the help file for sample using `?sample`. Now take a random sample of size 1 from the numbers 13 to 24 and report back the weight of the mouse represented by that row. Make sure to type `set.seed(1)` to ensure that everybody gets the same answer.


  
