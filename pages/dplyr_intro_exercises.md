---
layout: page
title:  dplyr exercises
---

## Exercises

For these exercises we will use a new dataset related to mammalian sleep. This data is described 
[here](http://docs.ggplot2.org/0.9.3.1/msleep.html). Download the CSV file from this location:



We are going to read in this data, then test your knowledge of key `dplyr` functions `select` and `filter`. We are also going to review two different _classes_: data frames and vectors.

1. Read in the `msleep_ggplot2.csv` file with the function `read.csv` and use the function `class` to determine what type of object is returned.



2. Now use the `filter` function to select only the primates. How many animals in the table are primates? Hint: the `nrow` function gives you the number of rows of a data frame or matrix.



3. What is the class of the object you obtain after subsetting the table to only include primates?



4. Now use the `select` function to extract the sleep (total) for the primates. What class is this object? Hint: use `%>%` to pipe the results of the  `filter` function to `select`.



5. Now we want to calculate the average of sleep for primates (the average of the numbers computed above). One challenge is that the `mean` function requires a vector so, if we simply apply it to the output above, we get an error. Look at the help file for `unlist` and use it to compute the desired average.



6. For the last exercise, we could also use the dplyr `summarize` function. We have not introduced this function, but you can read the help file and repeat exercise 5 this time using just `filter` and `summarize` to get the answer.




