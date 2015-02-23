---
layout: page
title: dplyr tutorial
---



### What is dplyr?

dplyr is a powerful R-package to transform and summarize tabular data with rows and columns. For another explanation of dplyr see the dplyr package vignette: [Introduction to dplyr](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html)

### Why is it useful?

The package contains a set of functions (or "verbs") that perform common data manipulation operations such as filtering for rows, selecting specific columns, re-ordering rows, adding new columns and summarizing data. 

In addition, dplyr contains a useful function to perform another common task which is the "split-apply-combine" concept.  We will discuss that in a little bit. 

### How does it compare to using base functions R?

If you are familiar with R, you are probably familiar with base R functions such as split(), subset(), apply(), sapply(), lapply(), tapply() and aggregate(). Compared to base functions in R, the functions in dplyr are easier to work with, are more consistent in the syntax and are targeted for data analysis around data frames instead of just vectors. 

### How do I get dplyr? 

To install dplyr


```r
install.packages("dplyr")
```

To load dplyr


```r
library(dplyr)
```

```
## Error in library(dplyr): there is no package called 'dplyr'
```

# Data: mammals sleep

The msleep (mammals sleep) data set contains the sleeptimes and weights for a set of mammals and is available in the dagdata repository on github. This data set contains 83 rows and 11 variables.  

Download the msleep data set in CSV format from [here](https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv), and then load into R:


```r
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- "msleep_ggplot2.csv"
if (!file.exists(filename)) download(url,filename)
msleep <- read.csv("msleep_ggplot2.csv")
head(msleep)
```

```
##                         name      genus  vore        order conservation
## 1                    Cheetah   Acinonyx carni    Carnivora           lc
## 2                 Owl monkey      Aotus  omni     Primates         <NA>
## 3            Mountain beaver Aplodontia herbi     Rodentia           nt
## 4 Greater short-tailed shrew    Blarina  omni Soricomorpha           lc
## 5                        Cow        Bos herbi Artiodactyla domesticated
## 6           Three-toed sloth   Bradypus herbi       Pilosa         <NA>
##   sleep_total sleep_rem sleep_cycle awake brainwt  bodywt
## 1        12.1        NA          NA  11.9      NA  50.000
## 2        17.0       1.8          NA   7.0 0.01550   0.480
## 3        14.4       2.4          NA   9.6      NA   1.350
## 4        14.9       2.3   0.1333333   9.1 0.00029   0.019
## 5         4.0       0.7   0.6666667  20.0 0.42300 600.000
## 6        14.4       2.2   0.7666667   9.6      NA   3.850
```

The columns (in order) correspond to the following: 

column name | Description
--- | ---
name | common name
genus | taxonomic rank
vore | carnivore, omnivore or herbivore?
order | taxonomic rank
conservation | the conservation status of the mammal
sleep\_total | total amount of sleep, in hours
sleep\_rem | rem sleep, in hours
sleep\_cycle | length of sleep cycle, in hours
awake | amount of time spent awake, in hours
brainwt | brain weight in kilograms
bodywt | body weight in kilograms


# Important dplyr verbs to remember

dplyr verbs | Description
--- | ---
`select()` | select columns 
`filter()` | filter rows
`arrange()` | re-order or arrange rows
`mutate()` | create new columns
`summarise()` | summarise values
`group_by()` | allows for group operations in the "split-apply-combine" concept


# dplyr verbs in action

The two most basic functions are `select()` and `filter()` which selects columns and filters rows, respectively. 

### Selecting columns using `select()`

Select a set of columns: the name and the sleep\_total columns. 


```r
sleepData <- select(msleep, name, sleep_total)
```

```
## Error in eval(expr, envir, enclos): could not find function "select"
```

```r
head(sleepData)
```

```
## Error in head(sleepData): object 'sleepData' not found
```

To select all the columns *except* a specific column, use the "-" (subtraction) operator (also known as negative indexing)


```r
head(select(msleep, -name))
```

```
## Error in head(select(msleep, -name)): could not find function "select"
```

To select a range of columns by name, use the ":" (colon) operator


```r
head(select(msleep, name:order))
```

```
## Error in head(select(msleep, name:order)): could not find function "select"
```

To select all columns that start with the character string "sl", use the function `starts_with()`


```r
head(select(msleep, starts_with("sl")))
```

```
## Error in head(select(msleep, starts_with("sl"))): could not find function "select"
```

Some additional options to select columns based on a specific criteria include

1. `ends_with()` = Select columns that end with a character string
2. `contains()` = Select columns that contain a character string
3. `matches()` = Select columns that match a regular expression
4. `one_of()` = Select columns names that are from a group of names



### Selecting rows using `filter()`

Filter the rows for mammals that sleep a total of more than 16 hours. 


```r
filter(msleep, sleep_total >= 16)
```

```
## Error in filter(msleep, sleep_total >= 16): object 'sleep_total' not found
```

Filter the rows for mammals that sleep a total of more than 16 hours *and* have a body weight of greater than 1 kilogram.


```r
filter(msleep, sleep_total >= 16, bodywt >= 1)
```

```
## Error in match.arg(method): object 'bodywt' not found
```

Filter the rows for mammals in the Perissodactyla and Primates taxonomic order


```r
filter(msleep, order %in% c("Perissodactyla", "Primates"))
```

```
## Error in match(x, table, nomatch = 0L): 'match' requires vector arguments
```

You can use the boolean operators (e.g. >, <, >=, <=, !=, %in%) to create the logical tests. 

# Pipe operator: %>%

Before we go any futher, let's introduce the pipe operator: %>%. dplyr imports this operator from another package (magrittr).  This operator allows you to pipe the output from one function to the input of another function. Instead of nesting functions (reading from the inside to the outside), the idea of of piping is to read the functions from left to right. 

Here's an example you have seen:


```r
head(select(msleep, name, sleep_total))
```

```
## Error in head(select(msleep, name, sleep_total)): could not find function "select"
```

Now in this case, we will pipe the msleep data frame to the function that will select two columns (name and sleep\_total) and then pipe the new data frame to the function `head()` which will return the head of the new data frame. 


```r
msleep %>% 
    select(name, sleep_total) %>% 
    head
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

You will soon see how useful the pipe operator is when we start to combine many functions.  

# Back to dplyr verbs in action

Now that you know about the pipe operator (%>%), we will use it throughout the rest of this tutorial. 


### Arrange or re-order rows using `arrange()`

To arrange (or re-order) rows by a particular column such as the taxonomic order, list the name of the column you want to arrange the rows by


```r
msleep %>% arrange(order) %>% head
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

Now, we will select three columns from msleep, arrange the rows by the taxonomic order and then arrange the rows by sleep\_total. Finally show the head of the final data frame


```r
msleep %>% 
    select(name, order, sleep_total) %>%
    arrange(order, sleep_total) %>% 
    head
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

Same as above, except here we filter the rows for mammals that sleep for 16 or more hours instead of showing the head of the final data frame


```r
msleep %>% 
    select(name, order, sleep_total) %>%
    arrange(order, sleep_total) %>% 
    filter(sleep_total >= 16)
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

Something slightly more complicated: same as above, except arrange the rows in the sleep\_total column in a descending order. For this, use the function `desc()`


```r
msleep %>% 
    select(name, order, sleep_total) %>%
    arrange(order, desc(sleep_total)) %>% 
    filter(sleep_total >= 16)
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```


### Create new columns using `mutate()`

The `mutate()` function will add new columns to the data frame. Create a new column called rem_proportion which is the ratio of rem sleep to total amount of sleep. 



```r
msleep %>% 
    mutate(rem_proportion = sleep_rem / sleep_total) %>%
    head
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

You can many new columns using mutate (separated by commas). Here we add a second column called bodywt_grams which is the bodywt column in grams. 


```r
msleep %>% 
    mutate(rem_proportion = sleep_rem / sleep_total, 
           bodywt_grams = bodywt * 1000) %>%
    head
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

### Create summaries of the data frame using `summarise()`

The `summarise()` function will create summary statistics for a given column in the data frame such as finding the mean. For example, to compute the average number of hours of sleep, apply the `mean()` function to the column sleep\_total and call the summary value avg\_sleep. 


```r
msleep %>% 
    summarise(avg_sleep = mean(sleep_total))
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

There are many other summary statistics you could consider such `sd()`, `min()`, `max()`, `median()`, `sum()`, `n()` (returns the length of vector), `first()` (returns first value in vector), `last()` (returns last value in vector) and `n_distinct()` (number of distinct values in vector). 


```r
msleep %>% 
    summarise(avg_sleep = mean(sleep_total), 
              min_sleep = min(sleep_total),
              max_sleep = max(sleep_total),
              total = n())
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

    
### Group operations using `group_by()`

The `group_by()` verb is an important function in dplyr. As we mentioned before it's related to concept of "split-apply-combine". We literally want to split the data frame by some variable (e.g. taxonomic order), apply a function to the individual data frames and then combine the output.   

Let's do that: split the msleep data frame by the taxonomic order, then ask for the same summary statistics as above. We expect a set of summary statistics for each taxonomic order. 


```r
msleep %>% 
    group_by(order) %>%
    summarise(avg_sleep = mean(sleep_total), 
              min_sleep = min(sleep_total), 
              max_sleep = max(sleep_total),
              total = n())
```

```
## Error in eval(expr, envir, enclos): could not find function "%>%"
```

