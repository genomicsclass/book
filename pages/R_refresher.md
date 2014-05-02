---
layout: page
title: R refresher
---


## Data summaries: summary, str

First we load an example data frame:


```r
rats <- data.frame(id = paste0("rat", 1:10), sex = factor(rep(c("female", "male"), 
    each = 5)), weight = c(2, 4, 1, 11, 18, 12, 7, 12, 19, 20), length = c(100, 
    105, 115, 130, 95, 150, 165, 180, 190, 175))
rats
```

```
##       id    sex weight length
## 1   rat1 female      2    100
## 2   rat2 female      4    105
## 3   rat3 female      1    115
## 4   rat4 female     11    130
## 5   rat5 female     18     95
## 6   rat6   male     12    150
## 7   rat7   male      7    165
## 8   rat8   male     12    180
## 9   rat9   male     19    190
## 10 rat10   male     20    175
```


The `summary` and `str` functions are two helpful functions for getting a sense of data. `summary` works on vectors or matrix-like objects (including data.frames). `str` works on an arbitrary R object and will compactly display the structure.


```r
summary(rats)
```

```
##        id        sex        weight          length   
##  rat1   :1   female:5   Min.   : 1.00   Min.   : 95  
##  rat10  :1   male  :5   1st Qu.: 4.75   1st Qu.:108  
##  rat2   :1              Median :11.50   Median :140  
##  rat3   :1              Mean   :10.60   Mean   :140  
##  rat4   :1              3rd Qu.:16.50   3rd Qu.:172  
##  rat5   :1              Max.   :20.00   Max.   :190  
##  (Other):4
```

```r
summary(rats$weight)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1.00    4.75   11.50   10.60   16.50   20.00
```

```r
str(rats)
```

```
## 'data.frame':	10 obs. of  4 variables:
##  $ id    : Factor w/ 10 levels "rat1","rat10",..: 1 3 4 5 6 7 8 9 10 2
##  $ sex   : Factor w/ 2 levels "female","male": 1 1 1 1 1 2 2 2 2 2
##  $ weight: num  2 4 1 11 18 12 7 12 19 20
##  $ length: num  100 105 115 130 95 150 165 180 190 175
```


## Aligning two objects: match, merge

We load another example data frame, with the original ID and another secretID. Suppose we want to sort the original data frame by the secretID.


```r
ratsTable <- data.frame(id = paste0("rat", c(6, 9, 7, 3, 5, 1, 10, 4, 8, 2)), 
    secretID = 1:10)
ratsTable
```

```
##       id secretID
## 1   rat6        1
## 2   rat9        2
## 3   rat7        3
## 4   rat3        4
## 5   rat5        5
## 6   rat1        6
## 7  rat10        7
## 8   rat4        8
## 9   rat8        9
## 10  rat2       10
```

```r
# wrong!
cbind(rats, ratsTable)
```

```
##       id    sex weight length    id secretID
## 1   rat1 female      2    100  rat6        1
## 2   rat2 female      4    105  rat9        2
## 3   rat3 female      1    115  rat7        3
## 4   rat4 female     11    130  rat3        4
## 5   rat5 female     18     95  rat5        5
## 6   rat6   male     12    150  rat1        6
## 7   rat7   male      7    165 rat10        7
## 8   rat8   male     12    180  rat4        8
## 9   rat9   male     19    190  rat8        9
## 10 rat10   male     20    175  rat2       10
```


`match` is a very useful function in R, which can give us this order, but it's easy to get its arguments mixed up. Remember that `match` gives you, for each element in the first vector, the index of the first match in the second vector. So typically the data.frame or vector you are reordering would appear as the second argument to `match`. It's always a good idea to check that you got it right, which you can do by using `cbind` to line up both data frames.


```r
match(ratsTable$id, rats$id)
```

```
##  [1]  6  9  7  3  5  1 10  4  8  2
```

```r
rats[match(ratsTable$id, rats$id), ]
```

```
##       id    sex weight length
## 6   rat6   male     12    150
## 9   rat9   male     19    190
## 7   rat7   male      7    165
## 3   rat3 female      1    115
## 5   rat5 female     18     95
## 1   rat1 female      2    100
## 10 rat10   male     20    175
## 4   rat4 female     11    130
## 8   rat8   male     12    180
## 2   rat2 female      4    105
```

```r
cbind(rats[match(ratsTable$id, rats$id), ], ratsTable)
```

```
##       id    sex weight length    id secretID
## 6   rat6   male     12    150  rat6        1
## 9   rat9   male     19    190  rat9        2
## 7   rat7   male      7    165  rat7        3
## 3   rat3 female      1    115  rat3        4
## 5   rat5 female     18     95  rat5        5
## 1   rat1 female      2    100  rat1        6
## 10 rat10   male     20    175 rat10        7
## 4   rat4 female     11    130  rat4        8
## 8   rat8   male     12    180  rat8        9
## 2   rat2 female      4    105  rat2       10
```


Or you can use the `merge` function which will handle everything for you. You can tell it the names of the columns to merge on, or it will look for columns with the same name.


```r
ratsMerged <- merge(rats, ratsTable, by.x = "id", by.y = "id")
ratsMerged[order(ratsMerged$secretID), ]
```

```
##       id    sex weight length secretID
## 7   rat6   male     12    150        1
## 10  rat9   male     19    190        2
## 8   rat7   male      7    165        3
## 4   rat3 female      1    115        4
## 6   rat5 female     18     95        5
## 1   rat1 female      2    100        6
## 2  rat10   male     20    175        7
## 5   rat4 female     11    130        8
## 9   rat8   male     12    180        9
## 3   rat2 female      4    105       10
```


## Analysis over groups: split, tapply, and dplyr libary

Suppose we need to calculate the average rat weight for each sex. We could start by splitting the weight vector into a list of weight vectors divided by sex. `split` is a useful function for breaking up a vector into groups defined by a second vector, typically a factor. We can then use the `lapply` function to calculate the average of each element of the list, which are vectors of weights.


```r
sp <- split(rats$weight, rats$sex)
sp
```

```
## $female
## [1]  2  4  1 11 18
## 
## $male
## [1] 12  7 12 19 20
```

```r
lapply(sp, mean)
```

```
## $female
## [1] 7.2
## 
## $male
## [1] 14
```


A shortcut for this is to use `tapply` and give the function which should run on each element of the list as a third argument:


```r
tapply(rats$weight, rats$sex, mean)
```

```
## female   male 
##    7.2   14.0
```


R is constantly being developed in the form of add-on packages, which can sometimes greatly simplify basic analysis tasks. A new library "dplyr" can accomplish the same task as above, and can be extended to many other more complicated operations. The "d" in the name is for data.frame, and the "ply" is because the library attempts to simplify tasks typically used by the set of functions: `sapply`, `lapply`, `tapply`, etc. Here is the same task as before done with the dplyr functions `group_by` and `summarise`:


```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
sexes <- group_by(rats, sex)
summarise(sexes, ave = mean(weight))
```

```
## Source: local data frame [2 x 2]
## 
##      sex  ave
## 1 female  7.2
## 2   male 14.0
```


With dplyr, you can chain operations using the `%.%` operator:


```r
rats %.% group_by(sex) %.% summarise(ave = mean(weight))
```

```
## Source: local data frame [2 x 2]
## 
##      sex  ave
## 1 female  7.2
## 2   male 14.0
```

