---
layout: page
title: Expressing design formula in R
---




# Expressing experimental designs using R formula

In this module, we will show how to use the two base R functions:

- `formula`
- `model.matrix`

...in order to produce design matrices for a variety of linear models.


```r
x <- c(1, 1, 2, 2)
f <- formula(~x)
f
```

```
## ~x
```



```r
model.matrix(f)
```

```
##   (Intercept) x
## 1           1 1
## 2           1 1
## 3           1 2
## 4           1 2
## attr(,"assign")
## [1] 0 1
```



```r
x <- factor(c(1, 1, 2, 2))
model.matrix(~x)
```

```
##   (Intercept) x2
## 1           1  0
## 2           1  0
## 3           1  1
## 4           1  1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
```



```r
x <- factor(c(1, 1, 2, 2, 3, 3))
model.matrix(~x)
```

```
##   (Intercept) x2 x3
## 1           1  0  0
## 2           1  0  0
## 3           1  1  0
## 4           1  1  0
## 5           1  0  1
## 6           1  0  1
## attr(,"assign")
## [1] 0 1 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
```

```r
model.matrix(~x, contrasts = list(x = "contr.sum"))
```

```
##   (Intercept) x1 x2
## 1           1  1  0
## 2           1  1  0
## 3           1  0  1
## 4           1  0  1
## 5           1 -1 -1
## 6           1 -1 -1
## attr(,"assign")
## [1] 0 1 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.sum"
```

```r
# ?contr.sum
```



```r
x <- factor(c(1, 1, 1, 1, 2, 2, 2, 2))
y <- factor(c("a", "a", "b", "b", "a", "a", "b", "b"))
model.matrix(~x + y)
```

```
##   (Intercept) x2 yb
## 1           1  0  0
## 2           1  0  0
## 3           1  0  1
## 4           1  0  1
## 5           1  1  0
## 6           1  1  0
## 7           1  1  1
## 8           1  1  1
## attr(,"assign")
## [1] 0 1 2
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
## 
## attr(,"contrasts")$y
## [1] "contr.treatment"
```



```r
model.matrix(~x + y + x:y)
```

```
##   (Intercept) x2 yb x2:yb
## 1           1  0  0     0
## 2           1  0  0     0
## 3           1  0  1     0
## 4           1  0  1     0
## 5           1  1  0     0
## 6           1  1  0     0
## 7           1  1  1     1
## 8           1  1  1     1
## attr(,"assign")
## [1] 0 1 2 3
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
## 
## attr(,"contrasts")$y
## [1] "contr.treatment"
```

```r
model.matrix(~x * y)
```

```
##   (Intercept) x2 yb x2:yb
## 1           1  0  0     0
## 2           1  0  0     0
## 3           1  0  1     0
## 4           1  0  1     0
## 5           1  1  0     0
## 6           1  1  0     0
## 7           1  1  1     1
## 8           1  1  1     1
## attr(,"assign")
## [1] 0 1 2 3
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
## 
## attr(,"contrasts")$y
## [1] "contr.treatment"
```



```r
x <- factor(c(1, 1, 2, 2))
model.matrix(~x)
```

```
##   (Intercept) x2
## 1           1  0
## 2           1  0
## 3           1  1
## 4           1  1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
```

```r
x <- relevel(x, "2")
model.matrix(~x)
```

```
##   (Intercept) x1
## 1           1  1
## 2           1  1
## 3           1  0
## 4           1  0
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"
```

```r
x <- factor(x, levels = c("1", "2"))
```



```r
z <- 1:4
model.matrix(~z)
```

```
##   (Intercept) z
## 1           1 1
## 2           1 2
## 3           1 3
## 4           1 4
## attr(,"assign")
## [1] 0 1
```

```r
model.matrix(~0 + z)
```

```
##   z
## 1 1
## 2 2
## 3 3
## 4 4
## attr(,"assign")
## [1] 1
```

```r
model.matrix(~z + I(z^2))
```

```
##   (Intercept) z I(z^2)
## 1           1 1      1
## 2           1 2      4
## 3           1 3      9
## 4           1 4     16
## attr(,"assign")
## [1] 0 1 2
```



```r
x <- 1:4
model.matrix(~x)
```

```
##   (Intercept) x
## 1           1 1
## 2           1 2
## 3           1 3
## 4           1 4
## attr(,"assign")
## [1] 0 1
```

```r
model.matrix(~x, data = data.frame(x = 5:8))
```

```
##   (Intercept) x
## 1           1 5
## 2           1 6
## 3           1 7
## 4           1 8
## attr(,"assign")
## [1] 0 1
```



