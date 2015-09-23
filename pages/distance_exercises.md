---
Title: Distance exercises
---

## Exercises


If you have not done so already, install the data package `tissueGeneExpression`: 


```r
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
```

The data represents RNA expression levels for eight tissues, each with several _biological replictes_. We call samples that we consider to be from the same population, such as liver tissue from different individuals, _biological replicates_: 


```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)
```

1. How many biological replicates for hippocampus?



2. What is the distance between samples 3 and 45?



3. What is the distance between gene `210486_at` and `200805_at`



4. If I run the command (don't run it!):

    
    ```r
    d = as.matrix( dist(e) )
    ```

    how many cells (number of rows times number of columns) will this matrix have?




5. Compute the distance between all pair of samples:

    
    ```r
    d = dist( t(e) )
    ```

    Read the help file for `dist`.

    How many distances are stored in `d`? Hint: What is the length of d?




6. Why is the answer to exercise 5 not `ncol(e)^2`?
    - A) R made a mistake there.
    - B) Distances of 0 are left out.
    - C) Because we take advantage of symmetry: only lower triangular matrix is stored thus only `ncol(e)*(ncol(e)-1)/2` values.
    - D)  Because it is equal`nrow(e)^2`



