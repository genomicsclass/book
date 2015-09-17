---
Title: Distance exercises
---

A>## Exercises
A>
A>
A>If you have not done so already, install the data package `tissueGeneExpression`: 
A>
A>
```r
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
```
A>
A>The data represents RNA expression levels for eight tissues, each with several _biological replictes_. We call samples that we consider to be from the same population, such as liver tissue from different individuals, _biological replicates_: 
A>
A>
```r
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)
```
A>
A>1. How many biological replicates for hippocampus?
A>
A>
A>
A>2. What is the distance between samples 3 and 45?
A>
A>
A>
A>3. What is the distance between gene `210486_at` and `200805_at`
A>
A>
A>
A>4. If I run the command (don't run it!)
A>
A>    
    ```r
    d = as.matrix( dist(e) )
    ```
A>
A>    how many cells (number of rows times number of columns) with this matrix have?
A>
A>
A>
A>
A>5. Compute the distance between all pair of samples:
A>
A>    
    ```r
    d = dist( t(e) )
    ```
A>
A>    Read the help file for `dist`.
A>
A>    How many distances are stored in `d`? (Hint: What is the length of d)? 
A>
A>
A>
A>
A>6. Why is the answer to 1.1.4 not `ncol(e)^2`?
A>    - A) R made a mistake there
A>    - B) Distances of 0 are left out
A>    - C) Because we take advantage of symmetry: only lower triangular matrix is stored thus only `ncol(e)*(ncol(e)-1)/2` values.
A>    -D.  Because it is equal`nrow(e)^2`
A>
A>
A>
