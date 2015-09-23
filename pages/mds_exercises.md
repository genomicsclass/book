---
layout: page
Title: MDS exercises
---

{pagebreak} 
 
## Exercises

1. Using the `z` we computed in exercise 4 of the previous exercises:

    
    ```r
    library(tissuesGeneExpression)
    data(tissuesGeneExpression)
    y = e - rowMeans(e)
    s = svd(y)
    z = s$d * t(s$v)
    ```
    
    we can make an mds plot:

    
    ```r
    library(rafalib)
    ftissue = factor(tissue)
    mypar2(1,1)
    plot(z[1,],z[2,],col=as.numeric(ftissue))
    legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
    ```

    Now run the function `cmdscale` on the original data:

    
    ```r
    d = dist(t(e))
    mds = cmdscale(d)
    ```

    What is the absolute value of the correlation between the first dimension of `z` and the first dimension in mds?


2. What is the absolute value of the  correlation between the second dimension of `z` and the second dimension in mds?


3. Load the following dataset:

    
    ```r
    library(GSE5859Subset)
    data(GSE5859Subset)
    ```

    Compute the svd and compute `z`.

    
    ```r
    s = svd(geneExpression-rowMeans(geneExpression))
    z = s$d * t(s$v)
    ```

    Which dimension of `z` most correlates with the outcome `sampleInfo$group` ?


4. What is this max correlation?


5. Which dimension of `z` has the second highest correlation with the outcome `sampleInfo$group`?


6. Note these measurements were made during two months:

    
    ```r
    sampleInfo$date
    ```

    We can extract the month this way:
    
    ```r
    month = format( sampleInfo$date, "%m")
    month = factor( month)
    ```

    Which dimension of `z` has the second highest correlation with the outcome `month`



7. What is this correlation?


8. (Advanced) The same dimension is correlated with both the group and the date. The following are also correlated:

    
    ```r
    table(sampleInfo$g,month)
    ```

    So is this first dimension related directly to group or is it related only through the month? Note that the correlation with month is higher. This is related to _batch effects_ which we will learn about later.


    In exercise 3 we saw that one of the dimensions was highly correlated to the `sampleInfo$group`. Now take the 5th column of $$\mathbf{U}$$ and stratify by the gene chromosome. Remove `chrUn` and make a boxplot of the values of $$\mathbf{U}_5$$ stratified by chromosome. 

    Which chromosome looks different from the rest? Copy and paste the name as it appears in `geneAnnotation`.


Given the answer to the last exercise, any guesses as to what `sampleInfo$group` represents?
