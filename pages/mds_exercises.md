---
Title: MDS exercises
---

{pagebreak} 
 
A>## Exercises
A>
A>1. Using the `z` we computed in exercise 4 of the previous exercises:
A>
A>    
    ```r
    library(tissuesGeneExpression)
    data(tissuesGeneExpression)
    y = e - rowMeans(e)
    s = svd(y)
    z = s$d * t(s$v)
    ```
A>    
A>    we can make an mds plot
A>
A>    
    ```r
    library(rafalib)
    ftissue = factor(tissue)
    mypar2(1,1)
    plot(z[1,],z[2,],col=as.numeric(ftissue))
    legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
    ```
A>
A>    Now run the function `cmdscale` on the original data
A>
A>    
    ```r
    d = dist(t(e))
    mds = cmdscale(d)
    ```
A>
A>    What is the absolute value of the correlation between the first dimension of `z` and the first dimension in mds?
A>
A>
A>2. What is the absolute value of the  correlation between the second dimension of `z` and the second dimension in mds?
A>
A>
A>3. Load the following dataset
A>
A>    
    ```r
    library(GSE5859Subset)
    data(GSE5859Subset)
    ```
A>
A>    Compute the svd and compute `z`
A>
A>    
    ```r
    s = svd(geneExpression-rowMeans(geneExpression))
    z = s$d * t(s$v)
    ```
A>
A>    Which dimension of `z` most correlates with the outcome `sampleInfo$group`
A>
A>
A>4. What is this max correlation?
A>
A>
A>5. Which dimension of `z` has the second highest correlates with the outcome `sampleInfo$group`?
A>
A>
A>6. Note these measurements were made during two months:
A>
A>    
    ```r
    sampleInfo$date
    ```
A>
A>    We can extract the month this way:
A>    
    ```r
    month = format( sampleInfo$date, "%m")
    month = factor( month)
    ```
A>
A>    Which dimension of `z` has the second highest correlates with the outcome `month`
A>
A>
A>
A>7. What is this correlation?
A>
A>
A>8. (Advanced) Note that the same dimension is correlated with both the group and the date. Not also that these are correlated:
A>
A>    
    ```r
    table(sampleInfo$g,month)
    ```
A>
A>    So is this first dimension related directly to group or is it related only through the month? Note that the correlation with month is higher. This is related to _batch effects_ which we will learn about later.
A>
A>
A>    In exercise 3 we saw that one the  5th dimension was highly correlate to the `sampleInfo$group`. Now take the 5th column of {$$}\mathbf{U}{/$$} and stratify by the gene chromosome. Remove `chrUn` and make a boxplot of the values of {$$}\mathbf{U}_5{/$$} stratified by chromosome. 
A>
A>    Which chromosome looks different from the rest? Copy and paste the name as it appears in `geneAnnotation`
A>Explanation:
A>
A>
A>Given the answer to the last exercise, any guesses as to what `sampleInfo$group` represents?
