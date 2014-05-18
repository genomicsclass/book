---
layout: page
title: Downloading data from GEO using GEOquery
---




## Example of how to download CEL files from GEO

## contributed by Stephanie Hicks

If the `GEOquery` R/Biocondcutor package is not installed, use `biocLite()` to install the package:  

```r
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
```


Load the `GEOquery` R/Bioconductor package: 

```r
library(GEOquery)
```



### Access the GEO Series Data
To access the GEO Sample (GSM), GEO Series (GSE) (lists of GSM files that together form a single experiment) or GEO Dataset (GDS), use the function `getGEO()` which returns a list of ExpressionSets: 

```r
### This will download a 20 Mb
gse <- getGEO("GSE21653", GSEMatrix = TRUE)
show(gse)
```



### Accessing raw data from GEO
If raw data such as .CEL files exist on GEO, you can easily access this dea using the `getGEOSuppFiles()` function.  The function takes in a GEO accession as the argument and will download all the raw data associated with that accession. By default the `getGEOSuppFiles()` function will create a directory within the current working directory to store the raw data.  Here, the file paths of the downloaded files (often with as a .tar extension) are stored in a data frame called `filePaths`. 


```r
filePaths = getGEOSuppFiles("GSE21653")
filePaths
```

From here you can use, for example, `ReadAffy()` to read in the CEL files.  


### Access GSE Data Tables from GEO
To access the phenotypic information about the samples, the best way is to use `getGEO()` function to obtain the GSE object and then extract the phenoData object from that.  Unfortunately this means downloadint the entire GSE Matrix file.  


```r
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])
```


Sometimes GSEs are include separate data tables with the sample information. If these exist, you can uuse the `getGSEDataTables()` function. For example here is the phenoData object from a different GSE accession GSE3494 with a Data Table. 

```r
df1 <- getGSEDataTables("GSE3494")
lapply(df1, head)
```

```
## [[1]]
##   INDEX (ID) p53 seq mut status (p53+=mutant; p53-=wt)
## 1    X101B88                                      p53+
## 2    X102B06                                      p53+
## 3    X104B91                                      p53+
## 4    X110B34                                      p53+
## 5    X111B51                                      p53+
## 6    X127B00                                      p53+
##   p53 DLDA classifier result (0=wt-like, 1=mt-like)
## 1                                                 1
## 2                                                 1
## 3                                                 0
## 4                                                 1
## 5                                                 1
## 6                                                 1
##   DLDA error (1=yes, 0=no) Elston histologic grade ER status PgR status
## 1                        0                      G3       ER-       PgR-
## 2                        0                      G3       ER+       PgR+
## 3                        1                      G3       ER+       PgR+
## 4                        0                      G2       ER+       PgR+
## 5                        0                      G3       ER+       PgR+
## 6                        0                      G3       ER+       PgR+
##   age at diagnosis tumor size (mm) Lymph node status
## 1               40              12               LN-
## 2               51              26               LN-
## 3               80              24               LN?
## 4               74              20               LN-
## 5               41              33               LN-
## 6               57              22               LN-
##   DSS TIME (Disease-Specific Survival Time in years)
## 1                                             11.833
## 2                                             11.833
## 3                                              3.583
## 4                                             11.667
## 5                                              7.167
## 6                                              4.667
##   DSS EVENT (Disease-Specific Survival EVENT; 1=death from breast cancer, 0=alive or censored )
## 1                                                                                             0
## 2                                                                                             0
## 3                                                                                             0
## 4                                                                                             0
## 5                                                                                             1
## 6                                                                                             1
## 
## [[2]]
##   GEO Sample Accession # Patient ID Affy platform
## 1               GSM79114    X100B08      HG-U133A
## 2               GSM79115    X101B88      HG-U133A
## 3               GSM79116    X102B06      HG-U133A
## 4               GSM79117    X103B41      HG-U133A
## 5               GSM79118    X104B91      HG-U133A
## 6               GSM79119    X105B13      HG-U133A
```

