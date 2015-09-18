
## Exercises

We will use the following dataset to demonstrate the use of permutations. 

```r
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
```

1. We will generate the following random variable based on a sample size of 10 and observe the following difference:

    
    ```r
    N=10
    set.seed(1)
    nonsmokers <- sample(bwt.nonsmoke , N)
    smokers <- sample(bwt.smoke , N)
    obs <- mean(smokers) - mean(nonsmokers)
    ```
  
    The question is weather this observed difference is statisticall significant. We do not want to rely on the assumptions needed for the normal or t-distribution approximations to hold. So instead we will use permutations. What we will do is reshuffle the data and recompute the mean. We can create one permuted sample with the following code:

    
    ```r
    dat <- c(smokers,nonsmokers)
    shuffle <- sample( dat )
    smokersstar <- shuffle[1:N]
    nonsmokersstar <- shuffle[(N+1):(2*N)]
    mean(smokersstar)-mean(nonsmokersstar)
    ```

    The last value is one observation from the null distribution we will construct. Set the seed at 1, and then repeat the permutation 1,000 times to create a null distribution. What is the permutation derived p-value for our observation?



2. Repeat the above exercise but instead of the differences in mean, consider the differences in median `obs <- median(smokers) - median(nonsmokers)`. What is the permutation based p-value?


