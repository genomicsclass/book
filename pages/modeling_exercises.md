---
Title: Modeling Exercises
---

## Exercises


1. Supposed you have an urn with blue and red balls. If $$N$$ balls at randaom with replacement (you put the ball back after you pick it) we can denote the outcomes as random variables $$X_1,\dots,X_N$$ that are 1 or 0. If the proportion of red ball is $$p$$ then the distribution of each of these is $$\mbox{Pr}(X_i=1)=p$$. 

    These are also called bernoulli trials. Note that these random variables are independent because we replace the balls. Flipping a coin is an example of this with $$p=0.5$$. 

    You can show that the mean and variance are $$p$$ and $$p(1-p)$$ respectively. The binomial distribution gives us the distribution of the sum $$S_N$$ of these random variables. What is the probability that we see $$k$$ red balls is given by 

    $$ \mbox{Pr}(S_N=k) = {N \choose k} p^k (1-p)^{N-k} $$

    In R the function `dbimom` gives you this result. The function `pbinom` gives us $$\mbox{Pr}(S_N\leq k)$$.

    This equation has many uses in the life sciences. We give some examples below

    The probability of conceiving a girl is 0.49. What is the probability that a familiy with 4 children has 2 girls and 2 boys (you can 


2. Use what you learned in Question 4.1.1. to answer these questions:

    What is the probability that a familiy with 10 children has 6 girls and 4 boys (you can assume no twins)?


3. The genome has 3 billion bases. About  20% are C, 20% are G, 30% are T and 30% are A. Suppose you take a random interval of 20 bases, what is the probability that the GC-content (proportion of Gs of Cs) is striclty above 0.5 in this interval.




4. The probability of winning the lottery is 1 in  175,223,510. If 20,000,000 people buy a ticket, what is the probability that more than on person wins?


5. We can show that the binomail approximation is approximately normal with $$N$$ is large and $$p$$ is not too close to 0 or 1. This means that 

    $$\frac{S_N - \mbox{E}(S_N)}{ \sqrt{ \mbox{Var}(S_N)}}  $$ 

    is approximately normal with mean 0 and SD 1. Using the results for sums of indepedent random variables we learned in a previos course, we can show that $$\mbox{E}(S_N) = Np$$ and $$\mbox{Var}(S_n)=Np(1-p)$$. 

    The genome has 3 billion bases. About  20% are C, 20% are G, 30% are T and 30% are A. Suppose you take a random interval of 20 bases, what is the exact probability that the GC-content (proportion of Gs of Cs) is greater than 0.35 and smaller or equal to 0.45 in this interval.



6. For the question above, what is the normal approximation to the probability?


    Note how close these are. 

7. Repeat question 4.1.3 but using an interval of 1000 bases. What is the difference (in absolute value) between the normal approximation and the exact distribution of the GC-content being greater than 0.35 and lesser or equal to 0.45?


8. The Cs in our genomes can be [_methylated_](http://en.wikipedia.org/wiki/DNA_methylation) or _unmethylated_. Suppose we have a large (millions) group of cells in which a proportion $$p$$ of a C of interest are methylated. We break up the DNA of these cells and randomely select pieces and end up with $$N$$ pieces that contain the C we care abut. This means that the probability of seing $$k$$ methylated Cs is binomial:

    
    ```r
    exact = dbinom(k,N,p)
    ```

    We can approximate this with the normal distribution:

    
    ```r
    a <- (k+0.5 - N*p)/sqrt(N*p*(1-p))
    b <- (k-0.5 - N*p)/sqrt(N*p*(1-p))
    approx = pnorm(a) - pnorm(b)
    ```

    Compute the difference `approx - exact` for 
    
    ```r
    N <- c(5,10,50,100,500)
    p <- seq(0,1,0.25)
    ```

    Compare the apprximation and exact probabilty of the proportion of Cs being $$p$$, $$k=1,\dots,N-1$$ plotting the exact versus approx for each $$p$$ and $$N$$ combination
    - A) The normal approximation works well when $$p$$ is close to 0.5 even for small $$N=10$$
    - B) The normal apprximation breaks down when $$p$$ is close to 0 or 1 even for large $$N$$
    - C) When $$N$$ is 100 all approximations are spot on.
    - D) When $$p=0.01$$ the approximation are terrible for $$N=5,10,30$$ and only ok for $$N=100$$


9. We saw in the previous question that when $$p$$ is very small the normal approximation breaks down. If $$N$$ is very large then we can use the Poisson approximation. 

    Earlier we computed of 1 or more people winning the lottery when the probability of winning was 1 in  175,223,510 and  20,000,000 people bought a ticket. Using the binmial we can compute the probability of exactly two people winning to be:
    
    ```r
    N <- 20000000
    p <- 1/175223510
    dbinom(2,N,p)
    ```
    
    If we were to use the normal approximation we would greatly underestimate this

    
    ```r
    a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
    b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
    pnorm(a) - pnorm(b)
    ```

    To use the Poisson approximation here use the rate $$\lambda= Np$$ representing the number of people per 20,000,000 that win the lottery. Note how much better the approximation is:

    
    ```r
    dpois(2,N*p)
    ```

    In this case it is practically the same because $$N$$ is very very large and $$Np$$ is not 0. These are the assumptions needed for the Poisson to work. 

    What is the Poisson approximation for more than one person winning?


10. Now we are going to explore if palindromes are over-represented in some part of the HCMV genome. Make sure you have the latest version of the `dagdata` library and then load the palindrom data from the Human cytomegalovirus genome:
    
    ```r
    library(dagdata)
    data(hcmv)
    ```

    These are the locations of palindromes on the genome of this virus:
    
    ```r
    library(rafalib)
    mypar()
    plot(locations,rep(1,length(locations)),ylab="",yaxt="n")
    ```

    These palindroms are quite rare, $$p$$ is very samll. If we break the genome into bins of 4000 basepairs, then we have $$Np$$ not so small and we might be able to use Poison to model the number of palindroms in each bin:

    
    ```r
    breaks=seq(0,4000*round(max(locations)/4000),4000)
    tmp=cut(locations,breaks)
    counts=as.numeric(table(tmp))
    ```

    So if our model is correct `counts` should follow a Poisson distribtion. The distribution seems about right:

    
    ```r
    hist(counts)
    ```

    Let $$X_1,\dots,X_n$$ be the random variables represneting counts then   
  
    $$
\Pr(X_i=k) = \lambda^k / k! \exp ( -\lambda)
    $$

    So to fully describe this distribtuion we need $$\lambda$$. For this we will use MLE.
To compute the Maximum Likelihood Estimate (MLE) we ask what is the probability of observing our data (which we denote with small caps):

    $$L(\lambda) = \mbox{Pr}(X_1=x_1 \mbox{ and } X_2=x2 \mbox{ and } \dots X_n=x_n)$$

    We assume that the $$X$$ are independent thus the probabilities multiply:

    $$L(\lambda) = \mbox{Pr}(X_1=x_1) \times \mbox{Pr}(X_2=x2) \times \dots \times \mbox{Pr}(X_n=x_n)$$

    Now we can write it in R. For example for $$\lambda=4$$ we have:
    
    ```r
    probs <- dpois(counts,4)
    likelihood <- prod(probs)
    likelihood
    ```

    Note it's a tiny number. It is usually more convinient to compute log-likelihoods

    
    ```r
    logprobs <- dpois(counts,4,log=TRUE)
    loglikelihood <- sum(logprobs)
    loglikelihood
    ```

    Now write a function that takes $$\lambda$$ and the vector of counts as input returns the log-likelihood. Compute this log-likelihood for `lambdas = seq(0,15,len=300)` and make a plot. 

    What value of `lambdas` maximizes the log-likelihood? 


    It turns out that we can work mathematically what $$\lambda$$ maximizes the likelihood using calculus. I turns out that the avearage of the counts is the MLE.

    
    ```r
    mean(counts)
    ```

11. The point of collecting this dataset was to try to determine if there is a region of the genome that has higher palindrome rate than expected. We can create a plot and see the counts per location:

    
    ```r
    library(dagdata)
    data(hcmv)
    breaks=seq(0,4000*round(max(locations)/4000),4000)
    tmp=cut(locations,breaks)
    counts=as.numeric(table(tmp))
    binLocation=(breaks[-1]+breaks[-length(breaks)])/2
    plot(binLocation,counts,type="l",xlab=)
    ```

    What is the center of the bin with the highest count?


    What is the maximum count? 



13. Now that we have identified the location with the largest palindrom count we want to know if we could see a value this big by chance?

    If $$X$$ is a Poisson random variable with rate

    
    ```r
    lambda = mean(counts[ - which.max(counts) ])
    ```

    What is the probability of seeing a count of 14 or more?


14. So we obtain a p-value smaller than 0.001 for a count of 14. 

    Why is it problematic to report this p-value as strong evidence of a location that is different?

    - A) Poisson in only an approximation
    - B) We selected the highest region out of 57 and need to adjust for multiple testing.
    - C) $$\lambda$$ is an estimate, a random variable, and we didn't take into account its variability
    - D)  We don't know the effect size



15. Use the Bonferonni correction to determine the p-value cut-off thtat guarantees a FWER of 0.05. What is this p-value cutoff ?



    Note that our observed p-value satistify the Bonferroni correction. 


16. Create a qq-plot to see if our Poisson model is a good fit:

    
    ```r
    ps <- (seq(along=counts) - 0.5)/length(counts)
    lambda <- mean( counts[ -which.max(counts)])
    poisq <- qpois(ps,lambda)
    plot(poisq,sort(counts))
    abline(0,1)
    ```

    How would you characterize this qq-plot
    - A) Poisson is a terrible approximation
    - B) Poisson is a very good approximation except for on point that we actually think is a region of interest
    - C) There are too many 1s in the data
    - D) A normal distribution provides a better approximation




17. Load the `tissuesGeneExpression` data library
    
    
    ```r
    library(tissuesGeneExpression)
    ```

    Now load this data and select the columns related to endometrium 

    
    ```r
    library(genefilter)
    y = e[,which(tissue=="endometrium")]
    ```

    This will give you a matrix `y` with 15 samples.

    Compute the across sample variance for the first three samples. Then make a qq-plot to see if the data follow a normal distribution. Which of the following is true?
    - A) With the exception of a handful of outliers, the data follow a normal distribution
    - B) The variance does not follow a normal distribution but taking the square root fixes this
    - C) The normal distribution is not usable here: the left tail in over estimated and the right tail is underestimated
    - D) The normal distribution fits the data almost perfectly



18. Now fit an F-distribtuion with 14 degrees of freedom using the `fitFDist` function in the `limma` package:


19. Now creat a qq-plot of the observed sample variances versus the F-distribution quantiles. Which of the following best describes the qq-plot?
    - A) The fitted F-distribution provides a perfect fit
    - B) If we exlucde the lowest 0.1% of the data, the F-distribution provides a good fit.
    - C) The normal distribution proived a better fit
    - D) If we exclude the highest 0.1% of the data, the F-distrubtion provides a good fit.
  


