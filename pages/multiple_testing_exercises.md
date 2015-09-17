---
Title: Multiple Testing Exercises
---

A>## Exercises
A>
A>With these exercises we hope to help you further grasp the concept that p-values are random variables and start laying the ground work for the development of procedures that control error rates. The calculations to compute error rates require us to understand the random behavior of p-values.
A>
A>We are going to ask you to perform some calculations related to introductory probability theory. One particular concept you need to grasp is statistical independence. You also will need to know that the probability of two random events that are statistically independent occurring is {$$}P( A \mbox{ and } B) = P(A)P(B){/$$}. Note that this is a consequence of the more general formula {$$}P(A \mbox{ and } B) = P(A) P(B | A ){/$$}.
A>
A>1. Assume the null is true and denote the p-value you would get if you ran a test as {$$}P{/$$}. Define the function {$$}f(x) = \mbox{Pr}(P>x){/$$} What does {$$}f(x){/$$} look like?
A>    - A) A uniform distribution.
A>    - B) The identity line.
A>    - C) A constant at 0.05.
A>    - D) P is not a random value
A>
A>
A>
A>2. In the previous exercises we saw how the probability of incorrectly rejecting the null for at least one of 20 experiments for which the null is true is well over 5%. Now let's consider a case in which we run thousands of tests as we would do in a high throughput experiment. 
A>
A>    We previously learned that under the null, the probability of a p-value < p is p. If we run 8,793 independent tests, what it the probability of incorrectly rejecting at least one of the null hypothesis? 
A>
A>
A>
A>
A>
A>3. Suppose we need to run 8,793 statistical tests and we want to make the probability of a mistake very small, say 5%. Use the answer to Exercise n 2 to determine how small we have to change the cutoff, previously 0.05, to lower our probability of at least one mistake to be 5%. 
A>
A>
A>
A>    The following exercises should help you understand the concept of an error controlling procedure. You can think of it as defnining a set of instructions, such as "reject all the null hypothesis for  for which p-values < 0.0001" or "reject the null hypothesis for the 10 features with smallest p-values". Then, knowing the p-values are random variables, we use statistical theory to compute how many mistakes, on average, will we make if we follow this procedure. More precisely we commonly bounds on these rates, meaning that we show that they are smaller than some predermined value.
A>
A>    As described in the video, we can compute different error rates. The FWER tells us the probability of having at least one false positive. The FDR is the expected rate of rejected null hypothesis.
A>
A>    Note 1: the FWER and FDR are not procedures but error rates. We will review procedures here and use Monte Carlo simulations to estimate their error rates.
A>
A>    Note 2: We sometimes use the colloquial term "pick genes that" meaning "reject the null hypothesis for genes that".
A>
A>4. So we have learned about the family wide error rate FWER. This is the probability of incorrectly rejecting the null at least once. Using the notation in the video this probability is written like this: {$$}\mbox{Pr}(V>0){/$$}. 
A>
A>    What we want to do in practice is choose a _procedure_ that guarantees this probability is smaller than a predetermined value such as 0.05. Here we keep it general and instead of 0.05 we use {$$}\alpha{/$$} 
A>
A>    We have already learned that the procedure "pick all the genes with p-value <0.05" fails  miserably as we have seen that {$$}Pr(V>0) \approx 1{/$$}. So what else can we do?
A>
A>    The Bonferroni procedure assumes we have computed p-values for each test and asks what constant {$$}k{/$$} should we pick so that the procedure "pick all genes with p-value less than {$$}k{/$$} has {$$}\mbox{Pr}(V>0) = 0.05{/$$}. And we typically want to be conservative rather than lenient, so we accept a procedure that has {$$}\mbox{Pr}(V>0) \leq 0.05{/$$}. 
A>
A>    So the first result we rely on is that this probability is largest when all the null hypotheses are true:
A>
A>    {$$}\mbox{Pr}(V>0) \leq \mbox{Pr}(V>0 | \mbox{all nulls are true}){/$$}
A>
A>    or using the notation in the video:
A>
A>    {$$}\mbox{Pr}(V>0) \leq \mbox{Pr}(V>0 | m_1=0){/$$}
A>
A>    In an earlier assessment we showed that if the tests are independent then 
A>
A>    {$$}\mbox{Pr}(V>0 | m_1) = 1-(1-k)^m{/$$}
A>
A>    And we pick {$$}k{/$$} so that {$$}1-(1-k)^m = \alpha \implies k = 1-(1-\alpha)^{1/m}{/$$}
A>
A>    Now, this requires the tests to be independent. The Bonferroni procedure does not make this assumption and as we previously we sets {$$}k=\alpha/m{/$$} and shows that with this for this procedure {$$}Pr(V>0) \leq \alpha{/$$}. 
A>
A>
A>    In R define 
A>    
    ```r
    alphas <- seq(0,0.25,0.01)
    ```
A>  
A>    Make a plot of {$$}\alpha/m{/$$} and {$$}1-(1-\alpha)^{1/m}{/$$} for various values of m>1. 
A>
A>    Which procedure is more conservative Bonferroni's or Sidek's?
A>    - A) The are the same
A>    - B) Bonferroni's
A>    - C) Depends on m
A>    - D) Sidek's
A>
A>
A>
A>5.  To simulate the p-value results of, say a 8,792 t-tests for which the null is true we don't actual have to generate the original data. As we learned in class we can generate p-values for a uniform distribution like this: 'pvals <- runif(8793,0,1)`. 
A>
A>    Using what we have learned set the cutoff using the Bonferroni correction and report back the FWER. Set the seed at 1 and run 10,000 simulation.
A>
A>
A>
A>6. Using the same seed repeat the above for Sidek's cutoff
A>
A>
A> 
A>    In the following questions assessment we will define error controlling procedures for experimental data. We will make list of genes based on q-values. We will also assess your understanding of false positives rates and false negative rates by asking you to create a Monte Carlo simulation.
A>
A>7. Load the gene expression data
A>
A>    
    ```r
    library(GSE5859Subset)
    data(GSE5859Subset)
    ```
A>
A>    We are interested in comparing gene expression between the two groups defined in the sampleInfo table. 
A>
A>    Compute a p-value for each gene using the function rowttests from the genefilter package.
A>
A>    
    ```r
    library(genefilter)
    ?rowttests
    ```
A>
A>    How many genes have p-values smaller than 0.05?
A>
A>
A>
A>
A>8. Apply the Bonferroni correction to achieve a FWER of 0.05. How many genes are called significant under this procedure?
A>
A>
A>
A>9. Note that the FDR is a property of a list of features, not each specific feature. The q-value relates FDR to individual feature. To define the q-value we order features we tested by p-value then compute the FDRs for a list with the most significant, the two most significant, the three most significant, etc... The FDR of the list with the, say, {$$}m{/$$} most significant tests is defined as the q-value of the {$$}m{/$$}-th most significant feature. In other words, the q-value of a feature, is the FDR of the biggest list that includes that gene.
A>
A>    In R, we can compute the q-value of using the p.adjust function with the FDR option. Read the help file for p.adjust and the for our gene expression dataset compute the how many genes achieve an FDR < 0.05
A>
A>
A>
A>
A>10. Now use the qvalue function, in the Bioconductor qvalue package, to estimate q-values using the procedure described by Storey. How many genes have q-values below 0.05?
A>
A>
A>
A>11. Read the help file for qvalue and report the estimated proportion of genes for which the null hypothesis is true {$$}\pi_0=m_0/m{/$$}
A>
A>
A>
A>12. Note that we have the number of genes passing the q-value <0.05 threshold is larger with the qvalue function than the p.adjust difference. Why is this the case? Make a plot of the ratio of these two estimates to help answer the question.
A>    - A) One of the two procedures is flawed
A>    - B) The two functions are estimating different things
A>    - C) The qvalue function estimates the proportion of genes for which the null hypothesis is true and provides a less conservative estimate
A>    - D) The qvalue function estimates the proportion of genes for which the null hypothesis is true and provides a more conservative estimate
A>
A>
A>
A>    To get an idea of how {$$}\pi_0{/$$} is estimated, note that if we look at the histogram, {$$}\pi_0{/$$} roughly tells us the proportion that looks about uniform:
A>
A>
```r
hist(pvals,breaks=seq(0,1,len=21))
expectedfreq <- length(pvals)/20 #per bin
abline(h=expectedfreq*qvalue(pvals)$pi0,col=2,lty=2)
```
A>
A>
A>    
A>13. Note that this and the remaining questions are more advanced questions.
A>. Create a Monte Carlo Simulation in which you simulate measurements from 8,793 genes for 24 samples, 12 cases and 12 controls. 
A>
A>    
    ```r
    n <- 24
    m <- 8793
    mat <- matrix(rnorm(n*m),m,n)
    ```
A>
A>    Now for 100 genes, there is a difference of 1 between cases and controls:
A>
A>    
    ```r
    delta <- 1
    positives <- 500
    mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
    ```
A>
A>    So the null hypothesis is true for 8793-100 genes.
A>
A>    Run this experiment 1,000 with a Monte Carlo simulation. For each instance compute p-values using a t-test and report keep track of the number of false positives and false negatives. Compute the false positive rate and false negative rates if we use Bonferroni, q-values from p.adjust, and  q-values from qvalue function. Set the seed to 1 for each simulations. What is false positive rate for Bonferroni?
A>
A>
A>
A>14. What are the false negative rate for Bonferroni?
A>
A>
A>15. What are the false negative rate for `p.adjust`?
A>
A>
A>16. What are the false negative rate for `p.adjust`?
A>
A>
A>17. What are the false negative rate for `qvalues`?
A>
A>
A>18. What are the false negative rate for `qvalues`?
A>
A>
