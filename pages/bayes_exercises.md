---
layout: page
title: Bayes Exercises
---

{pagebreak} 
 
## Exercises


1. A test for cystic fibrosis has an accuracy of 99%. Specifically, we mean that:

    $$\mbox{Prob}(+|D)=0.99, \mbox{Prob}(-|\mbox{no } D)=0.99$$

    The cystic fibrosis rate in the general population is 1 in 3,900, $$\mbox{Prob}(D)=0.00025$$

    If we select a random person and they test positive, what is probability that they have cystic fibrosis $$\mbox{Prob}(D|+)$$ ? Hint: use Bayes Rule. 

    $$
    \mbox{Pr}(A|B)  =  \frac{\mbox{Pr}(B|A)\mbox{Pr}(A)}{\mbox{Pr}(B)}
    $$




2. (Advanced) First download some baseball statistics.

    
    ```r
    tmpfile <- tempfile()
    tmpdir <- tempdir()
    download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
    ##this shows us files
    filenames <- unzip(tmpfile,list=TRUE)
    players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
    unlink(tmpdir)
    file.remove(tmpfile)
    ```

    We will use the `dplyr`, which you can read about [here](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) to 
obtain data from 2010, 2011, and 2012, with more than 500 at bats (AB >= 500).

    
    ```r
    dat <- filter(players,yearID>=2010, yearID <=2012) %>% mutate(AVG=H/AB) %>% filter(AB>500)
    ```

    What is the average of these batting averages?


3. What is the standard deviation of these batting averages?


4. Use exploratory data analysis to decide which of the following distributions approximates our AVG:
    - A) Normal.
    - B) Poisson.
    - C) F-distribution.
    - D) Uniform.


5. It is April and after 20 at bats, José Iglesias is batting .450 (which is very good). We can think of this as a binomial distribution with 20 trials, with probability of success $$p$$. Our sample estimate of $$p$$ is .450. What is our estimate of standard deviation? Hint: This is the sum that is binomial divided by 20.



6. The Binomial is approximated by normal, so our sampling distribution is approximately normal with mean $$Y=0.45$$ and SD $$\sigma=0.11$$. Earlier we used a baseball database to determine that our prior distribution is Normal with mean $$\mu=0.275$$ and SD $$\tau=0.027$$. We also saw that this is the posterior mean prediction of the batting average. 

    What is your Bayes prediction for the batting average going forward?
   
    $$
\begin{align*}
\mbox{E}(\theta|y) &= B \mu + (1-B) Y\\
&= \mu + (1-B)(Y-\mu)\\
B &= \frac{\sigma^2}{\sigma^2+\tau^2}
\end{align*}
    $$


