---
Title: Bayes Exercises
---

A>## Exercises
A>
A>
A>1. A test for cystic fibrosis has an accuracy of 99%. Specifically we mean that 
A>
A>    {$$}\mbox{Prob}(+|D)=0.99, \mbox{Prob}(-|\mbox{no } D)=0.99{/$$}
A>
A>    The cystic fibrosis rate in the general population is 1 in 3,900, {$$}\mbox{Prob}(D)=0.00025{/$$}
A>
A>    If we select random person and they test postive what is probability that they have cysttic fibrosos {$$}\mbox{Prob}(D|+){/$$} ? Hint: use Bayes Rule: 
A>
A>    {$$}
A>    \mbox{Pr}(A|B)  =  \frac{\mbox{Pr}(B|A)\mbox{Pr}(A)}{\mbox{Pr}(B)}
A>    {/$$}
A>
A>
A>
A>
A>2. (Advanced) First download some basaball statistics
A>
A>    
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
A>
A>    We will use the `dplyr`, which you can read about  [here](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) to 
A>obtain data from 2010, 2011, and 2012 and with more than 500 at bats (AB >= 500).
A>
A>    
    ```r
    dat <- filter(players,yearID>=2010, yearID <=2012) %>% mutate(AVG=H/AB) %>% filter(AB>500)
    ```
A>
A>    What is the a averge of these batting averges?
A>
A>
A>3. What is the standard deviation of these batting averages?
A>
A>
A>4. Use exploratory data analysis to decide which of the following distributions approximates our AVG
A>    - A) Normal
A>    - B) Poisson
A>    - C) F-distribution
A>    - D) Uniform
A>
A>
A>5. It is April and after 20 at bats, Jose Iglesias is batting .450 (this is very good). We can think of this as a binomial distribution with 20 trial with probability of success {$$}p{/$$}. Our sample estimate of {$$}p{/$$} is .450. What is our estimate of standard devitaion? Hint: This is sum that is binomial divided by 20
A>
A>
A>
A>6. The Binomial is approximated by normal so our sampling distribution is approximately normal wih mean {$$}Y=0.45{/$$} and SD {$$}\sigma=0.11{/$$}. Earlier we used a baseball database to determine that our prior distribution is Normal with mean {$$}\mu=0.275{/$$} and SD {$$}\tau=0.027{/$$}. We also saw that this is the posterior mean prediction of the batting average. 
A>
A>    What is your Bayes prediction for the batting average going forward.
A>   
A>    {$$}
A>\begin{align*}
A>\mbox{E}(\theta|y) &= B \mu + (1-B) Y\\
A>&= \mu + (1-B)(Y-\mu)\\
A>B &= \frac{\sigma^2}{\sigma^2+\tau^2}
A>\end{align*}
A>    {/$$}
A>
A>
