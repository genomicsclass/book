---
Title: Counfounding Exercises
---

## Exercises

Load the admissions data from the `dagdata` package (which is avaialbe from the genomicsclass repository):


```r
library(dagdata) 
data(admissions)
```

Familiarize yourself with this table:


```r
admissions
```

1. Let's compute the proportion of men who were accepted:

    
    ```r
    index = which(admissions$Gender==1)
    accepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
    applied = sum(admissions$Number[index])
    accepted/applied
    ```

    What is the proportion of women that were accepted?


2. Now that we have observed different acceptance rates between genders, test for the significance of this result.

    If you perform an independence test, what is the p-value?



    This difference actually led to a [lawsuit](http://en.wikipedia.org/wiki/Simpson%27s_paradox#Berkeley_gender_bias_case).

    Now notice that looking at the data by major, the differences disappear. 

    
    ```r
    admissions
    ```

    How can this be? This is referred to as Simpson's Paradox. In the following questions we will try to decipher why this is happening.

3. We can quantify how "hard" a major is using the percent of students that were accepted. Compute the percent that were accepted (regardless of gender) to each major and call this vector `H`

    Which is the hardest major? 




4. What proportion gets in for this major?


5. For men, what is the correlation between the number of applications across majors and `H`


6. For women, what is the correlation between the number of applications across majors and `H`


7. Given the answers to the above, which best explains the differences in admission percentages when we combine majors
    - A) We made a coding mistake when computing the overall admissions percentages.
    - B) There were more total number of women applications which made the denominator much bigger.
    - C) There is confounding between gender and preference for "hard" majors: females are more likely to apply to harder majors.
    - D) The sample size for the individual majors was not large enough to draw the correct conclusion.


