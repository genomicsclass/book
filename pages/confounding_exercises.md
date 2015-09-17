---
Title: Counfounding Exercises
---

A>## Exercises
A>
A>Load the admissions data from the `dagdata` package (which is avaialbe from the genomicsclass repository):
A>
A>
```r
library(dagdata) 
data(admissions)
```
A>
A>Familiarize yourself with this table:
A>
A>
```r
admissions
```
A>
A>1. Let's compute the proportion of men who were accepted:
A>
A>    
    ```r
    index = which(admissions$Gender==1)
    accepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
    applied = sum(admissions$Number[index])
    accepted/applied
    ```
A>
A>    What is the proportion of women that were accepted?
A>
A>
A>2. Now that we have observed different acceptance rates between genders, test for the significance of this result.
A>
A>    If you perform an independence test, what is the p-value?
A>
A>
A>
A>    This difference actually led to a [lawsuit](http://en.wikipedia.org/wiki/Simpson%27s_paradox#Berkeley_gender_bias_case).
A>
A>    Now notice that looking at the data by major, the differences disappear. 
A>
A>    
    ```r
    admissions
    ```
A>
A>    How can this be? This is referred to as Simpson's Paradox. In the following questions we will try to decipher why this is happening.
A>
A>3. We can quantify how "hard" a major is using the percent of students that were accepted. Compute the percent that were accepted (regardless of gender) to each major and call this vector `H`
A>
A>    Which is the hardest major? 
A>
A>
A>
A>
A>4. What proportion gets in for this major?
A>
A>
A>5. For men, what is the correlation between the number of applications across majors and `H`
A>
A>
A>6. For women, what is the correlation between the number of applications across majors and `H`
A>
A>
A>7. Given the answers to the above, which best explains the differences in admission percentages when we combine majors
A>    - A) We made a coding mistake when computing the overall admissions percentages.
A>    - B) There were more total number of women applications which made the denominator much bigger.
A>    - C) There is confounding between gender and preference for "hard" majors: females are more likely to apply to harder majors.
A>    - D) The sample size for the individual majors was not large enough to draw the correct conclusion.
A>
A>
