---
Title: Expressing Design Formula Exercises
---

A>## Exercises
A>
A>Suppose we have an experiment with the following design: on three different days, we perform an experiment with two treated and two control samples. We then measure some outcome Y_i, and we want to test the effect of condition, while controlling for whatever differences might have occured due to the the different day (maybe the temperature in the lab affects the measuring device). Assume that the true condition effect is the same for each day (no interaction between condition and day). We then define factors in R for 'day' and for 'condition'.
A>
A>
A>|condition/day |  A |  B  | C|
A>|---|---|---|---|---|
A>|treated    |  2 |   2 |   2 |
A>|control    |  2 |   2 |  2 |
A>
A>1. Given the factors we have defined above, and not defining any new ones, which of the following R formula will produce a design matrix (model matrix) that let's us analyze the effect of condition, controlling for the different days:
A>    - A) `~ day + condition` 
A>    - B) `~ condition  ~ day` 
A>    - C) `~ A + B + C + control + treated`  
A>    - D) `~ B + C + treated`
A>    
A>
A>
A>Remember that using the `~` and the names for the two variables we want in the model will produce a design matrix controlling for all levels of day and all levels of condition. We do not use the levels in the design formula.
