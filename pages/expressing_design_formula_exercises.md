---
layout: page
title: Expressing Design Formula Exercises
---

## Exercises

Suppose we have an experiment with the following design: on three different days, we perform an experiment with two treated and two control units. We then measure some outcome $$Y_i$$, and we want to test the effect of treatment as well the effects of different days (perhaps the temperature in the lab affects the measuring device). Assume that the true condition effect is the same for each day (no interaction between condition and day). We then define factors in R for `day` and for `condition`.


|condition/day |  A |  B  | C|
|---|---|---|---|---|
|treatment    |  2 |   2 |   2 |
|control    |  2 |   2 |  2 |

1. Given the factors we have defined above and without defining any new ones, which of the following R formula will produce a design matrix (model matrix) that lets us analyze the effect of condition, controlling for the different days?
    - A) `~ day + condition` 
    - B) `~ condition  ~ day` 
    - C) `~ A + B + C + control + treated`  
    - D) `~ B + C + treated`
    


Remember that using the `~` and the names for the two variables we want in the model will produce a design matrix controlling for all levels of day and all levels of condition. We do not use the levels in the design formula.
