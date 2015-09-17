---
Title: Introduction Exercesis
---


A>## Exercises
A>
A>Here we include some refresher questions. 
A>If you haven't done so already, install the library UsingR
A>
A>
```r
install.packages("UsingR")
```
A>
A>Then once you load it you have access to Galton's father and son heights:
A>
A>
```r
data("father.son",package="UsingR")
```
A>1. What is the average height of the sons (don't round off)?<<
A>
A>
A>
A>2. One of the defining features of regression is that we stratify one variable based on others. In Statistics we use the verb "condition". For example, the linear model for son and father heights answers the question how tall do I expect a son to be if I condition on his father being x inches. The regression line answers this question for any {$$}x{/$$}.
A>
A>    Using the `father.son` dataset described above, we want to know the expected height of sons if we condition on the father being 71 inches. Create a list of son height's for sons that have fathers with height of 71 inches (round to the nearest inch).
A>
A>    What is the mean of the son heights for fathers that have a height of 71 inches (don't round off your answer)? (Hint: use the function round() on the fathers' heights)
A>
A>
A>    
A>3. We say a statistical model is a linear model when we can write it as a linear combination of parameters and known covariates plus random error terms. In the choices below, Y represents our observations, time t is our only covariate, unknown parameters are represented with letters a,b,c,d and measurment error is represented by the letter e. Note that if t is known, then any transformation of t is also known. So, for example, both {$$}Y=a+bt +\varepsilon{/$$} and {$$}Y=a+b f(t) + \varepsilon{/$$} are linear models. Which of the following **cannot** be written as a linear model?
A>  - A) {$$}Y = a + bt + \varepsilon{/$$}
A>  - B) {$$}Y = a + b \cos(t) + \varepsilon{/$$}
A>  - C) {$$}Y = a + b^t + \varepsilon{/$$}
A>  - D) {$$}Y = a + b t + c t^2 + d t^3 + \varepsilon{/$$}
A>  
A>
A>
A>3. Supposed you model the relationship between weight and height across individuals with a linear model. You assume that the height of individuals for a fixed weight x follows a liner model {$$}Y = a + b x + \varepsilon{/$$}. Which of the following do you feel best describes what e represents?
A>  - A) Measurement error: scales are not perfect.
A>  - B) Within individual random fluctuations: you don't weigh the same in the morning as the afternoon.
A>  - C) Round off error introduced by the computer we use to analyze the data.
A>  - D) Between individual variability: people of the same height vary in their weight.
A>  
A>  
A>
