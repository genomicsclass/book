---
title: EDA Exercises
---




{pagebreak}

A>## Exercises
A>
A>![ ](images/R/exploratory_data_analysis_exercises-tmp-hist_exercise-1.png) 
A>
A>1. Given the above histogram: how many people are between the ages of 35 and 45?
A>
A>
A>2. The `InsectSprays` data set is included in R. The dataset reports the counts of insects in agricultural experimental units treated with different insecticides. Make a boxplot and determine which insecticide appears to be most effective. 
A>
A>
A>
A>3.  Download and load [this](http://courses.edx.org/c4x/HarvardX/PH525.1x/asset/skew.RData) dataset into R.
A>Use exploratory data analysis tools to determine which two columns are different from the rest. What is the column that is positively skewed? 
A>
A>
A>
A>
A>
A>4. Which is the column that is negatively skewed?
A>
A>
A>5. Let's consider a random sample of finishers from the New York City Marthon in 2002.  This data set can be found in the UsingR package. Load the library and then load the nym.2002 data set. 
A>
A>    
    ```r
    library(dplyr)
    data(nym.2002, package="UsingR")
    ```
A>
A>    Use boxplots and histograms to compare the finishing times of males and females. Which of the following best describes the difference.
A>    - A) Males and females have the same distribution.
A>    - B) Most males are faster than most women.
A>    - C) Male and females have similar right skewed distributions with the former, 20 minues shifted to the left.
A>    - D) Both distribution are normally distributed with a difference in mean of about 30 minutes.
A>  
A>
A>
A>
A>6. Use `dplyr` to create two new data frames `males` and `females` with the data for each gender. For  males what is the Pearson correlation between age and time to finish? 
A>
A>
A>7. For females, what is the Pearson correlation between age and time to finish? 
A>
A>
A>8. If we interprete these correlations without visualizing the data we would conclude that the older we get the slower we run marathons, regardless of gender. Look at scatter plots and boxplots of times stratified by age groups (20-25, 25-30, etc..). After examing the data what is a more reasonable conclusion?
A>    - A) Finish times are constant up until about our 40s then we get slower
A>    - B) On average, finsih times go up by about 7 minutes every five years
A>    - C) The optimal age to run a marathon is 20-25
A>    - D) Coding errors never happen: a five year old boy completed the 2012 NY city marathon
A>
A>
A>
A>
A>9. When is it appropriate to use pie charts or donut charts?
A>  - A)  When you are hungry 
A>  - B) To compare percentages  
A>  - C) To compare values that add up to 100% 
A>  - D) Never
A>
A>
A>10. The use of pseudo-3D plots in the literature mostly adds:
A>  - A) Pizzazz 
A>  - B) The ability to see three dimensional data
A>  - C) Abilitiy to discover
A>  - D) Confusion
A>
A>
