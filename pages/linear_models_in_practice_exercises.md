---
Title: Linear Models in Practice Exercises.
---

## Exercises

The function `lm` can be used to fit a simple, two group linear model. The test statistic from a linear model is equivalent to the test statistic we get when we perform a t-test with the equal variance assumption. Though the linear model in this case is equivalent to a t-test, we will soon explore more complicated designs, where the linear model is a useful extension (confounding variables, testing contrasts of terms, testing interactions, testing many terms at once, etc.).

Here we will review the mathematics on why these produce the same test statistic and therefore p-value.

We already know that the numerator of the t-statistic in both cases is the difference between the average of the groups, so we only have to see that the denominator is the same. Of course, it makes sense that the denominator should be the same, since we are calculating the standard error of the same quantity (the difference) under the same assumptions (equal variance), but here we will show equivalence of the formula.

In the linear model, we saw how to calculate this standard error using the design matrix $$X$$ and the estimate of $$\sigma^2$$ from the residuals. The estimate of $$\sigma^2$$ was the sum of squared residuals divided by $$N - p$$, where $$N$$ is the total number of samples and $$p$$ is the number of terms (an intercept and a group indicator, so here $$p=2$$).

In the t-test, the denominator of the t-value is the standard error of the difference. The t-test formula for the standard error of the difference, if we assume equal variance in the two groups, is the square root of the variance:

$$ \frac{1}{1/N_x + 1/N_y}  
\frac{  \sum_{i=1}^{N_x} (X_i - \mu_x)^2  + \sum_{i=1} (Y_i - \mu_y)^2  }{N_x + N_y - 2}
$$


Here $$N_x$$ is the number of samples in the first group and $$N_y$$ is the number of samples in the second group.

If we look carefully, the second part of this equation is the sum of squared residuals, divided by $$N - 2$$.

All that is left to show is that the entry in the second row, second column of $$(X^\top X)^{-1}$$ is  $$(1/N_x + 1/N_y)$$

1. You can make a design matrix `X` for a two group comparison, either using `model.matrix` or simply with:

    
    ```r
    X <- cbind(rep(1,Nx + Ny),rep(c(0,1),c(Nx, Ny)))
    ```

    In order to compare two groups, where the first group has `Nx=5` samples and the second group has `Ny=7` samples, what is the element in the 1st row and 1st column of $$X^\top X$?




2. The other entries of $$X^\top X$$ are all the same. What is this number?



Now we just need to invert the matrix to obtain $$(X^\top X)^{-1}$$. The formula for matrix inversion for a 2x2 matrix is actually relatively easy to memorize:

$$ \,
\begin{pmatrix}
a&b\\
c&d
\end{pmatrix}^{-1} = \frac{1}{ad - bc}
\begin{pmatrix}
d&-b\\
-c&a
\end{pmatrix}
$$

The element of the inverse in the 2nd row and the 2nd column is the element which will be used to calculate the standard error of the second coefficient of the linear model. This is $$a / (ad - bc) $$. And for our two group comparison, we saw that $$a = N_x + N_y$$ and the $$b = c = d = N_y$$. So it follows that this element is:

$$
\frac{N_x + N_y}{(N_x + N_y) N_y - N_y N_y}
$$

which simplifies to:

$$
\frac{N_x + N_y}{N_x N_y} = 1/N_y + 1/N_x
$$
