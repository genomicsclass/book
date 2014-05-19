---
layout: page
title: Footnotes for Data Analysis for Genomics
---

# Week 1

## Transformation and robust summaries
### Robust statistics

Robust Statistics, Peter. J. Huber and Elvezio M. Ronchetti, Wiley, 2009.

Introduction to Robust Estimation and Hypothesis Testing, Rand R. Wilcox, 2012.

Robust Statistics: The Approach Based on Influence Functions, Frank R. Hampel, Elvezio M. Ronchetti, Peter J. Rousseeuw, Werner A. Stahel


----

# Week 2

## Basic Bioconductor infrastructure
For more information about the `GenomicRanges` package, check out the PLOS Comp Bio paper, which the authors of GenomicRanges published:

<http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>

Also the software vignettes have a lot of details about the functionality. Check out "An Introduction to Genomic Ranges Classes":

<http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf>

All of the vignette PDFs are available here:

<http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html>


----

# Week 3

## Basic inference for microarray data
Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>


----

## Statistical inference
C. Kendziorski, R. A. Irizarry, K.-S. Chen, J. D. Haag, and M. N. Gould, "On the utility of pooling biological samples in microarray experiments", PNAS, 2005. <http://www.pnas.org/content/102/12/4252.long>

Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>


----

## Rank tests
C. Kendziorski, R. A. Irizarry, K.-S. Chen, J. D. Haag, and M. N. Gould, "On the utility of pooling biological samples in microarray experiments", PNAS, 2005. <http://www.pnas.org/content/102/12/4252.long>


----

# Week 4

## Modeling
Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>


----

## Normalization
### loess

W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of Statistical Models in S eds J.M. Chambers and T.J. Hastie, Wadsworth & Brooks/Cole.

### Quantile normalization

Bolstad BM, Irizarry RA, Astrand M, Speed TP. "A comparison of normalization methods for high density oligonucleotide array data based on variance and bias." Bioinformatics. 2003. <http://www.ncbi.nlm.nih.gov/pubmed/12538238>

### Variance stabilization

For microarray:

Wolfgang Huber, Anja von Heydebreck, Holger Sültmann, Annemarie Poustka and Martin Vingron, "Variance stabilization applied to microarray data calibration and to the quantification of differential expression" Bioinformatics, 2002. <http://bioinformatics.oxfordjournals.org/content/18/suppl_1/S96.short>

B.P. Durbin, J.S. Hardin, D.M. Hawkins and D.M. Rocke, "A variance-stabilizing transformation for gene-expression microarray data", Bioinformatics. 2002. <http://bioinformatics.oxfordjournals.org/content/18/suppl_1/S105>

Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka, Martin Vingron, "Parameter estimation for the calibration and variance stabilization of microarray data" Stat Appl Mol Biol Genet, 2003. <http://dx.doi.org/10.2202/1544-6115.1008>

For NGS read counts:

Simon Anders and Wolfgang Huber, "Differential expression analysis for sequence count data", Genome Biology, 2010. <http://genomebiology.com/2010/11/10/r106>

General discussion of variance stabilization:

Robert Tibshirani, "Variance Stabilization and the Bootstrap" Biometrika, 1988. <http://www.jstor.org/discover/10.2307/2336593>


----

# Week 5

## Distance lecture
### MDS references

Wikipedia: <http://en.wikipedia.org/wiki/Multidimensional_scaling>

- Cox, T. F. and Cox, M. A. A. (2001) Multidimensional Scaling. Second edition. Chapman and Hall.

- Gower, J. C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53, 325–328.

- Torgerson, W. S. (1958). Theory and Methods of Scaling. New York: Wiley.

### Hierarchical clustering

Wikipedia: <http://en.wikipedia.org/wiki/Hierarchical_clustering>

a subset of the most recent references from `?hclust`:

- Legendre, P. and L. Legendre (2012). Numerical Ecology, 3rd English ed. Amsterdam: Elsevier Science BV.

- Murtagh, F. and Legendre, P. (2013). Ward's hierarchical agglomerative clustering method: which algorithms implement Ward's criterion? Journal of Classification (in press).


----

# Week 6

## Batch adjustment
### Principal Components Analysis (PCA)

Jolliffe, Ian. Principal component analysis. John Wiley & Sons, Ltd, 2005.

Dunteman, George H. Principal components analysis. No. 69. Sage, 1989.


----

# Week 7

## Gene set testing
Masuno K, Haldar SM, Jeyaraj D, Mailloux CM, Huang X, Panettieri RA Jr, Jain MK, Gerber AN., "Expression profiling identifies Klf15 as a glucocorticoid target that regulates airway hyperresponsiveness". Am J Respir Cell Mol Biol. 2011.

<http://www.ncbi.nlm.nih.gov/pubmed/21257922>

The ROAST method

Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK. "ROAST: rotation gene set tests for complex microarray experiments". Bioinformatics. 2010.

<http://www.ncbi.nlm.nih.gov/pubmed/20610611>


----

## Multiple testing
Yoav Benjamini and Yosef Hochberg, "Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing". Journal of the Royal Statistical Society. 1995.

<http://www.jstor.org/discover/10.2307/2346101>


----

## Using limma for microarray analysis
Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>


----

# Week 8

