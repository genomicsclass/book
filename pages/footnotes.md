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
### Methods within the limma package

Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK. "ROAST: rotation gene set tests for complex microarray experiments". Bioinformatics. 2010.

<http://www.ncbi.nlm.nih.gov/pubmed/20610611>

Di Wu and Gordon K. Smyth, "Camera: a competitive gene set test accounting for inter-gene correlation" Nucleic Acids Research, 2012.

<http://nar.oxfordjournals.org/content/40/17/e133>

### GSEA

Subramanian A1, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP, "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles" Proc Natl Acad Sci U S A. 2005.

<http://www.ncbi.nlm.nih.gov/pubmed/16199517>

### Correlation within gene sets

William T. Barry, Andrew B. Nobel, and Fred A. Wright, "A statistical framework for testing functional categories in microarray data" Ann. Appl. Stat, 2008.

<http://projecteuclid.org/euclid.aoas/1206367822>

William Barry has a package `safe` in Bioconductor for gene set testing with resampling.

<http://www.bioconductor.org/packages/release/bioc/html/safe.html>

Daniel M Gatti, William T Barry, Andrew B Nobel, Ivan Rusyn and Fred A Wright, "Heading Down the Wrong Pathway: on the Influence of Correlation within Gene Sets", BMC Genomics, 2010.

<http://www.biomedcentral.com/1471-2164/11/574#B24>

### Gene sets and power

The following article points out an issue with gene set testing: the power to detect differential expression for an individual gene depends on the number of NGS reads which align to that gene, which depends on the transcript length among other factors.

Alicia Oshlack* and Matthew J Wakefield, "Transcript length bias in RNA-seq data confounds systems biology", Biology Direct, 2009.

<http://www.biologydirect.com/content/4/1/14>

### The dataset used in this lab

Masuno K, Haldar SM, Jeyaraj D, Mailloux CM, Huang X, Panettieri RA Jr, Jain MK, Gerber AN., "Expression profiling identifies Klf15 as a glucocorticoid target that regulates airway hyperresponsiveness". Am J Respir Cell Mol Biol. 2011.

<http://www.ncbi.nlm.nih.gov/pubmed/21257922>


----

## Multiple testing
Yoav Benjamini and Yosef Hochberg, "Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing". Journal of the Royal Statistical Society. 1995.

<http://www.jstor.org/discover/10.2307/2346101>


----

## Using limma for microarray analysis
Smyth GK, "Linear models and empirical bayes methods for assessing differential expression in microarray experiments". Stat Appl Genet Mol Biol. 2004 <http://www.ncbi.nlm.nih.gov/pubmed/16646809>


----

# Week 8

## ChIP-seq analysis
### Model-based Analysis for ChIP-Seq (MACS)

Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. "Model-based Analysis of ChIP-Seq (MACS)". Genome Biol. 2008.

<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2592715/>

Software: 

<http://liulab.dfci.harvard.edu/MACS/>

### Motif finding

Wikipedia's article on DNA sequence motifs: <http://en.wikipedia.org/wiki/Sequence_motif>

A non-comprehensive list of software for motif finding:

- [MEME/DREME](http://meme.nbcr.net/meme/)

- [RSAT peak-motifs](http://rsat.ulb.ac.be/peak-motifs_form.cgi)

- [motifRG (Bioconductor)](http://www.bioconductor.org/packages/release/bioc/html/motifRG.html)

- [rGADEM (Bioconductor)](http://www.bioconductor.org/packages/release/bioc/html/rGADEM.html)

A survey of motif finding algorithms: <http://www.biomedcentral.com/1471-2105/8/S7/S21/>


----

## NGS read counting
### Methods for counting reads which overlap features

Bioconductor packages:

- `summarizeOverlaps` in the `GenomicAlignments` package

<http://www.bioconductor.org/packages/devel/bioc/html/GenomicAlignments.html>

- `featureCounts` in the `Rsubread` package

Liao Y, Smyth GK, Shi W., "featureCounts: an efficient general purpose program for assigning sequence reads to genomic features." Bioinformatics. 2014

<http://www.ncbi.nlm.nih.gov/pubmed/24227677>

<http://bioinf.wehi.edu.au/featureCounts/>

- `easyRNAseq` package

Delhomme N1, Padioleau I, Furlong EE, Steinmetz LM. "easyRNASeq: a bioconductor package for processing RNA-Seq data." Bioinformatics. 2012.

<http://www.ncbi.nlm.nih.gov/pubmed/22847932>

<http://www.bioconductor.org/packages/release/bioc/html/easyRNASeq.html>

Command line tools: 

- `htseq-count`, a program in the `htseq` Python package

Simon Anders, Paul Theodor Pyl, Wolfgang Huber.

HTSeq — A Python framework to work with high-throughput sequencing data

bioRxiv preprint (2014), doi: [10.1101/002824](http://dx.doi.org/10.1101/002824)

<http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>

- `bedtools` <https://code.google.com/p/bedtools/>

- `bedops` <https://code.google.com/p/bedops/>


----

## RNA-seq analysis
### Introduction

Mortazavi A, Williams BA, McCue K, Schaeffer L, Wold B., "Mapping and quantifying mammalian transcriptomes by RNA-Seq", Nat Methods. 2008.

<http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1226.html>

John C. Marioni, Christopher E. Mason, Shrikant M. Mane, Matthew Stephens, and Yoav Gilad, "RNA-seq: An assessment of technical reproducibility and comparison with gene expression arrays" Genome Res. 2008.

<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2527709/>

Trapnell C, Williams BA, Pertea G, Mortazavi AM, Kwan G, van Baren MJ, Salzberg SL, Wold B, Pachter L.,  "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation", Nature Biotechnology, 2010.

<http://www.nature.com/nbt/journal/v28/n5/full/nbt.1621.html>

<http://cufflinks.cbcb.umd.edu/>

### Hammer et al

Hammer P, Banck MS, Amberg R, Wang C, Petznick G, Luo S, Khrebtukova I, Schroth GP, Beyerlein P, Beutler AS. "mRNA-seq with agnostic splice site discovery for nervous system transcriptomics tested in chronic pain." Genome Res. 2010

<http://www.ncbi.nlm.nih.gov/pubmed?term=20452967>

### ReCount

Frazee AC, Langmead B, Leek JT. "ReCount: a multi-experiment resource of analysis-ready RNA-seq gene count datasets". BMC Bioinformatics 12:449 <http://www.ncbi.nlm.nih.gov/pubmed/22087737>

### Negative Binomial methods for differential expression of count data

All the following methods are available on Bioconductor:

- `edgeR`

Mark D. Robinson, Davis J. McCarthy, and Gordon K. Smyth, "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data" Bioinformatics 2010.

<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/>

- `DESeq` (the latest version is a separate package, `DESeq2`)

Simon Anders and Wolfgang Huber, "Differential expression analysis for sequence count data", Genome Biology 2010.

<http://genomebiology.com/2010/11/10/r106>

- `DSS`

Hao Wu, Chi Wang, Zhijin Wu, "A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data" Biostatistics 2013.

<http://biostatistics.oxfordjournals.org/content/14/2/232>

### Transformation followed by linear model methods

`voom` in the `limma` Bioconductor package

Charity W Law, Yunshun Chen, Wei Shi and Gordon K Smyth, "voom: precision weights unlock linear model analysis tools for RNA-seq read counts", Genome Biology. 2014.

<http://genomebiology.com/2014/15/2/R29>

### Resampling-based methods

`SAMseq` in the `samr` package on CRAN

Jun Li and Robert Tibshirani, "Finding consistent patterns: A nonparametric approach for identifying differential expression in RNA-Seq data", Stat Methods Med Res. 2013.

<http://smm.sagepub.com/content/22/5/519.short>

### Incorporating isoform-abundance

- `Cuffdiff` (the latest version is `Cuffdiff2`)

Trapnell C, Hendrickson DG, Sauvageau M, Goff L, Rinn JL, Pachter L., "Differential analysis of gene regulation at transcript resolution with RNA-seq" Nat Biotechnol. 2013.

<http://www.ncbi.nlm.nih.gov/pubmed/23222703>

- `BitSeq`

Peter Glaus, Antti Honkela, and Magnus Rattray, "Identifying differentially expressed transcripts from RNA-seq data with biological variation", Bioinformatics. 2012.

<http://bioinformatics.oxfordjournals.org/content/28/13/1721>


----

