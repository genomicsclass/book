---
layout: page
title: Analyzing DNA methylation data
---




In this unit we will show an example of analyzing methylation data. We will use colon cancer data from TCGA. The data was created with the Illumina 450K array and we have already processed the raw data to create matrix with methylation measurements. The script that creates these ojects is here: https://github.com/genomicsclass/labs/blob/master/Rscripts/read_tcga_meth.R

Let's begin by loading the data

```r
# devtools::install_github("genomicsclass/coloncancermeth")
library(S4Vectors)
library(coloncancermeth)
data(coloncancermeth)
```

We know have three tables one containing the methylation data, one with information about the samples or columns of the data matrix, and granges object with the genomic location of the CpGs represetned in the rows of the data matrix


```r
dim(meth) ##this is the methylation data
```

```
## [1] 485512     26
```

```r
dim(pd) ##this is sample information
```

```
## Loading required package: IRanges
```

```
## [1]  26 105
```

```r
length(gr)
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## [1] 485512
```

The `pd` object includes clinical information. One the columns tells us if the sample is from colon cancer or from normal tissue


```r
colnames(pd)
```

```
##   [1] "bcr_patient_barcode"                                              
##   [2] "bcr_sample_uuid"                                                  
##   [3] "bcr_sample_barcode"                                               
##   [4] "NCNNCT_OthMethONSP"                                               
##   [5] "current_weight"                                                   
##   [6] "days_to_collection"                                               
##   [7] "days_to_sample_procurement"                                       
##   [8] "freezing_method"                                                  
##   [9] "initial_weight"                                                   
##  [10] "intermediate_dimension"                                           
##  [11] "longest_dimension"                                                
##  [12] "lymphatic_invasion"                                               
##  [13] "margins_involved"                                                 
##  [14] "method_of_sample_procurement"                                     
##  [15] "number_regional_lymphnodes_exam"                                  
##  [16] "number_regional_lymphnodes_pos"                                   
##  [17] "oct_embedded"                                                     
##  [18] "pathology_report_uuid"                                            
##  [19] "primary_or_metastatic_status"                                     
##  [20] "sample_type"                                                      
##  [21] "sample_type_id"                                                   
##  [22] "shortest_dimension"                                               
##  [23] "time_between_clamping_and_freezing"                               
##  [24] "time_between_excision_and_freezing"                               
##  [25] "venous_invasion"                                                  
##  [26] "verification_by_bcr"                                              
##  [27] "vial_number.sample"                                               
##  [28] "bcr_patient_barcode.tumor"                                        
##  [29] "tumor_necrosis_percent"                                           
##  [30] "tumor_nuclei_percent"                                             
##  [31] "tumor_weight"                                                     
##  [32] "vial_number.tumor"                                                
##  [33] "bcr_patient_barcode.normal"                                       
##  [34] "days_to_normal_sample_procurement"                                
##  [35] "method_of_normal_sample_procurement"                              
##  [36] "normal_control_type"                                              
##  [37] "normal_tissue_anatomic_site"                                      
##  [38] "normal_tissue_proximity"                                          
##  [39] "vial_number"                                                      
##  [40] "ncedna_dna_conc"                                                  
##  [41] "ncedna_dna_qm"                                                    
##  [42] "ncedna_dna_qty"                                                   
##  [43] "ncedna_dna_vol"                                                   
##  [44] "patient.age_at_initial_pathologic_diagnosis"                      
##  [45] "patient.ajcc_cancer_staging_handbook_edition"                     
##  [46] "patient.anatomic_organ_subdivision"                               
##  [47] "patient.anatomic_site_colorectal"                                 
##  [48] "patient.bcr_patient_uuid"                                         
##  [49] "patient.braf_gene_analysis_performed"                             
##  [50] "patient.braf_gene_analysis_result"                                
##  [51] "patient.circumferential_resection_margin"                         
##  [52] "patient.colon_polyps_present"                                     
##  [53] "patient.date_of_form_completion"                                  
##  [54] "patient.date_of_initial_pathologic_diagnosis"                     
##  [55] "patient.days_to_birth"                                            
##  [56] "patient.days_to_death"                                            
##  [57] "patient.days_to_initial_pathologic_diagnosis"                     
##  [58] "patient.days_to_last_followup"                                    
##  [59] "patient.days_to_last_known_alive"                                 
##  [60] "patient.distant_metastasis_pathologic_spread"                     
##  [61] "patient.ethnicity"                                                
##  [62] "patient.gender"                                                   
##  [63] "patient.height"                                                   
##  [64] "patient.histological_type"                                        
##  [65] "patient.history_of_colon_polyps"                                  
##  [66] "patient.icd_10"                                                   
##  [67] "patient.icd_o_3_histology"                                        
##  [68] "patient.icd_o_3_site"                                             
##  [69] "patient.informed_consent_verified"                                
##  [70] "patient.kras_gene_analysis_performed"                             
##  [71] "patient.kras_mutation_codon"                                      
##  [72] "patient.kras_mutation_found"                                      
##  [73] "patient.loss_expression_of_mismatch_repair_proteins_by_ihc"       
##  [74] "patient.loss_expression_of_mismatch_repair_proteins_by_ihc_result"
##  [75] "patient.lymph_node_examined_count"                                
##  [76] "patient.lymphatic_invasion"                                       
##  [77] "patient.lymphnode_pathologic_spread"                              
##  [78] "patient.microsatellite_instability"                               
##  [79] "patient.non_nodal_tumor_deposits"                                 
##  [80] "patient.number_of_abnormal_loci"                                  
##  [81] "patient.number_of_first_degree_relatives_with_cancer_diagnosis"   
##  [82] "patient.number_of_loci_tested"                                    
##  [83] "patient.number_of_lymphnodes_positive_by_he"                      
##  [84] "patient.number_of_lymphnodes_positive_by_ihc"                     
##  [85] "patient.patient_id"                                               
##  [86] "patient.perineural_invasion_present"                              
##  [87] "patient.person_neoplasm_cancer_status"                            
##  [88] "patient.preoperative_pretreatment_cea_level"                      
##  [89] "patient.pretreatment_history"                                     
##  [90] "patient.primary_lymph_node_presentation_assessment"               
##  [91] "patient.primary_tumor_pathologic_spread"                          
##  [92] "patient.prior_diagnosis"                                          
##  [93] "patient.race"                                                     
##  [94] "patient.residual_tumor"                                           
##  [95] "patient.synchronous_colon_cancer_present"                         
##  [96] "patient.tissue_source_site"                                       
##  [97] "patient.tumor_stage"                                              
##  [98] "patient.tumor_tissue_site"                                        
##  [99] "patient.venous_invasion"                                          
## [100] "patient.vital_status"                                             
## [101] "patient.weight"                                                   
## [102] "Basename"                                                         
## [103] "Status"                                                           
## [104] "Tissue"                                                           
## [105] "Sex"
```

```r
table(pd$Status)
```

```
## 
## normal cancer 
##      9     17
```

```r
normalIndex <- which(pd$Status=="normal")
cancerlIndex <- which(pd$Status=="cancer")
```


Let's start by taking a quick look at the distribution of methylation measurements for the normal samples:


```r
i=normalIndex[1]
plot(density(meth[,i],from=0,to=1),main="",ylim=c(0,3),type="n")
for(i in normalIndex){
  lines(density(meth[,i],from=0,to=1),col=1)
}
### Add the cancer samples
for(i in cancerlIndex){
  lines(density(meth[,i],from=0,to=1),col=2)
}
```

![plot of chunk unnamed-chunk-4](figure/methylation-unnamed-chunk-4-1.png)

We are interested in finding regions of the genome that are different between cancer and normal samples. Furthermore, we want regions that are consistenly different therefore we can treat this as an inference problem. We can compute a t-statistic for each CpG:


```r
library(limma)
X<-model.matrix(~pd$Status)
fit<-lmFit(meth,X)
eb <- eBayes(fit)
```

A volcano plot reveals many differences:


```r
library(rafalib)
splot(fit$coef[,2],-log10(eb$p.value[,2]),xlab="Effect size",ylab="-log10 p-value")
```

![plot of chunk unnamed-chunk-6](figure/methylation-unnamed-chunk-6-1.png)

If we have reason to believe for DNA methylation to have an effect on gene expression a region of the genome needs to be affected, not just a single CpG, we should look beyond. Here is plot of the region surrounding the top hit:


```r
library(GenomicRanges)
i <- which.min(eb$p.value[,2])
middle <- gr[i,]
Index<-gr%over%(middle+10000)
cols=ifelse(pd$Status=="normal",1,2)
chr=as.factor(seqnames(gr))
pos=start(gr)

plot(pos[Index],fit$coef[Index,2],type="b",xlab="genomic location",ylab="difference")
```

![plot of chunk unnamed-chunk-7](figure/methylation-unnamed-chunk-7-1.png)

```r
matplot(pos[Index],meth[Index,],col=cols,xlab="genomic location")
```

![plot of chunk unnamed-chunk-7](figure/methylation-unnamed-chunk-7-2.png)

We can search for these regions explicitly instead of searching for single points, as explained by Jaffe and Irizarry (2012) [http://www.ncbi.nlm.nih.gov/pubmed/22422453]. 

If we are going to perform regional analysis we first have to define a region. But one issue is that not only do we have to separate the analysis by chromosome but that within each chromosome we usually have big gaps creating subgroups of regions to be analyzed.


```r
chr1Index <- which(chr=="chr1")
hist(log10(diff(pos[chr1Index])),main="",xlab="log 10 method")
```

![plot of chunk unnamed-chunk-8](figure/methylation-unnamed-chunk-8-1.png)

We can create groups in the following way.


```r
# BiocManager::install("bumphunter")
library(bumphunter)
cl=clusterMaker(chr,pos,maxGap=500)
table(table(cl)) ##shows the number of regions with 1,2,3, ... points in them
```

```
## 
##      1      2      3      4      5      6      7      8      9     10 
## 141457  18071  13227   6473   5144   3748   2517   2135   2029   1878 
##     11     12     13     14     15     16     17     18     19     20 
##   1792   1570   1269    933    684    472    337    240    181     99 
##     21     22     23     24     25     26     27     28     29     30 
##    113     62     57     42     36     39     26     17     21     19 
##     31     32     33     34     35     36     37     38     39     40 
##     12     12     12      7      7     12      3      9      4      5 
##     41     42     43     44     45     46     48     49     50     51 
##      7      7      9      7      7      8      3      5      5      2 
##     52     53     54     55     56     57     58     59     60     61 
##      1      5      2      4      1      1      3      2      1      4 
##     62     63     64     65     67     68     70     71     73     74 
##      3      4      1      2      2      1      2      1      2      2 
##     76     78     80     82     83     85     87     88     89     90 
##      2      1      3      2      1      2      1      2      1      1 
##     91     92     93    112    117    137    141    181 
##      1      1      1      2      1      1      1      1
```


Now let's consider two example regions:


```r
###Select the region with the smallest value
Index<- which(cl==cl[which.min(fit$coef[,2])])
matplot(pos[Index],meth[Index,],col=cols,pch=1,xlab="genomic location",ylab="methylation")
```

![plot of chunk unnamed-chunk-10](figure/methylation-unnamed-chunk-10-1.png)

```r
x1=pos[Index]
y1=fit$coef[Index,2]
plot(x1,y1,xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=0,lty=2)
abline(h=c(-.1,.1),lty=2)
```

![plot of chunk unnamed-chunk-10](figure/methylation-unnamed-chunk-10-2.png)

This region shows only a single CpG as different. In contrast, notice this region:


```r
Index=which(cl==72201) ##we know this is a good example from analysis we have already performed

matplot(pos[Index],meth[Index,],col=cols,pch=1,xlab="genomic location",ylab="methylation")
```

![plot of chunk unnamed-chunk-11](figure/methylation-unnamed-chunk-11-1.png)

```r
x2=pos[Index]
y2=fit$coef[Index,2]
plot(x2,y2,xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=0,lty=2)
abline(h=c(-.1,.1),lty=2)
```

![plot of chunk unnamed-chunk-11](figure/methylation-unnamed-chunk-11-2.png)

<a name="DMR"></a>

If we are interested in prioritizing regions over single points, we need an alternative approach. If we assume that the real signal is smooth, we could use statistical smoothing techniques such as loess. Here is an example two regions above


```r
lfit <- loess(y1~x1,degree=1,family="symmetric",span=1/2)
plot(x1,y1,xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=c(-.1,0,.1),lty=2)
lines(x1,lfit$fitted,col=2)
```

![plot of chunk unnamed-chunk-12](figure/methylation-unnamed-chunk-12-1.png)

```r
lfit <- loess(y2~x2,degree=1,family="symmetric",span=1/2)
plot(x2,y2,xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=c(-.1,0,.1),lty=2)
lines(x2,lfit$fitted,col=2)
```

![plot of chunk unnamed-chunk-12](figure/methylation-unnamed-chunk-12-2.png)


The bumphunter automates this procedure:


```r
res<-bumphunter(meth,X,chr=chr,pos=pos,cluster=cl,cutoff=0.1,B=0)
```

```
## [bumphunterEngine] Using a single core (backend: doSEQ, version: 1.4.7).
```

```
## [bumphunterEngine] Computing coefficients.
```

```
## [bumphunterEngine] Finding regions.
```

```
## [bumphunterEngine] Found 68682 bumps.
```

```r
tab<-res$table
```

We now have a list of regions instead of single points. Here we look at the region with the highest rank if we order by area:


```r
Index=(tab[1,7]-3):(tab[1,8]+3)
matplot(pos[Index],meth[Index,,drop=TRUE],col=cols,pch=1,xlab="genomic location",ylab="Methylation",ylim=c(0,1))
```

![plot of chunk unnamed-chunk-14](figure/methylation-unnamed-chunk-14-1.png)

```r
plot(pos[Index],res$fitted[Index,1],xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=c(-0.1,0,.1),lty=2)
```

![plot of chunk unnamed-chunk-14](figure/methylation-unnamed-chunk-14-2.png)

The function also allows from smoothing and permutation based inference for the regions. However, we do not recommend running the function with these options without the ability to parallelize. 

