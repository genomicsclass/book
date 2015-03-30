---
title: "Annotating genes and other structural and functional genomic elements"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---






# Programming with gene catalogs

The simplest way to get a detailed collection of gene models
with Bioconductor, for human, rat, or mouse,
is to use packages implementing the OrganismDb concept.
We'll illustrate with mouse.


```r
library(Mus.musculus)
Mus.musculus
```

```
## class: OrganismDb 
## Annotation resources:
## [1] "GO.db"                             
## [2] "org.Mm.eg.db"                      
## [3] "TxDb.Mmusculus.UCSC.mm10.knownGene"
## Annotation relationships:
##      xDbs           yDbs                                 xKeys     
## [1,] "GO.db"        "org.Mm.eg.db"                       "GOID"    
## [2,] "org.Mm.eg.db" "TxDb.Mmusculus.UCSC.mm10.knownGene" "ENTREZID"
##      yKeys   
## [1,] "GO"    
## [2,] "GENEID"
## For more details, please see the show methods for the component objects listed above.
```

The NCBI Entrez gene identifiers are keys in this database
interface.


```r
mk = keys(Mus.musculus, keytype="ENTREZID")
mk[1:5]
```

```
## [1] "11287" "11298" "11302" "11303" "11304"
```

The keys can be used to query for values of associated attributes.
The available attributes can be listed directly with "columns":


```r
columns(Mus.musculus)
```

```
##  [1] "GOID"         "TERM"         "ONTOLOGY"     "DEFINITION"  
##  [5] "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
##  [9] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"      
## [13] "CHRLOCEND"    "ENZYME"       "PATH"         "PMID"        
## [17] "REFSEQ"       "SYMBOL"       "UNIGENE"      "ENSEMBL"     
## [21] "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"     "UNIPROT"     
## [25] "GO"           "EVIDENCE"     "GOALL"        "EVIDENCEALL" 
## [29] "ONTOLOGYALL"  "MGI"          "CDSID"        "CDSNAME"     
## [33] "CDSCHROM"     "CDSSTRAND"    "CDSSTART"     "CDSEND"      
## [37] "EXONID"       "EXONNAME"     "EXONCHROM"    "EXONSTRAND"  
## [41] "EXONSTART"    "EXONEND"      "GENEID"       "TXID"        
## [45] "EXONRANK"     "TXNAME"       "TXCHROM"      "TXSTRAND"    
## [49] "TXSTART"      "TXEND"
```


```r
select(Mus.musculus, keys=mk[1:5],
  keytype="ENTREZID", columns=c("SYMBOL", "CHRLOC", "GENENAME"))
```

```
##   ENTREZID     CHRLOC CHRLOCCHR SYMBOL
## 1    11287 -128483567         6    Pzp
## 2    11298  116593687        11  Aanat
## 3    11302 -120007316        11   Aatk
## 4    11303  -53030789         4  Abca1
## 5    11304  122044460         3  Abca4
##                                              GENENAME
## 1                              pregnancy zone protein
## 2                  arylalkylamine N-acetyltransferase
## 3                apoptosis-associated tyrosine kinase
## 4 ATP-binding cassette, sub-family A (ABC1), member 1
## 5 ATP-binding cassette, sub-family A (ABC1), member 4
```

We can move directly to GRanges representations of
addresses using TxDb.  Here we generate a list indexed
by Entrez identifier.


```r
mt = transcriptsBy(Mus.musculus, by="gene")
mt
```

```
## GRangesList object of length 23725:
## $$100009600 
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames               ranges strand |     tx_id     tx_name
##          <Rle>            <IRanges>  <Rle> | <integer> <character>
##   [1]     chr9 [21062393, 21067925]      - |     31744  uc009veu.1
##   [2]     chr9 [21062393, 21075496]      - |     31745  uc033jjg.1
## 
## $$100009609 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand | tx_id    tx_name
##   [1]     chr7 [84940169, 84964009]      - | 26069 uc012fog.1
## 
## $$100009614 
## GRanges object with 1 range and 2 metadata columns:
##       seqnames               ranges strand | tx_id    tx_name
##   [1]    chr10 [77711446, 77712009]      + | 34226 uc011xhj.2
## 
## ...
## <23722 more elements>
## -------
## seqinfo: 66 sequences (1 circular) from mm10 genome
```

Visualization of gene models can occur easily
using the custom package.


```r
library(ph525x)
modPlot("Pzp", genome="mm10", annoResource=Mus.musculus, 
   collapse=FALSE, useGeneSym=FALSE)
```

```
## Loading required package: Gviz
```

![plot of chunk lkmodmm](figure/anno5genes-lkmodmm-1.png) 

Gene function and localization information can be retrieved using
various types of key.

```r
select(Mus.musculus, keys="Pzp",
  keytype="SYMBOL", columns=c("GO", "TERM"))
```

```
## Warning in .generateExtraRows(tab, keys, jointype): 'select' resulted in
## 1:many mapping between keys and return rows
```

```
## Warning in .generateExtraRows(tab, keys, jointype): 'select' resulted in
## 1:many mapping between keys and return rows
```

```
##   SYMBOL         GO EVIDENCE ONTOLOGY
## 1    Pzp GO:0004866      IDA       MF
## 2    Pzp GO:0004867      IEA       MF
## 3    Pzp GO:0005576      TAS       CC
## 4    Pzp GO:0005615      IEA       CC
## 5    Pzp GO:0007566      IGI       BP
## 6    Pzp GO:0010466      IEA       BP
## 7    Pzp GO:0010951      IDA       BP
## 8    Pzp GO:0030414      IEA       MF
## 9    Pzp GO:0032403      ISO       MF
##                                            TERM
## 1              endopeptidase inhibitor activity
## 2  serine-type endopeptidase inhibitor activity
## 3                          extracellular region
## 4                           extracellular space
## 5                           embryo implantation
## 6     negative regulation of peptidase activity
## 7 negative regulation of endopeptidase activity
## 8                  peptidase inhibitor activity
## 9                       protein complex binding
```

## BioMart

A vast collection of biological annotation can
be obtained from the Biomart servers.

A hierarchical interface is used.  We begin by
selecting a "mart" to use, and then a dataset
within the mart.  The mart instance is updated.


```r
library(biomaRt)
head(listMarts())
```

```
##               biomart                             version
## 1             ensembl        ENSEMBL GENES 78 (SANGER UK)
## 2                 snp    ENSEMBL VARIATION 78 (SANGER UK)
## 3 functional_genomics   ENSEMBL REGULATION 78 (SANGER UK)
## 4                vega                VEGA 58  (SANGER UK)
## 5       fungi_mart_25           ENSEMBL FUNGI 25 (EBI UK)
## 6 fungi_variations_25 ENSEMBL FUNGI VARIATION 25 (EBI UK)
```

```r
# m = useMart("ensembl")  # typical, but if biomart is down, use:
m=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
head(listDatasets(mart=m))
```

```
##                          dataset
## 1         oanatinus_gene_ensembl
## 2        cporcellus_gene_ensembl
## 3        gaculeatus_gene_ensembl
## 4         lafricana_gene_ensembl
## 5 itridecemlineatus_gene_ensembl
## 6        choffmanni_gene_ensembl
##                                  description version
## 1     Ornithorhynchus anatinus genes (OANA5)   OANA5
## 2            Cavia porcellus genes (cavPor3) cavPor3
## 3     Gasterosteus aculeatus genes (BROADS1) BROADS1
## 4         Loxodonta africana genes (loxAfr3) loxAfr3
## 5 Ictidomys tridecemlineatus genes (spetri2) spetri2
## 6        Choloepus hoffmanni genes (choHof1) choHof1
```

```r
m = useDataset("hsapiens_gene_ensembl", mart=m)
m
```

```
## Object of class 'Mart':
##  Using the ENSEMBL_MART_ENSEMBL BioMart database
##  Using the hsapiens_gene_ensembl dataset
```

We obtain data by issuing requests through filters.
Requests are framed using attributes.


```r
head(listFilters(m))
```

```
##              name     description
## 1 chromosome_name Chromosome name
## 2           start Gene Start (bp)
## 3             end   Gene End (bp)
## 4      band_start      Band Start
## 5        band_end        Band End
## 6    marker_start    Marker Start
```

```r
head(listAttributes(m))
```

```
##                    name           description
## 1       ensembl_gene_id       Ensembl Gene ID
## 2 ensembl_transcript_id Ensembl Transcript ID
## 3    ensembl_peptide_id    Ensembl Protein ID
## 4       ensembl_exon_id       Ensembl Exon ID
## 5           description           Description
## 6       chromosome_name       Chromosome Name
```

The query interface is the getBM function.  Here's an
example of getting three types of identifiers for ORMDL3
and its transcripts.


```r
getBM(attributes=c("ensembl_gene_id", "entrezgene", "ucsc"), 
    filters="hgnc_symbol", values="ORMDL3", mart=m)
```

```
##   ensembl_gene_id entrezgene       ucsc
## 1 ENSG00000172057      94103 uc002htj.2
## 2 ENSG00000172057      94103 uc002htk.2
## 3 ENSG00000172057      94103
```

# AnnotationHub

A recent and evolving resource for various annotation resources
is the AnnotationHub package.  The idea is to give convenient
access for programming to items like the UCSC genome browser
track set, or the datasets organized on [the epigenomic roadmap](http://www.roadmapepigenomics.org/).

To give a sense of the sort of data we need a convenient interface to,
consider this slice of information that forms a part of the
epigenomic road map.


```r
library(ph525x)
sydhTop()
```
Rows are transcription factors, columns are cell lines from
different organs and donors.
Our aim is to support
statistical analysis of binding patterns derived from ChIP-seq
experiments.

## Tracks from the epigenomic road map project

We begin by connecting to the hub.  Metadata are returned along
with the connection.


```r
library(AnnotationHub)
ah = AnnotationHub()
ah
```

```
## class: AnnotationHub 
## length: 10780 
## filters: none 
## hubUrl: http://annotationhub.bioconductor.org/ah 
## snapshotVersion: 3.0/1.6.0; snapshotDate: 2014-05-15
## hubCache: /Users/stvjc/.AnnotationHub
```

```r
head(names(ah))
```

```
## [1] "ensembl.release.69.fasta.ailuropoda_melanoleuca.ncrna.Ailuropoda_melanoleuca.ailMel1.69.ncrna.fa.rz"
## [2] "ensembl.release.69.fasta.ailuropoda_melanoleuca.pep.Ailuropoda_melanoleuca.ailMel1.69.pep.all.fa.rz"
## [3] "ensembl.release.69.fasta.anolis_carolinensis.cdna.Anolis_carolinensis.AnoCar2.0.69.cdna.all.fa.rz"  
## [4] "ensembl.release.69.fasta.anolis_carolinensis.ncrna.Anolis_carolinensis.AnoCar2.0.69.ncrna.fa.rz"    
## [5] "ensembl.release.69.fasta.anolis_carolinensis.pep.Anolis_carolinensis.AnoCar2.0.69.pep.all.fa.rz"    
## [6] "ensembl.release.69.fasta.bos_taurus.cdna.Bos_taurus.UMD3.1.69.cdna.all.fa.rz"
```

A flexible query interface is provided.  We know that "Sydh" is
used in the names of tracks related to the road map matrix shown
above.


```r
sydq = query(ah, "Sydh")
length(sydq)
```

```
## [1] 389
```

```r
head(names(sydq))
```

```
## [1] "goldenpath.hg19.encodeDCC.wgEncodeSydhTfbs.wgEncodeSydhTfbsHepg2Brca1a300IggrabPk.narrowPeak_0.0.1.RData"   
## [2] "goldenpath.hg19.encodeDCC.wgEncodeSydhHistone.wgEncodeSydhHistoneHct116H3k04me1UcdPk.narrowPeak_0.0.1.RData"
## [3] "goldenpath.hg19.encodeDCC.wgEncodeSydhHistone.wgEncodeSydhHistoneHct116H3k27acUcdPk.narrowPeak_0.0.1.RData" 
## [4] "goldenpath.hg19.encodeDCC.wgEncodeSydhHistone.wgEncodeSydhHistoneK562H3k27me3bUcdPk.narrowPeak_0.0.1.RData" 
## [5] "goldenpath.hg19.encodeDCC.wgEncodeSydhHistone.wgEncodeSydhHistoneK562H3k4me1UcdPk.narrowPeak_0.0.1.RData"   
## [6] "goldenpath.hg19.encodeDCC.wgEncodeSydhHistone.wgEncodeSydhHistoneK562H3k4me3bUcdPk.narrowPeak_0.0.1.RData"
```

A GRanges with information on binding of BRCA1 A300 (synthetic
peptide corresponding to a component of the BRCA1 gene) to 
DNA extracted from the HepG2 cell lines (liver-derived)
is:


```r
bh = ah[[ names(sydq)[1] ]]
bh
```

```
## GRanges object with 18329 ranges and 6 metadata columns:
##           seqnames                 ranges strand   |        name     score
##              <Rle>              <IRanges>  <Rle>   | <character> <integer>
##       [1]    chr20 [ 34329450,  34330798]      *   |           .      1000
##       [2]    chr20 [ 26188667,  26190662]      *   |           .      1000
##       [3]    chr20 [ 30945472,  30946282]      *   |           .      1000
##       [4]     chr2 [ 55844538,  55845276]      *   |           .      1000
##       [5]     chr2 [120516857, 120517730]      *   |           .      1000
##       ...      ...                    ...    ... ...         ...       ...
##   [18325]     chr6 [155281487, 155281866]      *   |           .       370
##   [18326]     chr5 [ 43105617,  43106142]      *   |           .       370
##   [18327]     chrY [ 13478010,  13478443]      *   |           .       370
##   [18328]    chr16 [ 12305337,  12306335]      *   |           .       370
##   [18329]     chr2 [ 43086126,  43086443]      *   |           .       370
##           signalValue    pValue    qValue      peak
##             <numeric> <numeric> <numeric> <integer>
##       [1]     7.96979  227.3722  222.8254       811
##       [2]     4.07638  192.4713  188.2255      1407
##       [3]    16.18614  163.8126  159.7429       395
##       [4]    11.94639  157.4212  153.4765       373
##       [5]    18.42133  145.2403  141.3924       427
##       ...         ...       ...       ...       ...
##   [18325]     2.13475   2.30744   2.02370        47
##   [18326]     2.07301   2.30744   2.02373       373
##   [18327]     1.57489   2.29858   2.01489       288
##   [18328]     1.45385   2.29533   2.01166       301
##   [18329]     1.99611   2.28660   2.00296       143
##   -------
##   seqinfo: 24 sequences from hg19 genome
```

## Resources from dbSNP

Here we will see
what is in the hub from dbSNP.


```r
qd = query(ah, "dbSNP")
head(names(qd))
```

```
## [1] "dbSNP.organisms.human_9606.VCF.clinvar_20130226.RData"            
## [2] "dbSNP.organisms.human_9606.VCF.common_and_clinical_20130226.RData"
## [3] "dbSNP.organisms.human_9606.VCF.ByChromosome.01.12156.ASW.RData"   
## [4] "dbSNP.organisms.human_9606.VCF.ByChromosome.01.12157.CHB.RData"   
## [5] "dbSNP.organisms.human_9606.VCF.ByChromosome.01.12158.CHD.RData"   
## [6] "dbSNP.organisms.human_9606.VCF.ByChromosome.01.12159.GIH.RData"
```

This indicates that representations of variants for HapMap populations
are available.  The data are in the [Variant Call Format](http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41) (VCF).
Let's retrieve one for population CEU, chr20.
We found the specific name by searching the full set of names(qd).


```r
aaa = ah[["dbSNP.organisms.human_9606.VCF.ByChromosome.20.1409.CEU.RData"]]
```

```
## 
## Attaching package: 'VariantAnnotation'
## 
## The following object is masked from 'package:base':
## 
##     tabulate
```

```r
dim(geno(aaa)$GT)
```

```
## [1] 68263   116
```

We'll learn more about Bioconductor's representation of VCF shortly.

## The AceView gene models


```r
allace = query(ah, "AceView")
names(allace)
```

```
## [1] "goldenpath.hg19.database.acembly_0.0.1.RData"
## [2] "goldenpath.hg18.database.acembly_0.0.1.RData"
## [3] "goldenpath.hg17.database.acembly_0.0.1.RData"
## [4] "goldenpath.mm9.database.acembly_0.0.1.RData"
```

```r
ace19 = ah[["goldenpath.hg19.database.acembly_0.0.1.RData"]]
ace19
```

```
## UCSC track 'acembly'
## UCSCData object with 259440 ranges and 5 metadata columns:
##            seqnames               ranges strand   |
##               <Rle>            <IRanges>  <Rle>   |
##        [1]     chr1 [66993147, 67171749]      +   |
##        [2]     chr1 [66999061, 67210057]      +   |
##        [3]     chr1 [66999253, 67109289]      +   |
##        [4]     chr1 [66999275, 67210767]      +   |
##        [5]     chr1 [66999276, 67210767]      +   |
##        ...      ...                  ...    ... ...
##   [259436]    chr22 [51222318, 51227614]      +   |
##   [259437]    chr22 [51222438, 51227612]      +   |
##   [259438]    chr22 [51223070, 51227614]      +   |
##   [259439]    chr22 [51223580, 51227607]      +   |
##   [259440]    chr22 [51237137, 51237759]      -   |
##                               name     score     itemRgb
##                        <character> <numeric> <character>
##        [1]            SGIP1.hAug10         0        <NA>
##        [2]            SGIP1.aAug10         0        <NA>
##        [3]            SGIP1.qAug10         0        <NA>
##        [4]            SGIP1.dAug10         0        <NA>
##        [5]            SGIP1.jAug10         0        <NA>
##        ...                     ...       ...         ...
##   [259436]       RPL23AP82.vaAug10         0        <NA>
##   [259437]        RPL23AP82.bAug10         0        <NA>
##   [259438]       RPL23AP82.vcAug10         0        <NA>
##   [259439]        RPL23AP82.cAug10         0        <NA>
##   [259440] setara.aAug10-unspliced         0        <NA>
##                           thick
##                       <IRanges>
##        [1] [67000042, 67171426]
##        [2] [67000042, 67208778]
##        [3] [67000042, 67109288]
##        [4] [67000042, 67208778]
##        [5] [67000042, 67160207]
##        ...                  ...
##   [259436] [51227615, 51227614]
##   [259437] [51222440, 51227409]
##   [259438] [51227615, 51227614]
##   [259439] [51223582, 51227409]
##   [259440] [51237760, 51237759]
##                                                      blocks
##                                               <IRangesList>
##        [1] [    1,   138] [ 6783,  6905] [98384, 98447] ...
##        [2] [    1,    30] [  869,   991] [92470, 92533] ...
##        [3] [    1,   103] [  677,   799] [92278, 92341] ...
##        [4] [    1,    81] [  655,   777] [92256, 92319] ...
##        [5] [    1,    80] [  654,   776] [92255, 92318] ...
##        ...                                              ...
##   [259436]           [   1,  132] [1284, 1404] [5003, 5297]
##   [259437]       [   1,   63] [1164, 1284] [4741, 4789] ...
##   [259438]           [   1,   45] [ 532,  652] [4251, 4545]
##   [259439]           [   1,  142] [3599, 3647] [3741, 4028]
##   [259440]                                         [1, 623]
##   -------
##   seqinfo: 93 sequences from hg19 genome
```

More work would be needed to isolate the AceView gene models, by
working on the `name` component of the mcols of the returned
GRanges.

## Building your own annotation resources

It has long been recognized that labs will generate their own
annotation of various forms for organisms that may not be
covered by Bioconductor's efforts.  See the
[AnnotationForge](http://www.bioconductor.org/packages/release/bioc/html/AnnotationForge.html) vignette, and post to the support site if this
is not sufficient.
