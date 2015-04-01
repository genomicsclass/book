---
title: "Genomic annotation in Bioconductor: Overview"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---





# Basic annotation resources and their discovery

In this document we will review Bioconductor's facilities for
handling and annotating genomic sequence.  We'll look at
reference genomic sequence, transcripts and genes, and
conclude with gene pathways.  Keep in mind that our ultimate aim
is to use annotation information to help produce reliable
interpretations of genomic experiments.  A basic objective of
Bioconductor is to make it easy to incorporate
information on genome structure and function 
into statistical analysis procedures.

## A simple hierarchy of annotation concepts

Bioconductor includes many different types of genomic annotation.
We can think of these annotation resources in a hierarchical structure.

- At the base is the reference genomic sequence for an organism.
This is always arranged into chromosomes, specified by linear
sequences of nucleotides.

- Above this is the organization of chromosomal sequence into
regions of interest.  The most prominent regions of interest are
genes, but other structures like SNPs or CpG sites are
annotated as well.  Genes have internal structure,
with parts that are transcribed and parts that are not,
and "gene models" define the ways in which
these structures are labeled and laid out in genomic coordinates.

- Above this is the organization of genes or gene products into
groups with shared structural or functional properties.  Examples
include pathways, groups of genes found together in cells, or
identified as cooperating in biological processes.

## Discovering available reference genomes

Bioconductor's collection of annotation packages brings
all elements of this hierarchy into a programmable environment.
Reference genomic sequences are managed using the infrastructure
of the Biostrings and BSgenome packages, and the `available.genomes`
function lists the reference genome build for humans and
various model organisms now available.


```r
library(Biostrings)
ag = available.genomes()
length(ag)
```

```
## [1] 71
```

```r
head(ag)
```

```
## [1] "BSgenome.Alyrata.JGI.v1"                
## [2] "BSgenome.Amellifera.BeeBase.assembly4"  
## [3] "BSgenome.Amellifera.UCSC.apiMel2"       
## [4] "BSgenome.Amellifera.UCSC.apiMel2.masked"
## [5] "BSgenome.Athaliana.TAIR.04232008"       
## [6] "BSgenome.Athaliana.TAIR.TAIR9"
```

## Reference build versions are important

Note that the genome sequence packages have long names
that include build versions.  It is very important to avoid
mixing coordinates from different reference builds.
We will show later how to convert genomic coordinates of
features between different reference builds, using the UCSC
"liftOver" utility interfaced to R in the `rtracklayer` package.

# A reference genomic sequence for H. sapiens

The reference sequence for *Homo sapiens* is acquired by installing
and attaching
a single package.  This is in contrast to downloading and parsing
FASTA files.  The package defines an object `Hsapiens`
that is the source of chromosomal sequence, but when
evaluated on its own
provides a report of the origins of the sequence data that
it contains.


```r
library(BSgenome.Hsapiens.UCSC.hg19)
```

```
## 
## Attaching package: 'BSgenome.Hsapiens.UCSC.hg19'
## 
## The following object is masked from 'package:BSgenome.Hsapiens.NCBI.GRCh38':
## 
##     Hsapiens
```

```r
Hsapiens
```

```
## Human genome
## | 
## | organism: Homo sapiens (Human)
## | provider: UCSC
## | provider version: hg19
## | release date: Feb. 2009
## | release name: Genome Reference Consortium GRCh37
## | 93 sequences:
## |   chr1                  chr2                  chr3                 
## |   chr4                  chr5                  chr6                 
## |   chr7                  chr8                  chr9                 
## |   chr10                 chr11                 chr12                
## |   chr13                 chr14                 chr15                
## |   ...                   ...                   ...                  
## |   chrUn_gl000235        chrUn_gl000236        chrUn_gl000237       
## |   chrUn_gl000238        chrUn_gl000239        chrUn_gl000240       
## |   chrUn_gl000241        chrUn_gl000242        chrUn_gl000243       
## |   chrUn_gl000244        chrUn_gl000245        chrUn_gl000246       
## |   chrUn_gl000247        chrUn_gl000248        chrUn_gl000249       
## | (use 'seqnames()' to see all the sequence names, use the '$' or '[['
## | operator to access a given sequence)
```

We acquire a chromosome's sequence using the `$` operator.

```r
Hsapiens$chr17
```

```
##   81195210-letter "DNAString" instance
## seq: AAGCTTCTCACCCTGTTCCTGCATAGATAATTGC...GGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGT
```

# The transcripts and genes for a reference sequence

The `TxDb` family of packages and data objects manages
information on transcripts and gene models.


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
```

```
## Loading required package: GenomicFeatures
```

```r
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene # abbreviate
txdb
```

```
## TxDb object:
## | Db type: TxDb
## | Supporting package: GenomicFeatures
## | Data source: UCSC
## | Genome: hg19
## | Organism: Homo sapiens
## | UCSC Table: knownGene
## | Resource URL: http://genome.ucsc.edu/
## | Type of Gene ID: Entrez Gene ID
## | Full dataset: yes
## | miRBase build ID: GRCh37
## | transcript_nrow: 82960
## | exon_nrow: 289969
## | cds_nrow: 237533
## | Db created by: GenomicFeatures package from Bioconductor
## | Creation time: 2014-09-26 11:16:12 -0700 (Fri, 26 Sep 2014)
## | GenomicFeatures version at creation time: 1.17.17
## | RSQLite version at creation time: 0.11.4
## | DBSCHEMAVERSION: 1.0
```

We can use `genes()` to get the addresses of genes using 
Entrez Gene IDs.


```r
ghs = genes(txdb)
ghs
```

```
## GRanges object with 23056 ranges and 1 metadata column:
##         seqnames                 ranges strand   |     gene_id
##            <Rle>              <IRanges>  <Rle>   | <character>
##       1    chr19 [ 58858172,  58874214]      -   |           1
##      10     chr8 [ 18248755,  18258723]      +   |          10
##     100    chr20 [ 43248163,  43280376]      -   |         100
##    1000    chr18 [ 25530930,  25757445]      -   |        1000
##   10000     chr1 [243651535, 244006886]      -   |       10000
##     ...      ...                    ...    ... ...         ...
##    9991     chr9 [114979995, 115095944]      -   |        9991
##    9992    chr21 [ 35736323,  35743440]      +   |        9992
##    9993    chr22 [ 19023795,  19109967]      -   |        9993
##    9994     chr6 [ 90539619,  90584155]      +   |        9994
##    9997    chr22 [ 50961997,  50964905]      -   |        9997
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
```

# The NCBI Entrez Gene annotation maps

Packages named org.*.eg.db collect information at the gene level
with links to location, protein product identifiers, KEGG pathway and
GO terms, PMIDs of papers mentioning genes, and to
identifiers for other annotation resources.


```r
library(org.Hs.eg.db)
```

```
## Loading required package: DBI
```

```r
keytypes(org.Hs.eg.db) # columns() gives same answer
```

```
##  [1] "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
##  [5] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"      
##  [9] "CHRLOCEND"    "ENZYME"       "MAP"          "PATH"        
## [13] "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"     
## [17] "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"    
## [21] "UNIPROT"      "GO"           "EVIDENCE"     "ONTOLOGY"    
## [25] "GOALL"        "EVIDENCEALL"  "ONTOLOGYALL"  "OMIM"        
## [29] "UCSCKG"
```

```r
head(select(org.Hs.eg.db, keys="ORMDL3", keytype="SYMBOL", 
   columns="PMID"))
```

```
##   SYMBOL     PMID
## 1 ORMDL3 11042152
## 2 ORMDL3 12093374
## 3 ORMDL3 12477932
## 4 ORMDL3 14702039
## 5 ORMDL3 15489334
## 6 ORMDL3 16169070
```

# A unified, self-describing approach

The OrganismDb packages simplify access to annotation.
Queries that succeed against TxDb, and org.[Nn].eg.db
can be directed at the OrganismDb object.


```r
library(Homo.sapiens)
```

```
## Loading required package: OrganismDbi
## Loading required package: GO.db
```

```r
class(Homo.sapiens)
```

```
## [1] "OrganismDb"
## attr(,"package")
## [1] "OrganismDbi"
```

```r
tx = transcripts(Homo.sapiens)
keytypes(Homo.sapiens)
```

```
##  [1] "GOID"         "TERM"         "ONTOLOGY"     "DEFINITION"  
##  [5] "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
##  [9] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"      
## [13] "CHRLOCEND"    "ENZYME"       "MAP"          "PATH"        
## [17] "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"     
## [21] "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"    
## [25] "UNIPROT"      "GO"           "EVIDENCE"     "GOALL"       
## [29] "EVIDENCEALL"  "ONTOLOGYALL"  "OMIM"         "UCSCKG"      
## [33] "GENEID"       "TXID"         "TXNAME"       "EXONID"      
## [37] "EXONNAME"     "CDSID"        "CDSNAME"
```

```r
columns(Homo.sapiens)
```

```
##  [1] "GOID"         "TERM"         "ONTOLOGY"     "DEFINITION"  
##  [5] "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
##  [9] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"      
## [13] "CHRLOCEND"    "ENZYME"       "MAP"          "PATH"        
## [17] "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"     
## [21] "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"    
## [25] "UNIPROT"      "GO"           "EVIDENCE"     "GOALL"       
## [29] "EVIDENCEALL"  "ONTOLOGYALL"  "OMIM"         "UCSCKG"      
## [33] "CDSID"        "CDSNAME"      "CDSCHROM"     "CDSSTRAND"   
## [37] "CDSSTART"     "CDSEND"       "EXONID"       "EXONNAME"    
## [41] "EXONCHROM"    "EXONSTRAND"   "EXONSTART"    "EXONEND"     
## [45] "GENEID"       "TXID"         "EXONRANK"     "TXNAME"      
## [49] "TXCHROM"      "TXSTRAND"     "TXSTART"      "TXEND"
```

