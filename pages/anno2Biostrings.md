---
title: "Chromosomes and their substructures 2: Biostrings"
author: "Vince"
date: "March 19, 2015"
output: html_document
layout: page
toc: yes
---






# Biostrings: basic infrastructure for computing on sequences

## Construction, sets, restricted alphabets

Very large strings like chromosome sequences receive
special handling in Bioconductor.  We use a general container
class called `BString` for "big" strings that are
distringuished from R character vectors in that BStrings a) obey
different rules for copying and b) do not contain multiple
strings (see the man page for BString).  Classes `DNAString`
and `AAString` have restrictions on the characters that can be
managed in instances.


```r
library(Biostrings)
bdemo = BString("BCDEF")
ddemo = try(DNAString("BCDEF"))
cat(ddemo)
```

```
## Error in .Call2("new_XString_from_CHARACTER", classname, x, start(solved_SEW),  : 
##   key 69 (char 'E') not in lookup table
```

```r
ademo = try(AAString("BCDEF"))
```

Efficient management of multiple strings employs classes with
"Set" as suffix.

```r
ddem2 = DNAStringSet(c("ACTG", "GTCAG"))
ddem2
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]     4 ACTG
## [2]     5 GTCAG
```

The restrictions on contents of genomic strings are defined
in constant vectors in `Biostrings`.  For example

```r
AA_ALPHABET
```

```
##  [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T"
## [18] "W" "Y" "V" "U" "O" "B" "J" "Z" "X" "*" "-" "+" "."
```

```r
IUPAC_CODE_MAP
```

```
##      A      C      G      T      M      R      W      S      Y      K 
##    "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT" 
##      V      H      D      B      N 
##  "ACG"  "ACT"  "AGT"  "CGT" "ACGT"
```

## Operations

There are over 200 functions defined in the Biostrings package,
all devoted to computation on sequence data.  Here's an
example illustrating basic notions.


```r
D = DNAString("ACTGACGTACGTAGGCTAGCGATCGATATACGATATACG")
translate(D)
```

```
##   13-letter "AAString" instance
## seq: TDVRRLAIDIRYT
```

```r
codons(D)
```

```
##   Views on a 39-letter DNAString subject
## subject: ACTGACGTACGTAGGCTAGCGATCGATATACGATATACG
## views:
##      start end width
##  [1]     1   3     3 [ACT]
##  [2]     4   6     3 [GAC]
##  [3]     7   9     3 [GTA]
##  [4]    10  12     3 [CGT]
##  [5]    13  15     3 [AGG]
##  ...   ... ...   ... ...
##  [9]    25  27     3 [GAT]
## [10]    28  30     3 [ATA]
## [11]    31  33     3 [CGA]
## [12]    34  36     3 [TAT]
## [13]    37  39     3 [ACG]
```

Notice that the output of codons is printed as a `Views` instance.
This is a very efficient approach to creating references to
subsequences of a sequence, without copying any data.
