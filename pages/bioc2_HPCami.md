---
layout: page
title: "Architecture: considerations on high performance computing with Bioconductor"
---





## Distributed concurrent computing with the Bioconductor AMI

Machines with lots of cores and lots of memory are common but
are relatively expensive and can be sources of contention.
For some problems we can code for deployment of large numbers of
of separate small machines to get the solution.  This introduces
managerial concerns, both for system control and communication and
for the management of results.  We will illustrate how the
BatchJobs package allows us to distribute work to members of a
cluster that we create, replicating a virtual machine to
as many workers as we would like to use.  The virtual machine
we will use is the Bioconductor Amazon Machine Instance (AMI).
You will need credentials (an access key and secret access key)
to use EC2, but you can set this up in Amazon's free tier so that
funds are not committed until a certain amount of free computation
resource has been consumed.  This course assumes that you are
able to interact with Amazon EC2 and will not provide tutorial
background on this.

The documentation on setting up the AMI for your use is very
extensive.  There are two basic approaches.  One uses the
Amazon EC2 dashboard to configure and deploy the cluster.
The one I will demonstrate uses a set of utilities developed
at MIT called [StarCluster](http://star.mit.edu/cluster).
See [Dan Tenenbaum's notes](http://bioconductor.org/help/bioconductor-cloud-ami/) for details on how to configure
starcluster for use with the AMI.  

### Checking that we use slaves when we ask

It is a good idea to verify that the cluster is actually
functioning.  We'll demonstrate the use of BatchJobs to
do this.  With BatchJobs we define a registry that will
hold information about our cluster requests and results.

```r
library(BatchJobs)  # transcript will not reflect AMI usage
```

```
## Loading required package: BBmisc
```

```
## 
## Attaching package: 'BBmisc'
```

```
## The following object is masked from 'package:Biostrings':
## 
##     collapse
```

```
## The following object is masked from 'package:IRanges':
## 
##     collapse
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     normalize
```

```
## Sourcing configuration file: '/Library/Frameworks/R.framework/Versions/3.3/Resources/library/BatchJobs/etc/BatchJobs_global_config.R'
```

```
## Sourcing configuration file: '/Users/stvjc/.BatchJobs.R'
```

```
## BatchJobs configuration:
##   cluster functions: Interactive
##   mail.from: 
##   mail.to: 
##   mail.start: none
##   mail.done: none
##   mail.error: none
##   default.resources: 
##   debug: FALSE
##   raise.warnings: FALSE
##   staged.queries: TRUE
##   max.concurrent.jobs: Inf
##   fs.timeout: NA
```

```r
reg1 = makeRegistry("reg1")
```

```
## Loading registry: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/registry.RData
```

```r
reg1
```

```
## Job registry:  reg1 
##   Number of jobs:  8 
##   Files dir: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files 
##   Work dir: /Users/stvjc/Teaching/EDX/labs/biocadv_6x 
##   Multiple result files: FALSE 
##   Seed: 628319000 
##   Required packages: BatchJobs
```
We'll define a function that asks the node on which it is running
what its name is.

```r
myfun = function(x) system("hostname", intern=TRUE)
myfun
```

```
## function(x) system("hostname", intern=TRUE)
```
We will now map this function to the cluster via
the registry.  The idea is that we will iterate over some
set of control values, and this defines a set of tasks
that will be passed to slaves.

```r
batchMap( reg1, myfun, 1:8 )
```

```
## Error in batchMap(reg1, myfun, 1:8): Registry is not empty!
```
We then submit the jobs.  We can submit them all at once or
can submit just a few to test.

```r
submitJobs(reg1, 1:3)
```

```
## Syncing registry ...
```

```
## Saving conf: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/conf.RData
```

```
## Submitting 3 chunks / 3 jobs.
```

```
## Cluster functions: Interactive.
```

```
## Auto-mailer settings: start=none, done=none, error=none.
```

```
## Writing 3 R scripts...
```

```
## 
## Attaching package: 'checkmate'
```

```
## The following object is masked from 'package:Biobase':
## 
##     anyMissing
```

```
## Loading registry: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/registry.RData
```

```
## Loading conf:
```

```
## 2016-02-29 10:00:17: Starting job on node PC001844.local.
```

```
## Auto-mailer settings: start=none, done=none, error=none.
```

```
## Setting work dir: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## ########## Executing jid=1 ##########
```

```
## Timestamp: 2016-02-29 10:00:17
```

```
## Setting seed: 628319000
```

```
## Writing result file: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/jobs/01/1-result.RData
```

```
## 2016-02-29 10:00:17: All done.
```

```
## Setting work back to: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## Memory usage according to gc:
```

```
## Loading registry: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/registry.RData
```

```
## Loading conf:
```

```
## 2016-02-29 10:00:18: Starting job on node PC001844.local.
```

```
## Auto-mailer settings: start=none, done=none, error=none.
```

```
## Setting work dir: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## ########## Executing jid=2 ##########
```

```
## Timestamp: 2016-02-29 10:00:18
```

```
## Setting seed: 628319001
```

```
## Writing result file: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/jobs/02/2-result.RData
```

```
## 2016-02-29 10:00:18: All done.
```

```
## Setting work back to: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## Memory usage according to gc:
```

```
## Loading registry: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/registry.RData
```

```
## Loading conf:
```

```
## 2016-02-29 10:00:18: Starting job on node PC001844.local.
```

```
## Auto-mailer settings: start=none, done=none, error=none.
```

```
## Setting work dir: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## ########## Executing jid=3 ##########
```

```
## Timestamp: 2016-02-29 10:00:18
```

```
## Setting seed: 628319002
```

```
## Writing result file: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg1-files/jobs/03/3-result.RData
```

```
## 2016-02-29 10:00:19: All done.
```

```
## Setting work back to: /Users/stvjc/Teaching/EDX/labs/biocadv_6x
```

```
## Memory usage according to gc:
```

```
## Sending 3 submit messages...
## Might take some time, do not interrupt this!
```

### Dealing with a missing package

We'd like to demonstrate distributed computation with the RNA-seq
data from the Zarnack paper.  Again we'll just work with the
BAM files of alignments to chr14.  However, this package was
not installed with the AMI.  This highlights the virtues of
defining persistent storage in the form of an EBS.  We have
configured starcluster to be aware of our EBS volume that is
named '/freshdata', and all the slaves have access to it.

We installed RNAseqData.HNRNPC.bam.chr14 in /freshdata/R_LIB and
we will need to set .libPaths appropriately to ensure access.

### Setup for counting reads on the cluster

We need to define a new registry and a new function.  Our aim
is to do the overlap summaries on the slaves, one per BAM file.

We'll need our table of regions for counting reads:

```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) = "chr14"
ebg = reduce(exonsBy(txdb, by="gene")) # don't double count
save(ebg, file="/freshdata/ebg.rda") #make visible to all slaves
```
We'll define the registry so that GenomicAlignments is
attached with each slave's session.

```r
reg2 = makeRegistry("reg2", packages=c("GenomicRanges", "GenomicAlignments"))
```

```
## Loading registry: /Users/stvjc/Teaching/EDX/labs/biocadv_6x/reg2-files/registry.RData
```
The function that counts reads must use .libPaths() appropriately.

```r
reader = function(x) {
 .libPaths(c(.libPaths(), "/freshdata/R_LIB"))
 load("/freshdata/ebg.rda")
 library(RNAseqData.HNRNPC.bam.chr14)
 fn = RNAseqData.HNRNPC.bam.chr14_BAMFILES[x]
 summarizeOverlaps(ebg, fn) 
}
```
Once all this is in place, we use

```r
batchMap(reg2, reader, 1:8)
submitJobs(reg2)
showStatus(reg2)
waitForJobs(reg2) # for example
final = do.call(cbind, c(loadResults(reg2), deparse.level=1))
```
and we have a SummarizedExperiment with all the reads counted
concurrently on the slaves.


