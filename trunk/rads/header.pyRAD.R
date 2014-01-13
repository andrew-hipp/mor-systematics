## handling pyRAD data
## A Hipp, 2012-09-18
## -- Sept 2012: accommodate new pyRAD format (GBS data)
## -- Dec 2012: new functions to help with exporting RAD data and blasting
## -- Jan 2013: consensus functions now use Biostrings instead of seqinr

## uncomment if not already installed:
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

## FUNCTIONS ARE BINNED INTO THE FOLLOWING FILES:
## read.pyRAD -- read files in
## write.pyRAD -- write pyRAD files out
## manipulate.pyRAD -- create consensus files
## subset.pyRAD -- subset and do some data interconversion
## swul -- a variety of functions for what used to be successive weighting under likelihood
## unique.multiPhylo -- fixes a bug in the ape function of the same name, and adds a useful attribute
## visualize.pyRAD -- generate distances and plots
## workflows.pyRAD -- workflows and wrappers around other functions

