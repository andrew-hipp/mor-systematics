parsimonyInformBipartition <- function(dat, bipartition, silent = TRUE) {
## calculates the parsimony informativeness of a locus / dataset for one bipartition
## dat = an object of class subset.pyRAD.loci
## bipartition = list of two sets of tips representing the bipartition

## 

## calculate per site as ((number tips with dominant nucleotide for set 1) 
##                       + (number of tips with dominant nucleotide for set 2)
##                       - (total tips if dominant in set 1 is the same as dominant in set 2))
##                       / total number of tips

## per locus, calculate so that it ranges from 0 to 1, where a zero is only when there are no variable sites, 
## and a 1 is when one or all variable nucleotides are perfect. Then, do stats either over all SNPs or just for first SNP in each locus

  dat.tips = XXX
  if(length(intersect(unique(unlist(bipartition)), dat.tips)) != length(unique(c(dat.tips), unlist(bipartition)))) {
    if(!silent) warning('data and bipartitions not identical; pruning as needed')
	