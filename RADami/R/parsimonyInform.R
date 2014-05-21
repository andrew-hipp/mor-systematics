parsimonyInformBipartition <- function(dat, bipartition, return.option = c('mean', 'first', 'all'), use.tidyNames = TRUE, cores = 2) {
## calculates the parsimony informativeness of a locus / dataset for one bipartition
## Arguments:
##   dat = an object of class subset.pyRAD.loci
##   bipartition = list of two sets of tips representing the bipartition
## currently conservative -- only uses unambiguous sites
## 

## calculate per site as ((number tips with dominant nucleotide for set 1) 
##                       + (number of tips with dominant nucleotide for set 2)
##                       - (total tips if dominant in set 1 is the same as dominant in set 2))
##                       / total number of tips
require(plyr)

## per locus, calculate so that it ranges from 0 to 1, where a zero is only when there are no variable sites, 
## and a 1 is when one or all variable nucleotides are perfect. Then, do stats either over all SNPs or just for first SNP in each locus
  if(!"subset.pyRAD.loci" %in% class(dat)) stop("this function requires an object of class subset.pyRAD.loci") 
  if(use.tidyNames) bipartition = lapply(bipartition, tidyNames)
  nucs = c('a', 'g', 'c', 't', 'A', 'G', 'C', 'T')
  ## nasty embedded function
  mat.stats <- function(datmat, bip) { 
	mat <- datmat[row.names(datmat) %in% bip, ]
	mat.sums <- lapply(apply(mat, 2, count), function(x) x[x$x %in% nucs, ])
	dom.mat <- cbind(do.call(rbind, lapply(mat.sums, function(y) y[which(y$freq == max(y$freq)), ])), total = sapply(mat.sums, function(x) sum(x$freq)))
	return(dom.mat)
	}
  do.it <- function(workingMat, option = c('mean', 'first', 'all')) {
    variable <- apply(as.character(phyDat(as.matrix(workingMat))),2, function(x) length(unique(x[x %in% nucs]))) > 1
    if(sum(variable) == 0) return(0) # even if we get past this without returning 0, there may be columns that have ambiguities
	if(use.tidyNames) row.names(workingMat) <- tidyNames(row.names(workingMat))
    dom.mat1 <- mat.stats(workingMat, bipartition[[1]])
	dom.mat2 <- mat.stats(workingMat, bipartition[[2]])
	sameDomFactor <- ifelse(dom.mat1$x == dom.mat2$x, -1 * (dom.mat1$total + dom.mat2$total), 0)
	stat <- (dom.mat1$freq + dom.mat2$freq - sameDomFactor) / dom.mat1$total + dom.mat2$total
	out <- switch(option[1], mean = mean(stat), first = stat[1], all = stat)
	return(out)
	}
  out <- mclapply(dat, return.option[1], mc.cores = cores)
  out