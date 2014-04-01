subset.pyRAD.loci <-
function(x, loci, taxa, format = 'DNAStringSet', reportInterval = 500, 
         nucVarType = c("verystrict", "strict", "relaxed"), use.tidyName = FALSE,
         snpsOnly = FALSE, cores = 1, ...) {
## Arguments:
##  x = pyRAD.loci object
##  loci = loci to use (by name)
##  taxa = taxa to use (by name)
##  format = only DNAStringSet export supported now
##  reportInterval = interval to use for reporting subsetting progress
##  nucVarType = either strict for requiring variability to be due to only unambiguous nucleotides (strict, verystrict) or allowing ambiguities to encode variability (relaxed), and insisting on parsimony informativeness (verystrict) or just variability (strict, relaxed) at at least one site
##  2014-03-18: added snpLocs to out to return where the snps are
  excludedNucs <- switch(nucVarType[1], verystrict = 5:17, strict = 5:17, relaxed = 15:17)
  nucThresh <- switch(nucVarType[1], verystrict = 2, strict = 1, relaxed = 1)
  if(use.tidyName) inds.vector <- tidyName(x$tips, ...) %in% tidyName(taxa, ...)
  else x$tips %in% taxa
  counter = 0
  start.time <- Sys.time()
  if(cores == 1) {
  out <- list(
    DNA = structure(vector('list', length(loci)), names = loci),
    variable = structure(logical(length(loci)), names = loci),
	ntaxa = structure(integer(length(loci)), names = loci),
	snpLocs = structure(vector('list', length(loci)), names = loci)
	)
    for(i in loci) {
	  counter <- counter + 1
	  if(counter / reportInterval - counter %/% reportInterval == 0) {
  	    message(paste('... subsetting', counter, 'of', length(loci), 
 	    '-- Estimated time remaining =', round(((Sys.time() - start.time) / counter) * (length(loci) - counter), 1), attr(Sys.time() - start.time, 'units')
  	    ))
	  }
      seq.index <- x$locus.index == i & inds.vector
	  out$DNA[[i]] <- DNAStringSet(x$seqs[seq.index])
	  names(out$DNA[[i]]) <- x$tips[seq.index]
	  out$snpLocs[[i]] <- which(apply(consensusMatrix(out$DNA[[i]])[-c(excludedNucs), ], 2, function(x) sum(x >= nucThresh) > 1))
	  out$variable[i] <- length(out$snpLocs[[i]] > 0)
	  out$ntaxa[i] <- sum(seq.index)
	  if(snpsOnly) {
	    tempDNA <- as.matrix(as.matrix(out$DNA[[i]])[, out$snpLocs[[i]]])
	    out$DNA[[i]] <- try(DNAStringSet(apply(tempDNA, 1, paste, collapse = '')))
	    if(class(out$DNA[[i]]) == 'try-error') {
	      message(paste('failed on locus', i))
		  out$DNA[[i]] <- NULL
		  }
	    } # close if
	  } # close i
	} # close cores = 1 branch
  else {
    do.it.dna <- function(i) {
	  seq.index <- x$locus.index == i & inds.vector
	  dna.temp <- DNAStringSet(x$seqs[seq.index])
	  names(dna.temp) <- x$tips[seq.index]
	  return(dna.temp)
	  } # close do.it.dna
	out <- list(DNA = mclapply(loci, do.it.dna, mc.cores = cores))
	names(out$DNA) <- loci
	out$snpLocs <- mclapply(out$DNA, function(i) which(apply(consensusMatrix(i)[-c(excludedNucs), ], 2, function(x) sum(x >= nucThresh) > 1)), mc.cores = cores)
	if(snpsOnly) {
	  do.it.snps <- function(i) {
	    tempDNA <- as.matrix(as.matrix(out$DNA[[i]])[, out$snpLocs[[i]]])
		snps.out <- try(DNAStringSet(apply(tempDNA, 1, paste, collapse = '')))
		if(class(snps.out) == 'try-error') message(paste('failed on locus', i))
		return(snps.out)
	    } # close do.it.snps
	  out$DNA <- mclapply(loci, do.it.snps, mc.cores = cores)
	  names(out$DNA) <- loci
	  out$DNA <- out$DNA[sapply(out$DNA, class) != 'try-error']
	  } # close if snpsOnly
	out$ntaxa <- sum(seq.index)
	out$variable <- sapply(out$snpLocs, function(i) length(i) > 0)
	} # close else -- the multicore branch
  class(out) <- 'subset.pyRAD.loci'
  return(out)
  }
