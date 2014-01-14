subset.pyRAD.loci <- function(x, loci, taxa, format = 'DNAStringSet', reportInterval = 20, ...) {
  ## only DNAStringSet export supported now
  out <- list(
    DNA = structure(vector('list', length(loci)), names = loci),
    variable = structure(logical(length(loci)), names = loci),
	ntaxa = structure(integer(length(loci)), names = loci)
	)
  inds.vector <- x$tips %in% taxa
  counter = 0
  start.time <- Sys.time()
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
	out$variable[i] <- any(apply(consensusMatrix(out$DNA[[i]])[-c(15:17), ], 2, function(x) sum(x > 0) > 1))
	out$ntaxa[i] <- sum(seq.index)
	}
  return(out)
  }
  
filter.by <- function(dat, taxa, threshold = 'all') {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat.mat <- dat$radSummary$inds.mat[taxa, ]
    else dat.mat <- dat$inds.mat[taxa, ] # this is the case if you pass in a summary.pyRAD.loci object
  if(threshold == 'all') threshold <- length(taxa)
  return(names(which(apply(dat.mat, 2, sum) == threshold)))
  }

locus.picker <- function(pyDat, minThreshold = 3, inds = row.names(pyDat$radSummary$inds.mat)) {
  out = names(which(colSums(pyDat$radSummary$inds.mat[inds, ]) >= minThreshold))
  }
  
