subset.pyRAD.loci <- function(x, loci, taxa, format = 'DNAStringSet', ...) {
  ## only DNAStringSet export supported now
  out <- list(
    DNA = structure(vector('list', length(loci)), names = loci),
    variable = structure(logical(length(loci)), names = loci)
	)
  inds.vector <- x$tips %in% taxa
  for(i in loci) {
    seq.index <- x$locus.index == i & inds.vector
	out$DNA[[i]] <- DNAStringSet(x$seqs[seq.index])
	names(out$DNA[[i]]) <- x$tips[seq.index == i]
	out$variable[i] <- any(apply(consensusMatrix(out$DNA[[i]])[-c(15:17), ], 2, function(x) sum(x > 0) > 1))
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

consensus.pyRAD <- function(pyIn, from = NA, to = NA, fastaNames = T, writeFile = 'rads.con.1_100.txt', ...) {
## use seqinr to generate a consensus sequence for each pyRAD locus
## 2013-01-04: updated to use Biostrings, which works better -- deleted arguments: method = 'majority', threshold = 0.001
  if(class(pyIn) != "pyRAD.loci") stop("pyRAD input required, from read.pyRAD")
  # require(seqinr)
  require(Biostrings)
  allLoci <- unique(as.character(pyIn$locus.index))
  allLoci <- allLoci[!allLoci == ""] ## this should probably part of read.pyRAD
  if(!is.na(from)) allLoci <- allLoci[from:to]
  seqs <- as.character(pyIn$seqs)
  loc.index <- as.character(pyIn$locus.index)
  out <- character(0)
  for (i in allLoci) {
    message(i) #only for debugging
	out <- c(out, consensusString(DNAStringSet(gsub("-", "N", seqs[loc.index == i])), ...))
	}
  # for(i in allLoci) out <- c(out, paste(consensus(str2mat(seqs[loc.index == i]), method, threshold), collapse = ""))
  if(fastaNames) allLoci <- paste(">", allLoci, sep = "")
  names(out) <- allLoci
  if(!is.na(writeFile)) write.table(out, writeFile, sep = "\n", quote = F, col.names = F)
  return(out)
  }

locus.picker <- function(pyDat, minThreshold = 3, inds = row.names(pyDat$radSummary$inds.mat)) {
  out = names(which(colSums(pyDat$radSummary$inds.mat[inds, ]) >= minThreshold))
  }
  
