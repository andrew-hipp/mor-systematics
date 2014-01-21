match.lnL.to.trees <- function(directory = 'getwd()', lnLprefix = 'RAxML_info.', lnLsuffix = '.lnL', treeIndexFile = 'tree.index.lines.txt', locus.names = NULL, ...) {
  treeIndex <- read.delim(paste(directory, '/', treeIndexFile, sep = ''), as.is = T, header = F, row.names = 1)
  if(is.null(locus.names)) locus.names <- row.names(treeIndex)
  lnL.list <- lapply(paste(directory, '/', lnLprefix, locus.names, lnLsuffix, sep = ''), get.raxml.treeLikelihoods)
  names(lnL.list) <- locus.names
  out.mat <- matrix(NA, nrow = length(locus.names), ncol = dim(treeIndex)[2], dimnames = list(locus.names, NULL))
  for(i in locus.names) {
    names(lnL.list[[i]]) <- unique(as.character(treeIndex[i,]))
	out.mat[i, ] <- lnL.list[[i]][as.character(treeIndex[i,])]
	}
  return(out.mat)
  }

get.raxml.siteLikelihoods <- function(x)  {
## gets likelihoods from the RAxML_perSiteLLs file
    lnL <- readLines(x)
    lnL <- strsplit(lnL[2:length(lnL)], "\t")
    names(lnL) <- unlist(lapply(lnL, function(x) x[1]))
    lnL <- unlist(lapply(lnL, function(x) x[2]))
    lnL <- t(sapply(strsplit(lnL, " "), as.numeric)) # this is a matrix with trees as rows, site lnL as columns
    return(lnL)
	}

get.raxml.treeLikelihoods <- function(x) {
## gets likelihoods from the RAxML_info file
	fileIn <- readLines(x)
	out <- as.double(sapply(grep("Tree ", fileIn, value=T), function(x) strsplit(x, ": ")[[1]][2]))
	names(out) <- as.character(1:length(out))
	out
	}

getLikelihoods.raxml <- function(dat, lnL = NA, method = c('by.locus', 'by.site'), missingSites = NA, which.loci, treeScoreFile = choose.files(multi = FALSE, caption = "Select RAxML site likelihoods file for trees")) {
  ## ARGUMENTS:
  ##   treeScoreFile - filenames for site likelihoods from RAxML outfile and a vector of locus assignments for each site
  ##   dat - pyRAD.loci object
  ## VALUE: a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## works: 2013-10-29
  ## 2014-01-21: added a 'by.site' option to accommodate the locus-by-locus analysis implemented now
  
  ## 1. read data -- result is a list of site likelihoods, one per tree
  if(method[1] == 'by.locus') {
    ## FILL IN
	}
  if(is.na(lnL[1])) lnL <- get.raxml.siteLikelihoods(treeScoreFile)
	
  # 1. book-keeping, so we can find trees and clusters
  message("Doing bookkeeping...")
  nLoci <- length(which.loci)
  loc.ranges <- cbind(c(1, cumsum(dat$radSummary$locus.lengths[which.loci]) + 1)[1:nLoci], cumsum(dat$radSummary$locus.lengths[which.loci]))
  row.names(loc.ranges) <- which.loci # otherwise locus names are offset by one
  clusterBP <- apply(loc.ranges, 1, function(x) x[1]:x[2]) # this is a list of numeric vectors, one per locus, giving all bp positions for that locus
  
  # 2. locus likelihoods for each locus and tree
  locusScores = matrix(0, nrow = dim(lnL)[1], ncol = nLoci, dimnames = list(row.names(lnL), names(clusterBP)))
  treeScores <- rowSums(lnL)
  for(treeNumber in c(1:length(treeScores))) {
	message(paste("Summing locus likelihoods on tree", treeNumber))
	# I think it would be more efficient to vectorize the next row by creating a vector of locus numbers, then splitting
	for(locusNumber in seq(nLoci)) locusScores[treeNumber, locusNumber] <- sum(lnL[treeNumber, clusterBP[[locusNumber]]])
	} #close treeNumber
  row.names(locusScores)[treeScores == max(treeScores)] <- names(treeScores)[treeScores == max(treeScores)] <- 'best'
  out <- list(locusScores = locusScores, treeScores = treeScores, locus.total = dat$radSummary$num.inds.per.locus)
  class(out) <- 'swulLikelihoods'
  return(out)
  }
  