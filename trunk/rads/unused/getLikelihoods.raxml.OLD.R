getLikelihoods.raxml <- function(dat, lnL = NA, method = c('by.site','by.locus'), missingSites = NA, which.loci, treeScoreFile = choose.files(multi = FALSE, caption = "Select RAxML site likelihoods file for trees")) {
  ## ARGUMENTS:
  ##   treeScoreFile - filenames for site likelihoods from RAxML outfile and a vector of locus assignments for each site
  ##   dat - pyRAD.loci object
  ## VALUE: a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## works: 2013-10-29
  ## 2014-01-21: added a 'by.site' option to accommodate the locus-by-locus analysis implemented now
  ## 2014-01-21, later in the day: tossing this function into the unused folder... it's unwieldy and unnecessary now.
  
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
  