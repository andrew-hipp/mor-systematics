getLikelihoods.paup <- function(rtreeScoreFile = choose.files(multi = FALSE, caption = "Select score file for random trees"), 
                           bestTreeScoreFile = choose.files(multi = FALSE, caption = "select score file for best tree"),
						   locusVector = clusterLens) {
  ## takes a filename for site likelihoods from PAUP outfile and a vector of locus assignments for each site; returns a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## note that the log-likelihood scores written to the file by PAUP are negative of the true lnL
  
  # 1. book-keeping, so we can find trees and clusters

  message("Doing bookkeeping...")
  clusterBP <- clusterSites(locusVector)
  bestTreeScores <- read.delim(bestTreeScoreFile) # data.frame with site likelihoods for best tree
  treeScores <- read.delim(rtreeScoreFile) # data.frame with all site likelihoods and tree likelihoods for rtrees
  endOfEachTree <- which(!is.na(rtreeScores$Tree)) # last row of each tree in rtreeScores
  begOfEachTree <- c(1, endOfEachTree[1:length(endOfEachTree) - 1] + 1) # first row of each tree in rtreeScores
  treeRows <- lapply(1:length(endOfEachTree), function(x) seq(from = begOfEachTree[x], to = endOfEachTree[x])) # list of rows for each tree in rtreeScores
  
  message("Finding clusters present...")
  clustersPresent <- which(sapply(clusterBP, function(x) sum(x %in% rtreeScores[treeRows[[1]], 'Site']) > 0))
  
  # 2. locus likelihoods for each locus and tree
  locusScores = matrix(0, nrow = length(treeRows) + 1, ncol = length(clustersPresent), dimnames = list(c(1:length(treeRows), 'best'), as.character(clustersPresent)))
  treeScores <- c(rtreeScores[endOfEachTree, 'X.lnL'], bestTreeScores[!is.na(bestTreeScores$Tree),'X.lnL'])
  names(treeScores) <- dimnames(locusScores)[[1]]  
  for(treeNumber in c(1:length(treeRows), -9)) {
    if(treeNumber == -9) {
	  tr <- bestTreeScores
	  treeNumber <- 'best'
	  message("Summing locus likelihoods on the best (ML) tree")
	  }
	else {
	  tr <- rtreeScores[treeRows[[treeNumber]], ]
	  message(paste("Summing locus likelihoods on tree", treeNumber))
	  }
	for(locusNumber in clustersPresent) {
	  locusScores[treeNumber, as.character(locusNumber)] <- 
	    sum(tr$X.lnL.1[which(tr$Site %in% clusterBP[[locusNumber]])])
	  } #close locusNumber
	} #close treeNumber
	
  out <- list(locusScores = locusScores, treeScores = treeScores, clustersPresent = clustersPresent)
  class(out) <- 'swulLikelihoods'
  return(out)
  }
  
 