## Functions to perform successive weighting using likelihood (SWUL)
## A Hipp, Sept 2010 (ahipp@mortonarb.org)

require(ape)
require(seqinr)

## Functions:
##  genTrees - generates trees and PAUP commands to get site likelihoods - DONE
##  getLikelihoods - takes a filename for site likelihoods from PAUP outfile and a vector of locus assignments for each site; returns a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods - DONE
##  rankLikelihoods - ranks likelihoods, counts likelihood steps, and calculates the goodness-of-fit for each locus to the optimal tree relative to the random tree set - DONE
##  swulData - generates a new dataset, cutting off either the best or worst sites by percentage

clusterSites <- function(locusVector) {
  clusterEnds <- cumsum(locusVector) # a vector of the position for each cluster first bp
  clusterBegs <- clusterEnds - locusVector + 1 # a vector of the position for each cluster last bp
  out <- lapply(1:length(locusVector), function(x) seq(from = clusterBegs[x], to = clusterEnds[x])) # list of base pairs in each cluster
  return(out)
  }

genTrees <- function(x, N = 20, filebase = 'trial', method = c('nni', 'random'), maxmoves = 3, perms = c(length(nni(x)), 100, 100), ...) {
  ## Arguments:
  ## x = phylo tree
  ## N = number of trees to generate
  ## filebase = file name base; a tree file (.tre) and paup command file (.nex) will be created for both
  ## method = method for generating trees
  ## maxmoves = maximum number of rearrangements per tree for nni or spr
  ## perms = number of permutations per maxmoves
  ## ... = additional arguments to pass along to rtree.phylo or rNNI
  ## works with nni, 12 nov 10
  treeset = vector('list', sum(perms[seq(maxmoves)]))
  if(method[1] == 'nni') {
    require(phangorn)
	for(i in seq(maxmoves)) {
	  if(i == 1) treesetTemp <- nni(x) 
	    else treesetTemp <- rNNI(x, i, perms[i])
	  for(j in 1:length(treesetTemp)) treeset[[j + sum(perms[0:(i-1)])]] <- treesetTemp[[j]]
      }	# end i
	} # end if	  
  else if(method[1] == 'random') treeset = rtree.phylo(x, N, ...)
  class(treeset) <- 'multiPhylo'
  write.nexus(treeset, file = paste(filebase, '.rtrees.tre', sep = ''))
  write.nexus(x, file = paste(filebase, '.optimal.tre', sep = ''))
  
  ## write paup block
  paup.out <- file(description = paste(filebase, '.paupBlock.nex', sep = ''), open = "w")
  writeLines("#NEXUS", paup.out)
  writeLines(paste("[PAUP block written using genTrees function in R,", date(), "]"), paup.out)
  writeLines("\nBEGIN PAUP;\n", paup.out)
  writeLines("  exclude constant /only;", paup.out)
  writeLines("  increase = auto autoclose = yes;", paup.out)
  writeLines("  cleartrees;", paup.out)
  writeLines(paste("  gettrees file =", paste(filebase, '.optimal.tre', sep = ''), ";"), paup.out)
  writeLines(paste("  lscores all / sitelikes = yes scorefile =", paste(filebase, '.optimal.scores', sep = ''), ";"), paup.out)
  writeLines("\n  cleartrees;", paup.out)
  writeLines(paste("  gettrees file =", paste(filebase, '.rtrees.tre', sep = ''), ";"), paup.out)
  writeLines(paste("  lscores all / sitelikes = yes scorefile =", paste(filebase, '.rtrees.scores', sep = ''), ";"), paup.out)
  writeLines("\nend;", paup.out)
  close(paup.out)
  return(treeset)
  }

getLikelihoods <- function(rtreeScoreFile = choose.files(multi = FALSE, caption = "Select score file for random trees"), 
                           bestTreeScoreFile = choose.files(multi = FALSE, caption = "select score file for best tree"),
						   locusVector = clusterLens) {
  ## takes a filename for site likelihoods from PAUP outfile and a vector of locus assignments for each site; returns a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## note that the log-likelihood scores written to the file by PAUP are negative of the true lnL
  
  # 1. book-keeping, so we can find trees and clusters

  message("Doing bookkeeping...")
  clusterBP <- clusterSites(locusVector)
  bestTreeScores <- read.delim(bestTreeScoreFile) # data.frame with site likelihoods for best tree
  rtreeScores <- read.delim(rtreeScoreFile) # data.frame with all site likelihoods and tree likelihoods for rtrees
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
  
rankLikelihoods <- function(x) {
  # ranks likelihoods, counts likelihood steps, and calculates the goodness-of-fit for each locus to the optimal tree relative to the random tree set
  # Arguments:
  #  x = a swulLikelihoods object
  rtreeVector <- dimnames(x$locusScores)[[1]][dimnames(x$locusScores)[[1]] != 'best'] # I'm sure there's a shorter way to do this....
  percentiles <- apply(x$locusScores, 2, function(x) sum(x[rtreeVector] > x[length(x)]) / (length(x) - 1)) # yields the percent of loci that have a higher likelihood on the ML tree than on the random set
  # categories <- apply(x$locusScores[rtreeVector, ], 2, function(x) sort(unique(x)))
  class(percentiles) <- 'rankLikelihoods'
  return(percentiles)
  }
  
besttree <- function(x, allow.doubles = TRUE, summaryByTree = FALSE, takeMax = FALSE) {
  # returns best tree for each locus
  if(class(x) != 'swulLikelihoods') warning('Not a swulLikelihoods object... what are you thinking?')
  if(!takeMax) {
  if(allow.doubles) bestTreeList <- apply(x[[1]], 2, function(x) which(x == min(x)))
  else bestTreeList = apply(x[[1]], 2, function(x) {
                     temp <- which(x == min(x))
                     ifelse(length(temp) == 1, temp, NA)
                     })
  } # end if not takeMax
  
  if(takeMax) {
  if(allow.doubles) bestTreeList <- apply(x[[1]], 2, function(x) which(x == max(x)))
  else bestTreeList = apply(x[[1]], 2, function(x) {
                     temp <- which(x == max(x))
                     ifelse(length(temp) == 1, temp, NA)
                     })
  } # end if takeMax
  
  if(summaryByTree) {
    treeSummary <- numeric(dim(x)[1]) 
    for(i in seq(length(treeSummary))) treeSummary[i] <- sum(sapply(bestTreeList, function(x) i %in% x))
    out <- list(bestTreeList = bestTreeList, treeSummary <- treeSummary)
    return(out)
    }
  else return(bestTreeList)
  }

plot.rankLikelihoods <- function(x) plot(sort(x), type = 'l')

swulData <- function(fasta = RADdat, lnlRanks = NULL, locusVector, filename = "swulRadsOut.nex", lociToRemove = names(which(rtree.200.likelihoods.best == 201)), locithreshold = 1, reverse = TRUE, excludeFile = TRUE) {
  # creates a new dataset, clipping off either the upper cut of loci (best or most poorly fit by the data)
  # arguments:
  #  fasta = fasta file in SeqFastadna format (from seqinr package)
  #  lnlRanks = object from rankLikelihoods function
  #  locusVector = numeric vector with number of sites in each locus, in same order as fasta file
  #  reverse = if T, then reverse successive weighting (best sites are removed); if F, then successive weighting (worst sites are removed)
  #  excludeFile = if T, then write a RAxML exclude file rather than a full nexus file
  #  updated 14 apr 11 to allow explicit subsetting -- note that the way the loci are saved, I am using the names vector to select the actual loci
  clusterBP <- clusterSites(locusVector)
  if(is.na(lociToRemove[1]))
  { if(!reverse) lnlRanks <- 1-lnlRanks
   lociToRemove <- names(lnlRanks)[lnlRanks>=threshold]
   }
  if(excludeFile) {
    cat(unlist(clusterBP[-as.numeric(lociToRemove)]), sep = " ", file = filename)
	return(lociToRemove)
	} #end if
  else {
    for(i in 1:length(fasta)) fasta[[i]] <- fasta[[i]][unlist(clusterBP[-as.numeric(lociToRemove)])]
    write.fasta.nexus(fasta, file.out = filename)
	return(fasta)
	} #end else
  }