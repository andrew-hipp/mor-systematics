## Functions to perform successive weighting using likelihood (SWUL)
## A Hipp, Sept 2010 (ahipp@mortonarb.org)
## Oct 2013: updated for PLoS ONE paper, revised to use RAxML

require(ape)
# require(seqinr)
require(phangorn)
  
## Functions:
##  genTrees - generates trees and PAUP commands to get site likelihoods - DONE
##  getLikelihoods - takes a filename for site likelihoods from PAUP outfile and a vector of locus assignments for each site; returns a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods - DONE
##  rankLikelihoods - ranks likelihoods, counts likelihood steps, and calculates the goodness-of-fit for each locus to the optimal tree relative to the random tree set - DONE
##  swulData - generates a new dataset, cutting off either the best or worst sites by percentage

genTrees <- function(x, N = 20, filebase = 'trial', method = c('nni', 'random'), maxmoves = 3, perms = c(length(nni(x)), as.integer(100-(length(nni(x))/2)), as.integer(100-(length(nni(x))/2))), software = c('raxml', 'paup'), ...) {
  ## Arguments:
  ## x = phylo tree
  ## N = number of trees to generate per nni / spr stratum
  ## filebase = file name base; a tree file (.tre) and paup command file (.nex) will be created for both
  ## method = method for generating trees
  ## maxmoves = maximum number of rearrangements per tree for nni or spr
  ## perms = number of permutations per maxmoves
  ## ... = additional arguments to pass along to rtree.phylo or rNNI
  ## works with nni, 12 nov 10
  if(class(x) != 'phylo') stop('This function requires a phylo object as its first argument')
  if(method[1] == 'nni') {
	for(i in seq(maxmoves)) {
	  message(paste('doing maxmoves', i))
	  if(i == 1) treeset <- nni(x) 
	  else treeset <- c(treeset, rNNI(x, i, perms[i]))
      }	# end i
	} # end if	  
  else if(method[1] == 'random') treeset = rtree.phylo(x, N, ...)
  if(software[1] == 'raxml') {
    message('writing raxml')
	treeset[2:(length(treeset) + 1)] <- treeset[1:length(treeset)]
	treeset[[1]] <- x
	write.tree(treeset, file = paste(filebase, '.trees.tre', sep = ''))
    # write.tree(x, file = paste(filebase, '.optimal.tre', sep = '')) ## no longer separating optimal from full trees
    }
  if(software[1] == 'raxml') {
    message('RAxML chosen as analysis software. Currently, you just need to run this on your own to get the site likelihoods.Try something like this:\n
	/home/andrew/code/raxml/standard-RAxML-7.7.2/raxmlHPC-PTHREADS-SSE3 -f g -T 10 -s d6m10.phy -m GTRGAMMA -z analysis.d6m10/RAxML_bestTree.d6m10.out -n d6m10.phy.reduced.siteLnL')
	}
  return(treeset)
  }

get.raxml.siteLikelihoods <- function(x)  {
    lnL <- readLines(x)
    lnL <- strsplit(lnL[2:length(lnL)], "\t")
    names(lnL) <- unlist(lapply(lnL, function(x) x[1]))
    lnL <- unlist(lapply(lnL, function(x) x[2]))
    lnL <- t(sapply(strsplit(lnL, " "), as.numeric)) # this is a matrix with trees as rows, site lnL as columns
    return(lnL)
	}
  
 
getLikelihoods.raxml <- function(dat, lnL = NA, info = NA, missingSites = NA, treeScoreFile = choose.files(multi = FALSE, caption = "Select RAxML site likelihoods file for trees"), infoFile = choose.files(multi = FALSE, caption = 'Select RAxML info file for analysis')) {
  ## this version of the getLikelihood function tosses 
  ## ARGUMENTS:
  ##   treeScoreFile - filenames for site likelihoods from RAxML outfile and a vector of locus assignments for each site
  ##   dat - pyRAD.loci object
  ## VALUE: a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## works: 2013-10-29
  
  ## 1. read data -- result is a list of site likelihoods, one per tree
  if(is.na(lnL[1])) lnL <- get.raxml.siteLikelihoods(treeScoreFile)
  ## following appears not to be needed...
  # if(is.na(info[1])) {
	# info <- readLines(infoFile)
	# missing.sites <- grep('Alignment column', info)
	# missing.sites <- as.numeric(sapply(info[missing.sites], function(x) strsplit(x, ' ')[[1]][5]))
	# lnL.2 <- matrix(NA, nrow = dim(lnL)[1], ncol = dim(lnL)[2] + length(missing.sites))
	# lnL.2[, -missing.sites] <- lnL
	# lnL.2[, missing.sites] <- 0
	# message(paste('Taxa read:', dim(lnL.2)[1]))
	# message(paste('Sites read:', dim(lnL.2)[2]))
	# message(paste('Target number of sites:', sum(dat$radSummary$locus.lengths)))
	# lnL <- lnL.2 # this is a matrix of site likelihoods, with 0s for undefined columns; okay because we care about rank order, not absolute lnL
	# }
	
  # 1. book-keeping, so we can find trees and clusters
  message("Doing bookkeeping...")
  nLoci <- length(dat$radSummary$locus.lengths)
  loc.ranges <- cbind(c(1, cumsum(dat$radSummary$locus.lengths) + 1)[1:nLoci], cumsum(dat$radSummary$locus.lengths))
  clusterBP <- apply(loc.ranges, 1, function(x) x[1]:x[2]) # this is a list of numeric vectors, one per locus, giving all bp positions for that locus
  
  # 2. locus likelihoods for each locus and tree
  locusScores = matrix(0, nrow = dim(lnL)[1], ncol = nLoci, dimnames = list(row.names(lnL), names(clusterBP)))
  treeScores <- apply(lnL, 1, sum)
  for(treeNumber in c(1:length(treeScores))) {
	message(paste("Summing locus likelihoods on tree", treeNumber))
	# I think it would be more efficient to vectorize the next row by creating a vector of locus numbers, then splitting
	for(locusNumber in seq(nLoci)) locusScores[treeNumber, locusNumber] <- sum(lnL[treeNumber, clusterBP[[locusNumber]]])
	} #close treeNumber
  row.names(locusScores)[treeScores == min(treeScores)] <- 'best'
  out <- list(locusScores = locusScores, treeScores = treeScores, clustersPresent = NA)
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