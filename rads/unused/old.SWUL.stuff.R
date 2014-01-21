## Functions to perform successive weighting using likelihood (SWUL)
## A Hipp, Sept 2010 (ahipp@mortonarb.org)
## Oct 2013: updated for PLoS ONE paper, revised to use RAxML
## Jan 2014: updates plotting for PLoS ONE reviews, updated to consider tree conclusiveness
## January 2014: SWUL is no longer an appropriate name for this stuff, as it's not successive weighting at all but rather data exploration;
##   -- also, many functions are just dross now and have been moved to "unused" folder
  
## Functions:
##  genTrees - generates trees and PAUP commands to get site likelihoods - DONE
##  getLikelihoods - takes a filename for site likelihoods from PAUP outfile and a vector of locus assignments for each site; returns a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods - DONE
##  rankLikelihoods - ranks likelihoods, counts likelihood steps, and calculates the goodness-of-fit for each locus to the optimal tree relative to the random tree set - DONE
##  swulData - generates a new dataset, cutting off either the best or worst sites by percentage

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

diff.swulLikelihoods <- function(x, ...) apply(x$locusScores, 2, function(z) abs(diff(range(z), ...)))