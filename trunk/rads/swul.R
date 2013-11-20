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

plot.swulLikelihoods <- function(x, scalar = 2, percentile = c(0.025, 0.975), output = c('jpg'), bw.scalar = 1, scale.by = c('numTaxa','sd', 1), ...) {
## each tree is a data point
## X: tree likelihood
## Y: 
  X <- x$treeScores
  star.best <- which(names(X) == 'best')
  bestTreeList <- apply(x$locusScores, 2, function(z) which(z > quantile(z, percentile[2])))
  bestTree <- unlist(bestTreeList)
  worstTreeList <- apply(x$locusScores, 2, function(z) which(z < quantile(z, percentile[1])))
  worstTree <- unlist(worstTreeList)
  Y.best <- sapply(1:length(X), function(z) sum(bestTree == z))
  Y.worst <- sapply(1:length(X), function(z) sum(worstTree == z))
  if(scale.by[1] == 'sd') dotSizes.best <- dotSizes.worst <- apply(x$locusScores, 1, sd) * scalar
  if(scale.by[1] == 'numTaxa') {
	dotSizes.best <- sapply(1:length(X), function(z) mean(x$locus.total[names(which(unlist(lapply(bestTreeList, function(y, locNum) {locNum %in% y}, locNum = z)) == T))], na.rm = TRUE)) 
	dotSizes.worst <- sapply(1:length(X), function(z) mean(x$locus.total[names(which(unlist(lapply(worstTreeList, function(y, locNum) {locNum %in% y}, locNum = z)) == T))], na.rm = TRUE)) 
	# dotSizes.best[is.na(dotSizes.best)] <- 1
	}
  if(class(scale.by) == 'numeric') dotSizes.best <- dotSizes.worst <- scale.by
  layout(matrix(1:3, 1, 3))
  plot(X, Y.best, cex = dotSizes.best, xlab = 'Tree log-likelihood', ylab = paste('Number of loci for which tree is above the',percentile[2],'quantile'), ylim = range(c(Y.best, Y.worst)), ...)
  plot(X, Y.worst, cex = dotSizes.worst, xlab = 'Tree log-likelihood', ylab = paste('Number of loci for which tree is below the',percentile[1], 'quantile'), ylim = range(c(Y.best, Y.worst)), ...) 
  plot(X, Y.best, cex = dotSizes.best * (bw.scalar / scalar), xlab = 'Tree log-likelihood', ylab = 'Number of loci for which tree is above quantile (B) or below quantile (W)', ylim = range(c(Y.best, Y.worst)), pch = 'B', ...)
  points(X, Y.worst, cex = dotSizes.worst * (bw.scalar / scalar), col = 'red', pch = 'W', ...)
  segments(X, Y.best, X, Y.worst, lty = 'dashed')
  out <- cbind(Y.best, Y.worst, lnL = x$treeScores, difference = Y.best-Y.worst)
  row.names(out) <- names(X)
  windows()
  plot(out[, c('lnL', 'difference')], pch = 21, col = 'grey', bg = 'black', cex = 2)
  abline(h = quantile(out[, c('difference')], c(0.025, 0.975)), lty = 'dashed')
  v.bad <- which(out[, c('difference')] < quantile(out[, c('difference')], c(0.025, 0.975))[1])
  v.good <- which(out[, c('difference')] > quantile(out[, c('difference')], c(0.025, 0.975))[2])
  text(out[v.good, c('lnL', 'difference')], labels = names(v.good), pos = 2)
  text(out[v.bad, c('lnL', 'difference')], labels = names(v.bad), pos = 4)
  out <- list(lnL.diff.mat = out, v.bad = v.bad, v.good = v.good, n.taxa.best = dotSizes.best, n.taxa.worst = dotSizes.worst)
  return(out)
  } 

plot.besties <- function(x, tr = tree.rooted.di.noDupes.noOG.nni, x.text = 0, y.text = 20.5, fileBase = "quantileLoci", includeDate = T, ...) {
## takes plot.swulLikelihood as input
  if(includeDate) fileBase = paste(fileBase, format(Sys.time(), ".%Y-%m-%d"), sep = '')

  plot6 <- function(trees, filename) {
    pdf(filename, width = 11, height = 8.5)
	layout(matrix(1:6, 2, 3))
    for(i in trees) {
      plot(ladderize(tr[[i]]))
	  text(x.text, y.text, paste("Tree ", i, ", lnL = ", round(x$lnL.diff.mat[i, "lnL"], 1), '; loci for:against = ', x$lnL.diff.mat[i, "Y.best"], ":", x$lnL.diff.mat[i, "Y.worst"],  sep = ''), pos = 4, ...)
	  }
	dev.off()
	}
  
  plot6(x$v.good, paste(fileBase, ".bestRatio.trees.out.pdf", sep = ''))
  plot6(x$v.good, paste(fileBase, ".worstRatio.trees.out.pdf", sep = ''))
  plot6(head(order(x$lnL.diff.mat[, 'lnL']), 6), paste(fileBase, ".worstScore.trees.out.pdf", sep = ''))
  plot6(tail(order(x$lnL.diff.mat[, 'lnL']), 6), paste(fileBase, ".bestScore.trees.out.pdf", sep = ''))
  plot6(tail(order(x$lnL.diff.mat[, 'Y.best']), 6), paste(fileBase, ".mostSupportingLoci.trees.out.pdf", sep = ''))
  plot6(tail(order(x$lnL.diff.mat[, 'Y.worst']), 6), paste(fileBase, ".mostRejectingLoci.trees.out.pdf", sep = ''))
 
}
  
getLikelihoods.raxml <- function(dat, lnL = NA, missingSites = NA, which.loci, treeScoreFile = choose.files(multi = FALSE, caption = "Select RAxML site likelihoods file for trees")) {
  ## this version of the getLikelihood function tosses 
  ## ARGUMENTS:
  ##   treeScoreFile - filenames for site likelihoods from RAxML outfile and a vector of locus assignments for each site
  ##   dat - pyRAD.loci object
  ## VALUE: a list with (1) a matrix of locus likelihoods (columns) by trees (rows) and (2) a vector of tree likelihoods
  ## works: 2013-10-29
  
  ## 1. read data -- result is a list of site likelihoods, one per tree
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