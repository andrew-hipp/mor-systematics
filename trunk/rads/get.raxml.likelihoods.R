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
  class(out.mat) <- 'partitionedRAD'
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