match.lnL.to.trees <- function(locus.names, lnLprefix = 'RAxML_info.', lnLsuffix = '.lnL', treeIndexFile = 'tree.index.lines.txt', directory = 'getwd()', ...) {
  treeIndex <- read.delim(treeIndexFile, as.is = T, header = F, row.names = 1)
  lnL.list <- lapply(paste(lnLprefix, locus.names, lnLsuffix, sep = ''), get.raxml.treeLikelihoods)
  names(lnL.list) <- locus.names
  out.mat <- matrix(NA, nrow = length(locus.names), ncol = dim(treeIndex)[1], dimnames = list(locus.names, NULL))
  
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
	out <- as.double(sapply(grep("Tree ", a, value=T), function(x) strsplit(x, ": ")[[1]][2]))
	names(out) <- as.character(1:length(out))
	out
	}
