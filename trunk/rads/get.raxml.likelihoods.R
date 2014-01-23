match.lnL.to.trees <- function(directory = 'getwd()', lnLprefix = 'RAxML_info.', lnLsuffix = '.lnL', treeIndexFile = 'tree.index.lines.txt', locus.names = NULL, ...) {
## updated 2014-01-23 to prune out files where no trees written -- TO TEST 
  treeIndex <- read.delim(paste(directory, '/', treeIndexFile, sep = ''), as.is = T, header = F, row.names = 1)
  if(is.null(locus.names)) locus.names <- row.names(treeIndex)
  logfile = file(format(Sys.time(), "match.lnL.%Y-%m-%d.log.txt"), open = "a")
  lnL.list <- lapply(paste(directory, '/', lnLprefix, locus.names, lnLsuffix, sep = ''), get.raxml.treeLikelihoods)
  close(logfile)
  names(lnL.list) <- locus.names
  raxml.worked <- names(which(lnL.list != 'FAIL')) ## added 2014-01-23
  lnL.list <- lnL.list[raxml.worked] ## added 2014-01-23
  treeIndex <- treeIndex[raxml.worked, ] ## added 2014-01-23
  locus.names <- raxml.worked ## added 2014-01-23
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

get.raxml.treeLikelihoods <- function(x, logfile = NA) {
## gets likelihoods from the RAxML_info file
## updated 2014-01-23 to deliver a fail if no trees written
	# cat(paste('working on file', x, '\n'), file = logfile)
	fileIn <- readLines(x)
	if(length(grep("Tree ", fileIn)) == 0) {
	  if(!is.na(logfile)) cat('... file', x, 'had no trees in it.', '\n')
	  return('FAIL')
	  }
	out <- as.double(sapply(grep("Tree ", fileIn, value=T), function(x) strsplit(x, ": ")[[1]][2]))
	names(out) <- as.character(1:length(out))
	out
	}