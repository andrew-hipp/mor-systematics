rank.partitionedRAD <- function(radMat, criterion = c('lnL.threshold', 'percentile'),
                                minTrees = 10, min.overall.diff.lnL = 1.5, discardDoublecounts = TRUE,
								threshold.lnL = 0) {

  X <- colSums(radMat)
  X <- x$treeScores
  # bestMat <- worstMat <- matrix(0, nrow = dim(radMat)[1], ncol = dim(radMat)[2], dimnames = dimnames(radMat))
  if(minTrees > 1) {
    nTrees <- apply(radMat, 1, function(x) length(unique(x)))
	radMat <- radMat[nTrees >= minTrees, ]
	}
  if(min.overall.diff.lnL > 0) {
    lnL.diff <- abs(diff(t(apply(radMat, 1, range))))
	radMat <- radMat[lnL.diff >= min.overall.diff.lnL, ]
	}
  if(criterion[1] == 'lnL.threshold') {
    bestMat <- t(apply(radMat, 1, function(x) abs(x - max(x)) <= threshold.lnL))
	worstMat <- t(apply(radMat, 1, function(x) abs(x - min(x)) <= threshold.lnL))
	doubleCountMat <- bestMat & worstMat
	}
  ## other options here
  if(discardDoubleCounts) {
    doubleCounts <- apply(doubleCountMat, 1, sum) > 0
	bestMat <- bestMat[!doubleCounts, ]
	worstMat <- worstMat[!doubleCounts, ]
	}
  out <- list(bestMat = bestMat, worstMat = worstMat, doubleCounts = doubleCounts)
  out
  }