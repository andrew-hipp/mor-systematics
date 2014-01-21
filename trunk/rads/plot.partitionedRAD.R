plot.partitionedRAD <- function(radMat, scalar = 1, criterion = c('best.worst', 'percentile'),
                                cutoff.lnL = 1.5, minTrees = 10, 
								opt = c('diff', 'both'), add.opt = T, ...) {
## each tree is a data point
## X: tree likelihood
## Y: 
  X <- colSums(radMat)
  X <- x$treeScores
  if(minTrees > 1) {
    nTrees <- apply(radMat, 1, function(x) length(unique(x)))
	radMat <- radMat[nTrees >= minTrees, ]
	}
  loc.scores <- x$locusScores[, diff(x) >= cutoff.lnL]
  bestTreeList <- apply(loc.scores, 2, function(z) which(z > quantile(z, percentile[2])))
  bestTree <- unlist(bestTreeList)
  worstTreeList <- apply(loc.scores, 2, function(z) which(z < quantile(z, percentile[1])))
  worstTree <- unlist(worstTreeList)
  Y.best <- sapply(1:length(X), function(z) sum(bestTree == z))
  Y.worst <- sapply(1:length(X), function(z) sum(worstTree == z))
  if(scale.by[1] == 'sd') dotSizes.best <- dotSizes.worst <- apply(loc.scores, 1, sd) * scalar
  if(scale.by[1] == 'max-min') {
	# THIS IS NOT RIGHT
	warning('max-min dot sizes is not currently implemented correctly! do not put any stock in it at all!')
	dotSizes.best <- abs(sapply(1:length(X), function(z) sum(c(1, -1) * range(loc.scores[, z])))) * scalar
	dotSizes.worst <- abs(sapply(1:length(X), function(z) sum(c(1, -1) * range(loc.scores[, z])))) * scalar
	}
  if(scale.by[1] == 'numTaxa') {
	dotSizes.best <- sapply(1:length(X), function(z) mean(x$locus.total[names(which(unlist(lapply(bestTreeList, function(y, locNum) {locNum %in% y}, locNum = z)) == T))], na.rm = TRUE)) * scalar
	dotSizes.worst <- sapply(1:length(X), function(z) mean(x$locus.total[names(which(unlist(lapply(worstTreeList, function(y, locNum) {locNum %in% y}, locNum = z)) == T))], na.rm = TRUE)) * scalar
	# dotSizes.best[is.na(dotSizes.best)] <- 1
	}
  out <- cbind(Y.best, Y.worst, lnL = x$treeScores, difference = Y.best-Y.worst)
  row.names(out) <- names(X)
  if(class(scale.by) == 'numeric') dotSizes.best <- dotSizes.worst <- scale.by
  layout(matrix(1:3, 1, 3))
  plot(X, Y.best, cex = dotSizes.best, xlab = 'Tree log-likelihood', ylab = paste('Number of loci for which tree is above the',percentile[2],'quantile'), ylim = range(c(Y.best, Y.worst)), pch = 21, bg = 'black', col = 'gray')
  if(add.opt) points(X[1], Y.best[1], cex = dotSizes.best, pch = 21, col = 'gray', bg = 'yellow')
  plot(X, Y.worst, cex = dotSizes.worst, xlab = 'Tree log-likelihood', ylab = paste('Number of loci for which tree is below the',percentile[1], 'quantile'), ylim = range(c(Y.best, Y.worst)), pch = 21, bg = 'red', col = 'black') 
  if(add.opt) points(X[1], Y.worst[1], cex = dotSizes.worst, pch = 21, col = 'black', bg = 'yellow')
  if(opt[1] == 'both') {
    plot(X, Y.best, cex = dotSizes.best * (bw.scalar / scalar), xlab = 'Tree log-likelihood', ylab = 'Number of loci for which tree is above quantile (B) or below quantile (W)', ylim = range(c(Y.best, Y.worst)), pch = 21, bg = 'black', col = 'gray')
    points(X, Y.worst, cex = dotSizes.worst * (bw.scalar / scalar), bg = 'red', pch = 21, col = 'black')
    segments(X, Y.best, X, Y.worst, lty = 'dashed')
	}
  else {
    plot(out[, c('lnL', 'difference')], pch = 21, col = 'grey', bg = 'black', cex = 2, xlab='Tree log-likelihood', ylab = 'Number of loci favoring minus number disfavoring tree')
    points(out[1, 'lnL'], out[1, 'difference'], pch = 21, col = 'grey', bg = 'yellow', cex = 2)
	abline(h = quantile(out[, c('difference')], c(0.025, 0.975)), lty = 'dashed')
    v.bad <- which(out[, c('difference')] < quantile(out[, c('difference')], c(0.025, 0.975))[1])
    v.good <- which(out[, c('difference')] > quantile(out[, c('difference')], c(0.025, 0.975))[2])
    text(out[v.good, c('lnL', 'difference')], labels = names(v.good), pos = 2)
    text(out[v.bad, c('lnL', 'difference')], labels = names(v.bad), pos = 4)
    out <- list(lnL.diff.mat = out, v.bad = v.bad, v.good = v.good, dot.sizes.best = dotSizes.best, dot.sizes.worst = dotSizes.worst)
	}
  return(out)
  } 
