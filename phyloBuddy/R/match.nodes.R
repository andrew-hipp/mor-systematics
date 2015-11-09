match.nodes <- function(guideTree, otherTrees, matchBoots = T, plotBoots = T, ...) {
require(phytools)
  if('phylo' %in% class(otherTrees)) otherTrees = list(tr1 = otherTrees)
  matched <- lapply(otherTrees, function(x) matchNodes(guideTree, x))
  match.mat <- do.call(cbind, c(guideTree = list(matched[[1]][, 1]), lapply(matched, function(x) x[, 2])))
  match.index <- which(apply(match.mat, 1, function(x) !any(is.na(x))))
  
  ## match up bootstraps
  if(matchBoots) {
    allTreesLabels <- lapply(c(guideTree = list(guideTree), otherTrees), function(x) c(x$tip.label, x$node.label))
    mat.boots <- match.mat[match.index, ]
    for(i in seq(dim(mat.boots)[2])) mat.boots[, i] <- as.integer(allTreesLabels[[i]][mat.boots[, i]])
    out <- list(mat.boots = mat.boots, mat.full = match.mat, mat.index = match.index)
    if(plotBoots) matplot(mat.boots[order(its.ets.matched$mat.boots[, 1], decreasing = T), ], ...)
    }
  else out <- list(mat.full = match.mat, mat.index = match.index)
  return(out)
  }

## EXAMPLE
# pdf('its.ets.bootComparison.2015-11-09.pdf', 5, 5)
# matplot(its.ets.matched$mat.boots[order(its.ets.matched$mat.boots[, 1], decreasing = T), ], type = 'l', lwd = c(2,0.5,0.5), lty = c('solid'), col = c('black', 'red', 'blue'), xlab = 'node (arbitrary)', ylab='Bootstrap')
# legend(x=0, y=20, legend=c('ITS + ETS (concatenated matrix)', 'ITS only', 'ETS only'), bty = 'n', lwd = c(2), lty = c('solid'), col = c('black', 'red', 'blue'), cex = 0.8)
# dev.off()
