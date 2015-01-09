compare.loci <-
function(taxa, radMat) {
  ## taxa is a pair of taxon vectors
  ## radMat is an inds.mat
  if(length(taxa[[1]]) > 1) v1 <- apply(radMat[taxa[[1]], ], 2, any)
  else v1 <- radMat[taxa[[1]], ]
  if(length(taxa[[2]]) > 1) v2 <- apply(radMat[taxa[[2]], ], 2, any)
  else v2 <- radMat[taxa[[2]], ]
  v.out <- v1 & v2
  out <- list(abs = sum(v.out), prop = sum(v.out) / sum(v1 | v2))
  return(out)
}

do.all.nodes <-
function(tr, rads) {
  to.do <- unique(tr$edge[, 1])
  out <- matrix(sapply(to.do, function(x) compare.loci(node.tips(tr, x), rads)), dimnames = list(to.do, 'sharedLoci'))
  out 
}

node.tips <-
function(tr, node) {
  to.do <- tr$edge[which(tr$edge[,1] == node), 2]
  out <- lapply(Descendants(tr, to.do, 'tips'), function(x) tr$tip.label[x])
  return(out)
}

do.a.rad.comparison.tree <- function(tr, rads, cex.scalar = 2, ...) {
  node.comparisons <- do.all.nodes(tr, rads)
  plot(tr)
  nodelabels(pie = nodes.comparisons$prop[,1], node = as.numeric(row.names(nodes.comparison$prop)), cex = cex.scalar * nodes.comparison$abs[, 1]/max(nodes.comparison$abs[, 1]), piecol = c('black', 'white'), ...)
  }