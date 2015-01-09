compare.loci <-
function(taxa, radMat, give = c('abs', 'prop')) {
  ## taxa is a pair of taxon vectors
  ## radMat is an inds.mat
  if(length(taxa[[1]]) > 1) v1 <- apply(radMat[taxa[[1]], ], 2, any)
  else v1 <- radMat[taxa[[1]], ]
  if(length(taxa[[2]]) > 1) v2 <- apply(radMat[taxa[[2]], ], 2, any)
  else v2 <- radMat[taxa[[2]], ]
  v.out <- v1 & v2
  if(give[1] == 'abs') out <- sum(v.out)
  if(give[1] == 'prop') out <- sum(v.out) / sum(v1 | v2)
  return(out)
}
do.all.nodes <-
function(tr, rads, ...) {
  to.do <- unique(tr$edge[, 1])
  out <- matrix(sapply(to.do, function(x) compare.loci(node.tips(tr, x), rads, ...)), dimnames = list(to.do, 'sharedLoci'))
  out 
}
node.tips <-
function(tr, node) {
  to.do <- tr$edge[which(tr$edge[,1] == node), 2]
  out <- lapply(Descendants(tr, to.do, 'tips'), function(x) tr$tip.label[x])
  return(out)
}
