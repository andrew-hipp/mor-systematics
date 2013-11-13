clusterSites <- function(locusVector) {
## 2013-10-23: this should be obsolete now, as read.pyRAD generates a matrix that can be used to infer this info
## left in in case it becomes handy to anyone
  clusterEnds <- cumsum(locusVector) # a vector of the position for each cluster first bp
  clusterBegs <- clusterEnds - locusVector + 1 # a vector of the position for each cluster last bp
  out <- lapply(1:length(locusVector), function(x) seq(from = clusterBegs[x], to = clusterEnds[x])) # list of base pairs in each cluster
  return(out)
  }

