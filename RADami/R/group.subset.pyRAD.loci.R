group.subset.pyRAD.loci <- function(dat, groups, mins = 10, loci = names(dat$DNA), use.tidyName = TRUE, ...) {
# arguments:
#  dat is a subset.pyRAD.loci object
#  groups is a list of individuals (can be named) in the groups
#  
# value: a matrix with the number of individuals in each group for each locus

  if(use.tidyName) {
    names(dat$DNA) <- tidyName(names(dat$DNA), ...)
	groups <- lapply(groups, tidyName, ...)
	}
  out <- mcmapply(dat$DNA, function(x) c(sapply(groups, function(y) sum(names(x) %in% y)), total = length(x)), mc.cores = 4)
  if(!is.na(mins[1])) leave.in <- apply(out, 1, function(x) all(x >= mins))
  else leave.in <- rep(TRUE, dim(out)[1])
  if(!is.na(loci[1])) out <- out[intersect(row.names(out), loci), ]
  out
  }
 