sort.strObj <- function(x) {
  x.clusters <- grep('cluster', names(x), value = T)
  for(i in x.clusters) x <- x[order(x[[i]]), ]
  return(x)
}
