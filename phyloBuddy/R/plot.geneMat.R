plot.geneMat <- function(x, tr = NA, genes = dimnames(x)[[2]][6:dim(x)[2]], sortByFreq = T, ...) {
  x <- t(x) # puts genes as rows, inds as columns
  x <- x[genes, ]
  if(sortByFreq) x <- x[names(sort(rowSums(x != ''), decreasing = T)), ]
  if(!is.na(tr[1])) {
    tr <- read.tree(text = write.tree(tr))
	inds <- tr$tip.label
	x <- x[, inds]
	layout(matrix(1:2, 2,1))
	} ## should add an option to subset by individuals
  nInds <- dim(x)[2]
  nGenes <- dim(x)[1]
  plot(1, xlim = c(1, nInds), ylim = c(1, nGenes), type = 'n')
  for(y in seq(nGenes)) {
    for(x in seq(nInds)) {
      if(x[y, x] != "") points(x,y, ...)
      }
	}
  if(!is.na(tr[1])) plot(tr, direction = 'upwards', show.tip.label = F)
  }