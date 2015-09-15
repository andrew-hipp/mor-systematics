plot.geneMat <- function(x, tr = NA, genes = colnames(x)[6:(dim(x)[2])], panes = c(1,3), margins = c(0,1,0.5,0), geneColors = c('red', 'black'), sortByFreq = T, label.cex = 0.5, ...) {
  x <- t(x) # puts genes as rows, inds as columns
  x <- x[genes, ]
  if(sortByFreq) x <- x[names(sort(rowSums(x != ''), decreasing = T)), ]
  if(!is.na(tr[1])) {
    tr <- read.tree(text = write.tree(tr))
	inds <- tr$tip.label
	x <- x[, inds]
	layout(matrix(1:2, 2,1), heights = panes)
	} ## should add an option to subset by individuals
  nInds <- dim(x)[2]
  nGenes <- dim(x)[1]
  col.genes <- rep(geneColors, nGenes/length(geneColors))
  par(mai=margins)
  plot(1, xlim = c(1, nInds), ylim = c(1, nGenes), type = 'n', xlab = '', ylab = '', axes = F)
  axis(2, at = seq(from = 1, to = nGenes-1, by = 2), labels = row.names(x)[seq(from = 1, to = nGenes-1, by = 2)], las = 2, cex.axis = label.cex, col.axis = col.genes[1])
  axis(2, at = seq(from = 2, to = nGenes, by = 2), labels = row.names(x)[seq(from = 2, to = nGenes, by = 2)], las = 2, cex.axis = label.cex, col.axis = col.genes[2])
  
  for(i in seq(nGenes)) {
    for(j in seq(nInds)) {
      if(x[i, j] != "") points(j,i, col = col.genes[i], ...)
      }
	}
  if(!is.na(tr[1])) {
    plot(tr, direction = 'upwards', show.tip.label = F, no.margin = F)
    }
  }