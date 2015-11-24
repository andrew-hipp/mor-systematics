plot.geneMat <- function(x, tr = NA, genes = colnames(x)[6:(dim(x)[2])], panes = c(5,1,2), margins = c(0,1,0.5,0), geneColors = c('red', 'black'), minGenes = 10, sortByFreq = T, label.cex = 0.7, pdfTitle = format(Sys.time(), 'geneMat.plot.%Y.%m.%d.%H%M.pdf'), pdfW = 10, pdfH = 10, ...) {
  require(ape)
  if(!is.na(pdfTitle)) pdf(pdfTitle, pdfW, pdfH)
  x <- t(x) # puts genes as rows, inds as columns
  x <- x[genes, ]
  if(sortByFreq) x <- x[names(sort(rowSums(x != ''), decreasing = T)), ]
  if(!is.na(tr[1])) {
    tr <- read.tree(text = write.tree(tr))
	inds <- tr$tip.label
	x <- x[, inds]
	layout(matrix(1:3, 3,1), heights = panes)
	} ## should add an option to subset by individuals
  x <- x[apply(x, 1, function(y) sum(y != '')) >= minGenes, ]
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
    geneSums <- apply(x, 2, function(y) sum(y != ''))
	plot(geneSums, type = 'p', pch = 19, axes = F, ylab = 'Genes sampled', cex.lab = 0.8)
	abline(h = mean(geneSums), lty = 'dashed')
	axis(2, at = seq(geneSums), labels = seq(geneSums), las = 2)
	plot(tr, direction = 'upwards', show.tip.label = F, no.margin = F)
    }
  if(!is.na(pdfTitle)) dev.off()
  invisible(0) # should return something more useful here
  }
