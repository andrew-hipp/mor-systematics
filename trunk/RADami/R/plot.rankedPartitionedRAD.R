plot.rankedPartitionedRAD <-
function(x, 
         tree.lnL.file = NULL,
		 fileprefix = NULL,
         widthScalar = .85,
         panels = c('bestMat', 'worstMat', 'doubleCountMat'),
         squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3),
         primeTreeColor = 'red',
         primeTreeCharacter = 19,
         filebase = 'DEFAULT',
         ...) {
  if(filebase == 'DEFAULT') filebase <- paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."),  paste(c('minT','rangeL','diffL','noDoubles'), x$params, collapse = "_", sep = ''), '.pdf', sep = '')
  if(class(x) != 'rankedPartitionedRAD') warning('Not the expected object class; this function may misbehave')
  if(is.null(tree.lnL.file)) {
    trees.lnL <- colSums(x$radMat)
	temp.xlab <- 'Tree log-likelihood, full data matrix'
	}
  else {
    trees.lnL <- get.raxml.treeLikelihoods(tree.lnL.file)
	temp.xlab <- 'Tree log-likelihood, summed over loci'
	}
  if(!is.null(fileprefix)) pdf(paste(fileprefix, filebase, sep = '.'), width = squareSize*length(panels)*widthScalar, height = squareSize)
  layout(matrix(seq(length(panels)), 1, length(panels)))
  for(i in panels) {
    temp.ylab = switch(i, 
	                   bestMat = 'Number of loci supporting tree',
					   worstMat = 'Number of loci disfavoring tree',
					   doubleCountMat = 'Loci overlapping, excluded'
					   )
	plot(trees.lnL, colSums(x[[i]]), xlab = temp.xlab, ylab = temp.ylab, ...)
    points(trees.lnL[1], colSums(x[[i]])[1], pch = primeTreeCharacter, col = primeTreeColor)
	}
  dev.off()
  return('done')
  }
