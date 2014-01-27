plot.rankedPartitionedRAD <- function(x, 
                                      fileprefix = NULL,
									  widthScalar = .8,
									  primeTreeColor = 'red',
									  primeTreeCharacter = 19,
									  filebase = paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."),  paste(c('minT','rangeL','diffL','noDoubles'), x$params, collapse = "_", sep = ''), '.pdf', sep = ''),
                                      panels = c('bestMat', 'worstMat', 'doubleCountMat'),
									  ...) {
  if(class(x) != 'rankedPartitionedRAD') warning('Not the expected object class; this function may misbehave')
  if(!is.null(fileprefix)) pdf(paste(fileprefix, filebase, sep = '.'), width = 3*length(panels)*widthScalar, height = 3)
  layout(matrix(seq(length(panels)), 1, length(panels)))
  for(i in panels) {
    temp.ylab = switch(i, 
	                   bestMat = 'Number of loci supporting tree',
					   worstMat = 'Number of loci disfavoring tree',
					   doubleCountMat = 'Loci overlapping, excluded'
					   )
	plot(x$trees.lnL, colSums(x[[i]]), xlab = 'Tree log-likelihood, summed over loci', ylab = temp.ylab, ...)
    points(x$trees.lnL[1], colSums(x[[i]])[1], pch = primeTreeCharacter, col = primeTreeColor)
	}
  dev.off()
  return('done')
  }

