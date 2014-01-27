plot.rankedPartitionedRAD <- function(x, 
                                      fileprefix = NULL,
									  widthScalar = .85,
									  panels = c('bestMat', 'worstMat', 'doubleCountMat'),
									  squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3),
									  primeTreeColor = 'red',
									  primeTreeCharacter = 19,
									  filebase = paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."),  paste(c('minT','rangeL','diffL','noDoubles'), x$params, collapse = "_", sep = ''), '.pdf', sep = ''),
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
	plot(colSums(x$radMat), colSums(x[[i]]), xlab = 'Tree log-likelihood, summed over loci', ylab = temp.ylab, ...)
    points(colSums(x$radMat)[1], colSums(x[[i]])[1], pch = primeTreeCharacter, col = primeTreeColor)
	}
  dev.off()
  return('done')
  }

