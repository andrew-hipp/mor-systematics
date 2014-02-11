plot.rankedPartitionedRAD <-
function(x, 
         tree.lnL.file = NULL,
		 fileprefix = NULL,
		 lnL.break = NULL,
		 regression = NULL,
         widthScalar = .85,
         panels = c('bestMat', 'worstMat', 'doubleCountMat'),
         squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3),
         primeTreeColor = 'red',
         primeTreeCharacter = 19,
         filebase = 'DEFAULT',
         ...) {
  if(filebase == 'DEFAULT') filebase <- paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."),  paste(c('minT','rangeL','diffL','noDoubles'), x$params, collapse = "_", sep = ''), '.pdf', sep = '')
  if(class(x) != 'rankedPartitionedRAD') warning('Not the expected object class; this function may misbehave')
  if(!is.null(lnL.break)) {
    break.out <- try(names(lnL.break) <- panels, silent = TRUE)
	if(class(break.out) == 'try-error') {
	  warning('lnL.break not equal in length to panels, so ignored')
	  lnL.break <- NULL
	  }
	}
  if(is.null(regression)) regression <- rep(FALSE, length(panels))
  break.out <- try(names(regression) <- panels, silent = TRUE)
  if(class(break.out) == 'try-error') {
	warning('regression not equal in length to panels, so ignored')
	regression <- NULL
	}
 if(is.null(tree.lnL.file)) {
    trees.lnL <- colSums(x$radMat)
	temp.xlab <- 'Tree log-likelihood, summed over loci'
	}
  else {
    trees.lnL <- get.raxml.treeLikelihoods(tree.lnL.file)
	temp.xlab <- 'Tree log-likelihood, full data matrix'
	}
  if(!is.null(fileprefix)) pdf(paste(fileprefix, filebase, sep = '.'), width = squareSize*length(panels)*widthScalar, height = squareSize)
  layout(matrix(seq(length(panels)), 1, length(panels)))
  panelCount <- 0
  for(i in panels) {
    panelCount <- panelCount + 1
	temp.ylab = switch(i, 
	                   bestMat = 'Number of loci supporting tree',
					   worstMat = 'Number of loci disfavoring tree',
					   doubleCountMat = 'Loci overlapping, excluded'
					   )
	mat.lnL <- colSums(x[[i]])
	if(!is.null(lnL.break)) {
	  trees.x <- trees.lnL[trees.lnL > lnL.break[panelCount]]
	  mat.y <- mat.lnL[trees.lnL > lnL.break[panelCount]]
	  }
	plot(trees.x, mat.y, xlab = temp.xlab, ylab = temp.ylab, type = 'n', ...)
	points(trees.x[-c(1)], mat.y[-c(1)], ...)
    points(trees.x[1], mat.y[1], pch = primeTreeCharacter, col = primeTreeColor)
	if(regression[i]) abline(lm(mat.y ~ trees.x))
	}
  if(!is.null(fileprefix)) dev.off()
  return('done')
  }
