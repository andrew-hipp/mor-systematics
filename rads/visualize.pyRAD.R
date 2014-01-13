locus.dist <- function(pyIn, proportional = TRUE, upper = TRUE, diagonal = TRUE) {
  if(class(pyIn) == 'pyRAD.loci') pyIn <- pyIn$radSummary$inds.mat
  numInds <- dim(pyIn)[1]
  numLoci <- ifelse(proportional, dim(pyIn)[2], 1)
  out <- matrix(NA, numInds, numInds, dimnames = list(row.names(pyIn), row.names(pyIn)))
  for(i in 1:numInds) {
    for(j in 1:i) {
	  out[i, j] <- sum(colSums(pyIn[c(i, j), ]) == 2) / numLoci
	  }}
  if(upper) out <- as.matrix(as.dist(out))
  if(diagonal) diag(out) <- apply(pyIn, 1, sum) / numLoci
  out
  }

plot.locus.dist <- function(locD, tr, trW = 3, plotW = 5, labelsW = 3, plotGap = 0.25, scalar = 1.5, barH = 1, point.pch = c(21,25), tips = NULL, diag.col = 'red', ...) {
 require(geiger)
 require(ape)
 if(barH > 0) layout(matrix(c(2,1,3),3,1), widths = rep(trW + plotW + labelsW + plotGap * 2, 3), heights = c(barH, plotW, trW), TRUE)
 orig.tips <- tr$tip.label
 names(orig.tips) <- tips
 if(!is.null(tips)) tr$tip.label <- tips
 tr <- read.tree(text = write.tree(tr)) # reorder tip labels to match topology
 orig.tips <- orig.tips[gsub("_", " ", tr$tip.label, fixed = T)] # reorder original tips according to new tree order
 locD <- locD[orig.tips, orig.tips]
 color.mat <- matrix('black', dim(locD)[1], dim(locD)[2])
 if(!is.null(diag.col)) diag(color.mat) <- diag.col
 nloci <- dim(locD)[1]
 plot(transform(tr, 'depth', trW), x.lim = c(0, trW + plotW + plotGap * 2 + labelsW), show.tip.label=F, no.margin = T)
 xy <- matrix(seq(nloci), nloci, nloci, byrow = TRUE)
 Xs <- plotW*(as.numeric(t(xy)) / nloci) + trW + plotGap
 points(Xs, as.numeric(xy), pch = point.pch[1], cex = as.numeric(locD) * scalar, bg = as.character(color.mat), col = 'black')
 if(labelsW > 0) text(rep((trW + plotW + plotGap * 2), length(tr$tip.label)), 1:length(tr$tip.label), gsub("_", " ", tr$tip.label, fixed = T), pos = 4, ...)
 if(barH > 0) {
   Xs.1 = Xs[1:length(tr$tip.label)]
   heights = apply(locD, 2, function(x) sum(x) / (dim(locD)[1]-1)) * barH 
   plot(1, xlim = c(0, trW + plotW + labelsW + plotGap * 2), ylim = c(0, barH * 1.1), type = 'n', axes = F) # just to move to next layout area
   segments(x0 = Xs.1, y0 = 0, x1 = Xs.1, y1 = heights, lwd = 10, lty = 1, lend = 2)
   # axis(side = 4, pos = tail(Xs, 1) + plotGap, las = 2, cex.axis = 0.6) # right side
   # axis(side = 2, pos = Xs[1] - plotGap, las = 2, cex.axis = 0.6) # left side
   text(Xs.1, heights, round(heights, 2), pos = 3, cex = 0.5) # numbers above bars
   }
 tr.x.max <- (length(tr$tip.label) / plotW) * (plotW + labelsW + plotGap)
 tr.x.min <- -(length(tr$tip.label) / plotW) * (trW + plotGap)
 legend.header.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap)
 legend.dots.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap * 2)
 legend.text.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap * 3)
 legend.y <- c(trW, trW * 0.6, trW * 0.2)
 print(c(tr.x.min, tr.x.max))
 plot(transform(tr, 'depth', trW), x.lim = c(tr.x.min, tr.x.max), direction = 'upwards', show.tip.label=F, no.margin = T)
 points(rep(legend.dots.x, 2), legend.y[2:3], pch = point.pch, bg = 'black', col = 'black', cex = c(scalar, scalar / 2))
 text(c(legend.header.x, rep(legend.text.x, 2)), legend.y, c('Proportion of total loci', '1.0', '0.5'), cex = c(0.8, 0.6, 0.6), pos = 4)
 out = invisible(list(Xs = Xs[1:length(tr$tip.label)], heights = apply(locD, 2, mean) * barH))
 return(out)
 }
  
overlap.report <- function(dat, repPattern = "_re", origPattern = "_h") {
## reports on how replicate individuals fall out in the loci
  if(class(dat) == 'summary.pyRAD.loci') dat <- dat$inds.mat
  reps <- grep(repPattern, row.names(dat), fixed = TRUE, value = TRUE)
  orig <- gsub(repPattern, origPattern, reps, fixed = TRUE)
  outCols <- c("Original", "Replicate", "Intersection", "Union", "Not intersection", "Overlap proportion", "Original loci replicated", "Increase, original to replicate")
  out <- matrix(NA, nrow = length(orig), ncol = length(outCols), dimnames = list(orig, outCols))
  for(i in 1:length(orig)) {
    message(paste("Doing names", orig[i]))
	out[orig[i], "Original"] <- sum(dat[orig[i],])
	out[orig[i], "Replicate"] <- sum(dat[reps[i],])
	out[orig[i], "Intersection"] <- sum(colSums(dat[c(orig[i], reps[i]), ]) == 2, na.rm = T)
	out[orig[i], "Union"] <- sum(colSums(dat[c(orig[i], reps[i]), ]) %in% 1:2, na.rm = T)
	}
  out[, "Not intersection"] <- out[, "Union"] - out[, "Intersection"]
  out[, "Overlap proportion"] <- round(out[, "Intersection"] / out[, "Union"], 3)
  out[, "Original loci replicated"] <- round(out[, "Intersection"] / out[, "Original"], 3)
  out[, "Increase, original to replicate"] <- round(out[, "Replicate"] / out[, "Original"], 3)
  return(out)
  }

lengths.report <- function(dat, numtodo = 10, reportInterval = 2000, high.mem = TRUE) {
## set numtodo to 0 if you want to do all loci
  if(class(dat) != 'pyRAD.loci') stop("This function runs on a pyRAD data object")
  last.lines <- dat$cons - 1
  num.loci <- length(last.lines)
  datSeqs <- as.character(dat$seqs)
  if(high.mem) block.lengths <- sapply(datSeqs[last.lines][1:ifelse(numtodo<1,num.loci,numtodo)], function(x) nchar(as.character(x)))
  else {
    block.lengths = integer(0)
	for(i in 1:num.loci) block.lengths = c(block.lengths, nchar(as.character(datSeqs[last.lines[i]])))
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', num.loci, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (num.loci - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }

	}
  names(block.lengths) <- dat$locus.index[last.lines]
  return(block.lengths)
  }

