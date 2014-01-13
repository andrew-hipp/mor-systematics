## handling pyRAD data
## A Hipp, 2012-09-18
## -- Sept 2012: accommodate new pyRAD format (GBS data)
## -- Dec 2012: new functions to help with exporting RAD data and blasting
## -- Jan 2013: consensus functions now use Biostrings instead of seqinr

## uncomment if not already installed:
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

rad2phy <- function(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], outfile = 'pyMat.out.phy', padding = 50) {
## makes a phylip-style data matrix from rad.mat output, limiting by individuals and loci
  if(class(pyDat) != "rad.mat") warning("I'm expecting output from rad.mat")
  outfile = file(outfile, "wt")
  open(outfile)
  cat(paste(length(inds), sum(sapply(pyDat[inds[1], loci], nchar)), "\n"), file = outfile) #header: number of individuals, number of bases
  for(i in inds) {
    message(paste("Writing DNA line for individual", i))
	cat(i, file = outfile)
	cat(paste(rep(" ", padding - nchar(i)), collapse = ""), file = outfile)
	cat(paste(pyDat[i, loci], collapse = ""), file = outfile)
	cat("\n", file = outfile) # endline
	}
  close(outfile)
  }

locus.picker <- function(pyDat, minThreshold = 3, inds = row.names(pyDat$radSummary$inds.mat)) {
  out = names(which(colSums(pyDat$radSummary$inds.mat[inds, ]) >= minThreshold))
  }
  
rad2mat <- function(pyDat, fill.N = TRUE) {
## shoves RAD sequence into an individuals x loci matrix
  if(!"radSummary" %in% names(pyDat)) pyDat$radSummary <- summary(pyDat) # calls summary if it wasn't done at read-time 
  loci <- dimnames(pyDat$radSummary$inds.mat)[[2]]
  inds <- dimnames(pyDat$radSummary$inds.mat)[[1]]
  out <- matrix(NA, length(inds), length(loci), dimnames = list(inds, loci))
  for(i in loci) {
    message(paste('Doing locus', i))
	out[, i] <- pyDat$radSummary$seqs.per.locus[[i]][inds]
    if(fill.N) out[inds[!inds %in% pyDat$radSummary$tips.per.locus[[i]]], i] <- paste(rep("N", pyDat$radSummary$locus.lengths[[i]]), collapse = "")
	}
  class(out) <- 'rad.mat'
  return(out)
}

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
  
consensus.pyRAD <- function(pyIn, from = NA, to = NA, fastaNames = T, writeFile = 'rads.con.1_100.txt', ...) {
## use seqinr to generate a consensus sequence for each pyRAD locus
## 2013-01-04: updated to use Biostrings, which works better -- deleted arguments: method = 'majority', threshold = 0.001
  if(class(pyIn) != "pyRAD.loci") stop("pyRAD input required, from read.pyRAD")
  # require(seqinr)
  require(Biostrings)
  allLoci <- unique(as.character(pyIn$locus.index))
  allLoci <- allLoci[!allLoci == ""] ## this should probably part of read.pyRAD
  if(!is.na(from)) allLoci <- allLoci[from:to]
  seqs <- as.character(pyIn$seqs)
  loc.index <- as.character(pyIn$locus.index)
  out <- character(0)
  for (i in allLoci) {
    message(i) #only for debugging
	out <- c(out, consensusString(DNAStringSet(gsub("-", "N", seqs[loc.index == i])), ...))
	}
  # for(i in allLoci) out <- c(out, paste(consensus(str2mat(seqs[loc.index == i]), method, threshold), collapse = ""))
  if(fastaNames) allLoci <- paste(">", allLoci, sep = "")
  names(out) <- allLoci
  if(!is.na(writeFile)) write.table(out, writeFile, sep = "\n", quote = F, col.names = F)
  return(out)
  }
  
blast.pyRAD <- function(pyConsensus, ...) {
## note written; placeholder for a function to help automate RAD blasts
}

read.pyRAD <- function(filename, reportInterval = 20000, breakLinesSeparate = FALSE, doSummary = TRUE, makeSeqMat = TRUE, ...) {
## reads the all.aligned file out of pyRAD, parses into names, loci, sequences
## updated with breakLinesSeparate in Oct 2012 because pyRAD switched to single-line summaries at the end of each aligned file
## updated 2012-11-16 to keep breaklines intact
## updated 2012-12-19 to get rid of all conversions to factors... apparently no longer needed for space considerations in R
  message("Reading data...")
  dat <- readLines(filename, ...)
  dat.breakLines <- dat.consensusLines <- grep("//", dat, fixed = TRUE) # this is slow, but only ca. 1 sec for data runs of 10s of thousands
  dat.breakLines.vector <- dat[dat.breakLines] #added 2012-11-16; ignores possibility of separate breakLines
  message("Splitting data...")
  dat.split <- strsplit(dat, " {1,100}") #uses whitespace to separate taxon names from sequences
  dat.names <- sapply(dat.split, function(x) x[1])
  dat.seqs <- sapply(dat.split, function(x) x[2])
  # dat.seqs[dat.breakLines] <- dat.breakLines.vector # shoves the consensus seqs back into the sequence vector, assuming breakLinesSeparate = F
  dat.firstLocusLines <- c(1, (dat.breakLines[1:(length(dat.breakLines)-1)] + 1))
  if(breakLinesSeparate) {
    dat.lastLocusLines <- dat.breakLines - 2
    dat.consensusLines <- dat.breakLines - 1
	}
  else dat.lastLocusLines <- dat.breakLines - 1
  dat.locus.index <- character(length(dat))
  locusCounter <- 1
  message("Assigning locus number...")
  start.time <- Sys.time()
  ## crazy slow! before vectorization:
  #for (i in 1:length(dat)) {
  #  if(i-1 %in% dat.breakLines) locusCounter <- locusCounter + 1
  #	dat.locus.index[i] <- paste("locus.", locusCounter, sep = "")
  #	if(i / reportInterval - i %/% reportInterval == 0) {
  #	   message(paste('...', i, 'of', length(dat.locus.index), 
  #	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (length(dat.locus.index) - i), attr(Sys.time() - start.time, 'units')
  #	   ))
  #	   }
  #	}

  ## after some vectorization:
  for (i in 1:length(dat.firstLocusLines)) {
    dat.locus.index[dat.firstLocusLines[i]:dat.lastLocusLines[i]] <- names(dat.breakLines.vector)[i] <- paste("locus.", i, sep = "")
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', length(dat.firstLocusLines), 
 	   '-- Estimated time remaining =', round(((Sys.time() - start.time) / i) * (length(dat.firstLocusLines) - i), 1), attr(Sys.time() - start.time, 'units')
  	   ))
	}
  }
  # dat.locus.index <- as.factor(dat.locus.index) # done only for memory considerations... slows things down in analysis unless you transform back to character
  out = list(tips = dat.names, seqs = dat.seqs, breaks = dat.breakLines, break.vectors = dat.breakLines.vector, cons = dat.consensusLines, locus.index = dat.locus.index, file.read = filename, timestamp = date())
  class(out) <- 'pyRAD.loci'
  if(doSummary) out$radSummary <- summary(out)
  return(out)
  }

summary.pyRAD.loci <- function(object, var.only = FALSE, ...) {
# Arguments:
#  object = a pyRAD.loci object
#  var.only = if T, only includes variable loci; as written, the function assumes a "*" if there is 
  message("\nDoing a pyRAD summary")
  reportInterval <- 2000 # this is just for screen reporting... only matters with really long files
  ## currently, locus.names includes a null (""), b/c the break lines have no locus name
  locus.names <- as.character(unique(object$locus.index)) # this slows things down by a factor of 2 or 3, but it seems to prevent a subscript-out-of-bounds error
  locus.names <- locus.names[locus.names != ""]
  ## REWRITE TO LOOK FOR BREAKLINES THAT HAVE * OR - IN THEM
  variable.loci <- locus.names[!is.na(object$seqs[object$breaks])] # note that this works with the pyRAD output we are currently using... should be checked
  if(var.only) locus.names <- variable.loci
  num.loci <- length(locus.names)
  tip.names <- as.character(unique(object$tips[-c(object$breaks, object$cons)]))
  message("Splitting tips by locus name...")
  tips.per.locus <- split(object$tips, object$locus.index)[locus.names]
  seqs.per.locus <- split(object$seqs, object$locus.index)[locus.names]
  num.inds.per.locus <- sapply(tips.per.locus, length)
  inds.mat <- matrix(NA, nrow = length(tip.names), ncol = num.loci, dimnames = list(tip.names, locus.names))
  message("Making tips matrix...")
  start.time <- Sys.time()
  ## is there some way to vectorize the following:
  for(i in seq(num.loci)) {
	temp <- try(tip.names %in% tips.per.locus[[locus.names[i]]])
	if(class(temp) != "try-error") {
	  inds.mat[ , locus.names[i]] <- temp
	  names(seqs.per.locus[[i]]) <- tips.per.locus[[i]]
	  }
	else(message(paste("Error occurred with locus", locus.names[i])))
    if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', num.loci, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (num.loci - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }
	 }

  out <- list(num.loci = num.loci, tips.per.locus = tips.per.locus, break.vectors = object$break.vectors, seqs.per.locus = seqs.per.locus, num.inds.per.locus = num.inds.per.locus, variable.loci = variable.loci, inds.mat = inds.mat, locus.lengths = lengths.report(object, 0))
  class(out) <- 'summary.pyRAD.loci'
  out
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

subset.pyRAD.loci <- function(x, loci, taxa, format = 'DNAStringSet', ...) {
  ## only DNAStringSet export supported now
  out <- list(
    DNA = structure(vector('list', length(loci)), names = loci),
    variable = structure(logical(length(loci)), names = loci)
	)
  inds.vector <- x$tips %in% taxa
  for(i in loci) {
    seq.index <- x$locus.index == i & inds.vector
	out$DNA[[i]] <- DNAStringSet(x$seqs[seq.index])
	names(out$DNA[[i]]) <- x$tips[seq.index == i]
	out$variable[i] <- any(apply(consensusMatrix(out$DNA[[i]])[-c(15:17), ], 2, function(x) sum(x > 0) > 1))
	}
  return(out)
  }
  
write.DNAStringSet <- function(x, format= 'phylip', padding = 30, filename = 'DNAStringSetOut.phy') {
  # writes a sequence matrix to phylip format
  x.width <- width(x)[1]
  x <- as.character(x)
  for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
  writeLines(c(paste(length(x), x.width), x), filename)
  return(0)
  }
  
filter.by <- function(dat, taxa, threshold = 'all') {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat.mat <- dat$radSummary$inds.mat[taxa, ]
    else dat.mat <- dat$inds.mat[taxa, ] # this is the case if you pass in a summary.pyRAD.loci object
  if(threshold == 'all') threshold <- length(taxa)
  return(names(which(apply(dat.mat, 2, sum) == threshold)))
  }

