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
