cross.2830.2893G1 <- list(parents = c('>2830D', '>2830D2'), 
                          f1s = c('>2830Dx2893G1A', '>2830Dx2893G1B','>2830Dx2893G1C','>2830Dx2893G1D', '>2830Dx2893G1E', '>2893Gx2830D1A', '>2893Gx2830D1C')
						  )

							
cross.2816.2893G1 <- list(parents = c('>2816', '>2830D2'),
                          f1s = c(">2816x2893G1E", ">2816x2893G1B", ">2816x2893G1C", ">2816x2893G1D", ">2816x2893G1A")
						  )
							
cross.allOffspring.2816 <- list(parents = c('>2816', '>2893G1'),
                                f1s = c(">2816x2893G1E", ">2816x2893G1B", ">2816x2893G1C", ">2816x2893G1D", ">2816x2893G1A",
								'>2830Dx2893G1A', '>2830Dx2893G1B','>2830Dx2893G1C','>2830Dx2893G1D', '>2830Dx2893G1E', '>2893Gx2830D1A', '>2893Gx2830D1C')
								)

cross.allOffspring.2830 <- list(parents = c('>2830D', '>2893G1'),
                                f1s = c(">2816x2893G1E", ">2816x2893G1B", ">2816x2893G1C", ">2816x2893G1D", ">2816x2893G1A",
								'>2830Dx2893G1A', '>2830Dx2893G1B','>2830Dx2893G1C','>2830Dx2893G1D', '>2830Dx2893G1E', '>2893Gx2830D1A', '>2893Gx2830D1C')
								)
								
hybrid.test <- function(dat, parents, f1s, unambiguousParents = TRUE, silent = TRUE)								 
### TO DO (3/14/2013, AH and AMELirio):
##    don't analyze Ns and -s
##    Check for variability after screening out everyone except for the parents and offspring
								 {
require(Biostrings)
  ## go through all loci and find for each (1) the parent genotypes, 
  ##(2) the expected F1 genotypes and their ratios, 
  ##(3) the observed F1 genotypes and their ratios,
  ##(4) the percent F1 genotypes that are the expected genotypes
  
  variableSiteCharacters <- c("-", "*") ## change this if Deren rewrites pyRAD to use other characters
  if(class(dat) != 'summary.pyRAD.loci') stop('summary of pyRAD data needed here!')
  parentNumbers <- which(dimnames(dat$inds.mat)[[1]] %in% parents)
  f1Numbers <- which(dimnames(dat$inds.mat)[[1]] %in% f1s)
  loci.to.use <- which(apply(dat$inds.mat, 2, function(x) all(x[parentNumbers] == T) & any(x[f1Numbers] == T))) 
  loci.to.use.names <- dimnames(dat$inds.mat)[[2]][loci.to.use]
  loci.to.use.names <- loci.to.use.names[loci.to.use.names %in% dat$variable.loci]
  rm(loci.to.use)
  matsOut <- structure(vector('list', length(loci.to.use.names)), .Names = loci.to.use.names)
  colNames <- character(0)
  for(locusCounter in loci.to.use.names) {
	if(!silent) message(paste("Doing", locusCounter))
	seqsMat <- t(as.matrix(sapply(as.character(dat$seqs.per.locus[[locusCounter]]), function(x) strsplit(x, split = "")[[1]]), byrow = T))   
	dimnames(seqsMat)[[1]] <- dat$tips.per.locus[[locusCounter]]
	seqLength <- dim(seqsMat)[2]
	seqConsensus <- substr(dat$break.vectors[locusCounter], nchar(dat$break.vectors[locusCounter]) - seqLength + 1, nchar(dat$break.vectors[locusCounter]))
    variable.sites <- which(strsplit(seqConsensus, "")[[1]] %in% variableSiteCharacters)
	if(length(variable.sites) == 0) {
	  if(!silent) message("No variable sites... moving on.")
	  next
	  }
	if(!silent) message(paste("Found", length(variable.sites), "variable sites"))
	seqsMat <- as.matrix(seqsMat[c(parents, dimnames(seqsMat)[[1]][dimnames(seqsMat)[[1]] %in% f1s]), variable.sites]) # only includes parents and children
	#browser()
	for (i in 1:dim(seqsMat)[2]) {
	  seqPossibilities <- character(0)
	  parent.sites <- strsplit(IUPAC_CODE_MAP[seqsMat[parents, i]], "")
	  for(j in 1:length(parent.sites[[1]])) {
	    for(k in 1:length(parent.sites[[2]])) {
		  seqPossibilities <- c(seqPossibilities, mergeIUPACLetters(paste(parent.sites[[1]][j], parent.sites[[2]][k], sep = "")))
		  }}
	  seqsMat <- cbind(seqsMat, c(NA, NA, seqsMat[3:dim(seqsMat)[1], i] %in% seqPossibilities))
	  }
	matsOut[[locusCounter]] <- seqsMat
	colNames <- c(colNames, paste(locusCounter, "_", seq(dim(seqsMat)[2] / 2), sep = ""))
	}
  browser()
  if(!silent) message(paste("Columns in summary matrix:", length(colNames)))
  summaryMat <- matrix(NA, nrow = length(c("differ", f1s)), ncol = length(colNames), dimnames = list(c("Parents differ", f1s), colNames))
  colCounter <- 1
  for(i in which(sapply(matsOut, function(x) !identical(x, NULL)))) {
    if(!silent) message(paste('DOING MATRIX', i, 'OF', length(matsOut)))
	for(j in ((dim(matsOut[[i]])[2] / 2) + 1):dim(matsOut[[i]])[2]) {
	  for(k in 3:dim(matsOut[[i]])[1]) {
	    if(!silent) message(paste("Doing summary matrix column", colCounter))
		if(!silent) message(paste(" -- working on row", dimnames(matsOut[[i]])[[1]][k]))
		summaryMat[dimnames(matsOut[[i]])[[1]][k], colCounter] <- as.logical(matsOut[[i]][k,j])
		if(any(matsOut[[i]][1:2,(j - dim(matsOut[[i]])[2] / 2)] == "N")) summaryMat[1, colCounter] <- FALSE
		else summaryMat[1, colCounter] <- matsOut[[i]][1,(j - dim(matsOut[[i]])[2] / 2)] != matsOut[[i]][2,(j - dim(matsOut[[i]])[2] / 2)]
		if(unambiguousParents) {
		  if(any(!matsOut[[i]][1:2,(j - dim(matsOut[[i]])[2] / 2)] %in% c('A', 'C', 'G', 'T'))) summaryMat[1, colCounter] <- FALSE
		  }
		}# close k
		colCounter <- colCounter + 1	  
	  }# close j
	}# close i
  parents.differ <- summaryMat[1, ]
  percent.compatible.with.cross <- apply(summaryMat[, parents.differ], 1, mean, na.rm = TRUE)
  total.scorable.loci <- apply(summaryMat[, parents.differ], 1, function(x) sum(!is.na(x)))
  percent.compatible.with.cross[1] <- sum(parents.differ)
  out = list(matsOut = matsOut, summaryMat = summaryMat, percent.compatible.with.cross = percent.compatible.with.cross, total.scorable.loci = total.scorable.loci)
  return(out)
  }
  