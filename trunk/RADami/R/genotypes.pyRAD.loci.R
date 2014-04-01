genotypes.pyRAD.loci <- function(dat, groups = list(lobatae = inds.lobatae, quercus = inds.quercus),
                                 loci = 'all', taxa = 'all', useSnps = c('first', 'all'), concat = c(FALSE, TRUE), 
								 use.tidyName = TRUE, na.rm = c('columns', 'rows', 'none'), maxAlleles = 2, 
								 tidyVals = c('-', '.','>', '_', ' ', 'oak'), sortByGroups = TRUE, cores = 8) {
##  Makes a dataframe of SNP calls from a pyRAD.loci object for export to hierfstat
##  arguments:
##    dat = currently requires a subset.pyRAD.loci object
##    groups = list of groups of individuals,
##    taxa = taxa to include
##    loci = loci to include
##    useSnps = whether to use first or all SNPs per RAD locus (not yet implemented)
##    concat = whether to concatenate loci or leave separated as in DAT; not currently implemented

  if(!'subset.pyRAD.loci' %in% class(dat)) stop('Currently, this function is written to require DNAStringSet output from subset.pyRAD.loci,\n
                                                 with only SNPs exported')
  if(taxa[1] != 'all') {
    dat$DNA <- lapply(dat$DNA, function(x) x[names(x) %in% taxa])
	dat$DNA <- dat$DNA[sapply(dat$DNA, length) > 0]
	}
  if(loci[1] != 'all') dat$DNA <- dat$DNA[loci]
  out <- structure(vector('list', length(dat$DNA)), names = names(dat$DNA))
  duplicated.members <- unlist(groups)[duplicated(unlist(groups))]
  if(length(duplicated.members) > 0) warning('Some individuals are duplicated between groups; excluding duplicates from export, including first')
  groups.vector <- structure(integer(length(unlist(groups)[!duplicated(unlist(groups))])), names = unlist(groups)[!duplicated(unlist(groups))])
  for(i in 1:length(groups)) groups.vector[groups[[i]]][!names(groups.vector[groups[[i]]]) %in% duplicated.members] <- i
  
## 3. Translate SNPs to genotypes
  do.this <- function(y) {
	assign('counter', counter + 1, envir = .GlobalEnv)
	message(paste('Doing data', counter))
	# if(counter == 469) browser()
	trans.dna <- t(apply(as.matrix(y), 1, function(x) IUPAC_CODE_MAP[x]))
	# if(is.null(trans.dna)) trans.dna <- cbind(trans.dna, dummy.column = rep('A', length(trans.dna)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('A', '1', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('C', '2', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('G', '3', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('T', '4', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) {x[nchar(x) == 1] <- paste(x[nchar(x) == 1], x[nchar(x) == 1], sep = ''); return(x)}))
	trans.dna <- as.matrix(trans.dna)[!apply(as.matrix(trans.dna), 1, function(x) any(nchar(x) > maxAlleles)), ]
	if(na.rm[1] == 'rows') trans.dna <- as.matrix(trans.dna)[!apply(as.matrix(trans.dna), 1, function(x) any(is.na(x))), ]
	if(na.rm[1] == 'columns') trans.dna <- as.matrix(trans.dna)[, !apply(as.matrix(trans.dna), 2, function(x) any(is.na(x)))]
	groupMembership <- groups.vector[match(tidyName(row.names(as.matrix(trans.dna)), tidyVals), tidyName(names(groups.vector), tidyVals))]
	if(is.null(dim(trans.dna))) dna.out <- as.data.frame(cbind(groupMembership = groupMembership, as.integer(trans.dna)))
	if(length(dna.out) == 0) return(0)
	else dna.out <- as.data.frame(cbind(groupMembership = groupMembership, t(apply(trans.dna,1,as.integer))))
	row.names(dna.out) <- row.names(trans.dna) # necessary?
	if(sortByGroups) dna.out <- dna.out[order(dna.out$groupMembership), ]
	return(dna.out)
	}
  
  assign("counter", 0, envir = .GlobalEnv)
  
  out <- lapply(dat$DNA, do.this)
  # out <- mclapply(dat$DNA, do.this, mc.cores = cores)
  out <- out[class(out) %in% c('data.frame', 'matrix')] # gets rid of anything that isn't a matrix or a data.frame
  out <- out[!apply(t(sapply(out, dim)), 1, function(x) sum(x == 0) > 0)] # gets rid of all the matrices in which some dimension == 0
  attr(out, 'groupMembership') <- t(sapply(out, function(w) sapply(1:2, function(x) sum(w$groupMembership == x))))
  out
}