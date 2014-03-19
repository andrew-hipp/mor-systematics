genotypes.pyRAD.loci <- function(dat, groups = list(lobatae = inds.lobatae, quercus = inds.quercus),
                                 loci = 'all', taxa = 'all', useSnps = c('first', 'all'), concat = c(FALSE, TRUE), 
								 use.tidyName = TRUE, na.rm = c('columns', 'rows', 'none'), maxAlleles = 2, 
								 tidyVals = c('-', '.','>', '_', ' ', 'oak'), sortByGroups = TRUE, ...) {
##  Makes a dataframe of SNP calls from a pyRAD.loci object for export to hierfstat
##  arguments:
##    dat = currently requires a subset.pyRAD.loci object
##    groups = list of groups of individuals,
##    taxa = taxa to include
##    loci = loci to include
##    useSnps = whether to use first or all SNPs per RAD locus (not yet implemented)
##    concat = whether to concatenate loci or leave separated as in DAT; not currently implemented
##    ... = additional arguments passed along to group.pyRAD.loci

  if(!'subset.pyRAD.loci' %in% class(dat)) stop('Currently, this function is written to require DNAStringSet output from subset.pyRAD.loci,\n
                                                 with only SNPs exported')
  out <- structure(vector('list', length(dat$DNA)), names = names(dat$DNA))
  duplicated.members <- unlist(groups)[duplicated(unlist(groups))]
  if(length(duplicated.members) > 0) warning('Some individuals are duplicated between groups; excluding duplicates from export, including first')
  groups.vector <- structure(integer(length(unlist(groups)[!duplicated(unlist(groups))])), names = unlist(groups)[!duplicated(unlist(groups))])
  for(i in 1:length(groups)) groups.vector[groups[[i]]][!names(groups.vector[groups[[i]]]) %in% duplicated.members] <- i
  
## 3. Translate SNPs to genotypes
  for(i in 1:length(dat$DNA)) {
    trans.dna <- t(apply(as.matrix(dat$DNA[[i]]), 1, function(x) IUPAC_CODE_MAP[x]))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('A', '1', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('C', '2', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('G', '3', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) gsub('T', '4', x)))
	trans.dna <- t(apply(trans.dna, 1, function(x) {x[nchar(x) == 1] <- paste(x[nchar(x) == 1], x[nchar(x) == 1], sep = ''); return(x)}))
	trans.dna <- trans.dna[!apply(trans.dna, 1, function(x) any(nchar(x) > maxAlleles)), ]
	if(na.rm[1] == 'rows') trans.dna <- trans.dna[!apply(trans.dna, 1, function(x) any(is.na(x))), ]
	if(na.rm[1] == 'columns') trans.dna <- trans.dna[, !apply(trans.dna, 2, function(x) any(is.na(x)))]
	groupMembership <- groups.vector[match(tidyName(row.names(trans.dna), tidyVals), tidyName(names(groups.vector), tidyVals))]
	out[[i]] <- as.data.frame(cbind(groupMembership = groupMembership, t(apply(trans.dna,1,as.integer))))
	row.names(out[[i]]) <- row.names(trans.dna)
	if(sortByGroups) out[[i]] <- out[[i]][order(out[[i]]$groupMembership), ]
	} # close i
  out <- out[!apply(t(sapply(out, dim)), 1, function(x) sum(x == 0) > 0)] # gets rid of all the matrices in which someone dimension == 0
  attr(out, 'groupMembership') <- t(sapply(out, function(w) sapply(1:2, function(x) sum(w$groupMembership == x))))
  out
}