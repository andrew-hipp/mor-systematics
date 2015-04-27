## tools for extracting datasets
## A Hipp, 2013-11-05
## 2014-04-09: added full dataset, added header for analysis file and path prefix for absolute paths, got rid of quotes around files

require(ape)

post.process.cariceae <- function(base.dir, out.dir = 'tree.files', tree.patt = 'RAxML_bipartitions.analysis.', split.patt = tree.patt, split.side = 2) {
  if(!out.dir %in% dir(base.dir)) dir.create(paste(base.dir, '/', out.dir, sep = ''))
  all.files = dir(base.dir, patt = tree.patt, recursive = TRUE)
  lapply(all.files, function(x) file.copy(from = paste(base.dir, '/', x, sep = ''), to = paste(base.dir, '/', out.dir, '/', strsplit(x, '/')[[1]][length(strsplit(x, '/')[[1]]) - 1], '_', strsplit(x, split.patt)[[1]][split.side], '.tre', sep = '')))
  }
  
make.all.cariceae.dna <- function(base.dir = getwd(), 
                                  rax.threads = 8,
						          path.prefix = '~/Dropbox/NSF-CAREX-1255901-core-files/PRIMARY.DATA/SEQUENCE.DATA/MORTON/Analyses/',
								  file.header = '#!/bin/sh',
								  analysis.file.suffix = '.sh',
								  og = 'Carex_siderosticta_Carex_MOR_2594',
								  rax.path = '~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-PTHREADS-AVX',
								  rax.call = '-f a -x 12345 -p 12345 -m GTRCAT -# 100',
								  ...) {
  rax.string <- paste(rax.path, rax.call, '-T', rax.threads)
  dat <- read.cariceae.data(NULL, NULL, NULL) # right now, just gets all data iteratively
  phy.files.string <- directories.string <- outfiles.string <- character(0)
  for(i in c(dat$sequence.owners, 'all.data.pooled')) {
    seq.count <- numeric(0)
	working.dir <- paste(base.dir, '/', i, format(Sys.time(), ".%Y-%m-%d_%H.%M"), sep = '')
	dir.create(working.dir)
	if(i == 'all.data.pooled') seq.dat <- dat
	  else seq.dat <- read.cariceae.data(dat, source.labs = i, additional = og, append.source = FALSE, tail.to = 4, patt = 1:3)
	for(j in seq(length(seq.dat$seqs))) {
	  fasta.file.out <- paste(working.dir, '/', names(seq.dat$seqs)[j], '.fas', sep = '')
	  phy.file.out <- paste(working.dir, '/', names(seq.dat$seqs)[j], '.phy', sep = '')
	  out.seq <- seq.dat$seqs[[j]]
	  if(dim(out.seq)[1] < 3) {
	    message(paste('fewer than 3 sequences in dataset', j, 'for provider', i, '-- not writing matrices'))
		next
		}
	  a = try(write.dna(out.seq, fasta.file.out, format = 'fasta'), silent = TRUE)
      if(class(a) == 'try-error') message(paste('error in writing dataset', j, 'on data provider', i, '-- presumably no data available, or labels mismatch')) 
	    else {
		  write.dna(seq.dat$seqs[[j]], phy.file.out, format = 'sequential', nbcol = -1, colsep = '')
		  phy.files.string <- c(phy.files.string, phy.file.out)
		  outfiles.string <- c(outfiles.string, ifelse(length(grep('ETS', phy.file.out)) > 0, 'ETS', 'ITS')) ## ASSUMES JUST ETS OR ITS RIGHT NOW
		  directories.string <- c(directories.string, working.dir)
		  seq.count <- c(seq.count, j)
		  } # close else
      } # close j
	  if(length(seq.count) > 1) {
	    phy.file.out <- paste(working.dir, '/combinedData.phy', sep = '')
		write.dna(do.call("cbind.DNAbin", c(seq.dat$seqs[seq.count], fill.with.gaps = TRUE)), phy.file.out, format = 'sequential', nbcol = -1, colsep = '')
		phy.files.string <- c(phy.files.string, phy.file.out)
		outfiles.string <- c(outfiles.string, 'combinedAnalysis')
		directories.string <- c(directories.string, working.dir)
	    } # close if
	  write.csv(seq.dat$dat, paste(working.dir, '/', i, '.metadata.csv', sep = ''))
	} # close i

	## now write the raxml batch file
  out <- paste(rax.string, ' -s ', path.prefix, phy.files.string, ' -w ', path.prefix, directories.string, ifelse(is.null(og), '', paste('-o', og)), ' -n analysis.', outfiles.string, seq(length(rax.string)), '.out', sep = '')
  if(file.header != '') out <- c(file.header, out)
  writeLines(out, con = paste(format(Sys.time(), "analysis.%Y-%m-%d_%H.%M"), analysis.file.suffix, sep = ''))
  return(0)
  }

read.cariceae.data <- function(read.dat.obj = NULL,
                               source.labs = 'ALL_SEQUENCES', additional = NULL,
							   select.by = c('pattern', 'grep', 'identity'),
							   source.col = 'CONTRIBUTOR', append.source = TRUE, tail.to = 3, patt = 1:3,
							   ) {
## 2015-04-27: Several fundamental changes:
##             * this function now goes to a specimen table and a series of DNA tables to get data
##             * only the DNA code is used to pull sequences from the fasta files
##             * the label is written on the fly based on a formula
  if(!is.null(read.dat.obj)) {
    dat.fasta <- read.dat.obj$dat.dna
	dat.specimens <- read.dat.obj$dat.specimens
	dat.extractions <- read.dat.obj$dat.extractions
	} # you can feed a read.cariceae.data object in to start the process
  else {
    ## get fasta data
	fasta.file.names <- choose.files(caption = 'Select fasta files', multi= T)
	dat.fasta <- lapply(fasta.file.names, read.dna, format= 'fasta')
	names(dat.fasta) <- sapply(strsplit(fasta.file.names, '\\', fixed = T), function(x) return(tail(x,1)))
	
	## get specimens data
    dat.specimens <- read.delim(choose.files(caption = 'Select specimens metadata table', multi = F), as.is = TRUE)
	
	## get extractions data
	dat.extractions <- lapply(choose.files(caption = "Select one or more extractions metadata tables", multi = TRUE), read.delim, as.is = TRUE)
	dat.extractions <- do.call(rbind, dat.extractions) ## may need to check columns for naming
	} # end else, the reading function if a read.dat.obj was not brought in
	
	
	#### STOPPED HERE
	
	
  metadata$seqName <- paste(metadata$TAXON, metadata$DNA_TUBE_LABELS, sep = '_')
  sequence.owners <- c('ALL_SEQUENCES', as.character(sort(unique(metadata[[source.col]]))))
  sequence.owners <- sequence.owners[sequence.owners != ''] # get rid of blanks
  if(!source.labs %in% sequence.owners) source.labs <- select.list(sequence.owners, multi = TRUE, title = 'Select source lab(s)')
  if(source.labs != 'ALL_SEQUENCES') {
    # seqs.to.use <- unique(unlist(sapply(source.labs, grep, x = metadata$SOURCE.LAB..owner.of.DNA., value = TRUE)))
	if(select.by[1] %in% c('grep', 'pattern')) seqs.to.use <- metadata$DNA_TUBE_LABELS[metadata[[source.col]] %in% source.labs]
    if(select.by[1] == 'identity') seqs.to.use <- metadata$seqName[metadata[[source.col]] %in% source.labs]
	if(!is.null(additional)) seqs.to.use <- unique(c(seqs.to.use, additional)) # does not currently check whether additional is in fasta files
	for(i in names(fasta)) {
	  ## fasta[[i]] <- fasta[[i]][tidyName(row.names(fasta[[i]])) %in% tidyName(seqs.to.use), ]
	  fasta.tube.codes.extracted <- sapply(row.names(fasta[[i]]), function(x) paste(tail(strsplit(x, '_')[[1]], tail.to)[patt], collapse = "_"))
	  if(select.by[1] %in% c('grep', 'identity')) fasta[[i]] <- fasta[[i]][unique(unlist(sapply(tidyName(seqs.to.use), grep, x = tidyName(row.names(fasta[[i]]))))), ]
      if(select.by[1] == 'pattern') fasta[[i]] <- fasta[[i]][tidyName(fasta.tube.codes.extracted) %in% tidyName(seqs.to.use) , ]
	  } # close i	
	metadata <- metadata[metadata[[source.col]] %in% source.labs, ]
	}
  if(append.source) {
    for(i in names(fasta)) {
      fasta.tube.codes.extracted <- sapply(row.names(fasta[[i]]), function(x) paste(tail(strsplit(x, '_')[[1]], 3), collapse = "_"))
	  row.names(fasta[[i]]) <- paste(row.names(fasta[[i]]), metadata[[source.col]][match(tidyName(fasta.tube.codes.extracted), tidyName(metadata$DNA_TUBE_LABELS))], sep = '_')
	  }
	}
  out <- list(seqs = fasta, dat = metadata, sequence.owners = sequence.owners[-1])
  class(out) <- "cariceae.data"
  out
  }