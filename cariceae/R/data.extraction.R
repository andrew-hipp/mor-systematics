## tools for extracting datasets
## A Hipp, 2013-11-05
## 2014-04-09: added full dataset, added header for analysis file and path prefix for absolute paths, got rid of quotes around files
## 2015-04-27: rewriting to use specimen and extraction databases

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

dna.to.spm <- function(x, dnaDat, 
                         col.spm = 'SPMCODE',
                         col.extraction = 'MOR_DNA_TUBE_NO_CODE',
						 col.spmDNAindex = 'Specimen_ID_CODE',
						 use.tidyName = TRUE, 
						 ...) {
## x is a vector of extraction codes
  dnaRows <- match(tidyName(x, ...), tidyName(dnaDat[[col.extraction]], ...))
  spmID <- dnaDat[dnaRows, col.spmDNAindex]
  return(spmID)
  }
  
  
read.carex.data <- function(read.dat.obj = NULL,
                               additional = NULL,
							   select.by = c('pattern', 'grep'),
							   exclude.permissionNotGranted = TRUE,
							   dna.out = c('phylip', 'fasta'),
							   source.labs = 'ALL_SEQUENCES', 
							   col.owner = 'Ownership_of_Sequence', 
							   col.taxon = 'TAXA-Current_determination',
							   col.tubeNo = 'Original_Tube_No',
							   col.spm = 'SPMCODE',
							   col.permission = 'PERMISSION_TO_USE',
							   col.extraction = 'MOR_DNA_TUBE_NO_CODE',
							   col.spmDNAindex = 'Specimen_ID_CODE',
							   change.tip.labels = TRUE, tail.to = 3, patt = 1:3,
							   tip.label = c('TAXA.Current_determination', 'Ownership_of_Sequence', 'Original_Tube_No', 'Country', 'SPMCODE'),
							   tip.spaceSub = "_",
							   tip.delim = "|"
							   ) {
## 2015-04-27: Several fundamental changes:
##             * this function now goes to a specimen table and a series of DNA tables to get data
##             * only the DNA code is used to pull sequences from the fasta files
##             * the label is written on the fly based on a formula
## 2015-07-06: updated to take CSV or TSV files
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
    spec.file.name <- choose.files(caption = 'Select a single specimens metadata table', multi = FALSE)
	if(tail(strsplit(spec.file.name, '.', fixed = TRUE)[[1]], 1) %in% c('csv', 'CSV')) dat.specimens <- read.csv(spec.file.name, as.is = TRUE)
	else dat.specimens <- read.delim(spec.file.name, as.is = TRUE)
	
	## get extractions data
	dat.extractions <- lapply(choose.files(caption = "Select one or more extractions metadata tables", multi = TRUE), read.delim, as.is = TRUE)
	# sapply(dat.extractions, function(x) print(names(x)))
	dat.extractions <- do.call(rbind, dat.extractions) ## may need to check columns for naming
	} # end else, the reading function if a read.dat.obj was not brought in
		
  dat.specimens$seqName <- apply(dat.specimens[tip.label], 1, function(x) paste(gsub(" ", tip.spaceSub, x), collapse = tip.delim))
  
  ## these lines just subset by ownership -- NOT UPDATED 2015.04.29
  sequence.owners <- c('ALL_SEQUENCES', sort(unique(toupper(dat.specimens[[col.owner]]))))
  sequence.owners <- sequence.owners[sequence.owners != ''] # get rid of blanks
  if(!source.labs %in% sequence.owners) source.labs <- select.list(sequence.owners, multi = TRUE, title = 'Select source lab(s)')
  spms.to.use <- dat.specimens[which(dat.specimens[[col.owner]] %in% source.labs), col.spm]
  if(source.labs != 'ALL_SEQUENCES') {
    seqs.to.use <- dat.extractions[which(dat.extractions[[col.spmDNAindex]] %in% spms.to.use), col.extraction]
    if(!is.null(additional)) seqs.to.use <- unique(c(seqs.to.use, additional)) # does not currently check whether additional is in fasta files
	for(i in names(dat.fasta)) {
	  fasta.tube.codes.extracted <- sapply(row.names(dat.fasta[[i]]), function(x) paste(tail(strsplit(x, '_')[[1]], tail.to)[patt], collapse = "_"))
	  if(select.by[1] == 'grep') dat.fasta[[i]] <- dat.fasta[[i]][unique(unlist(sapply(tidyName(seqs.to.use), grep, x = tidyName(row.names(dat.fasta[[i]]))))), ]
      if(select.by[1] == 'pattern') dat.fasta[[i]] <- dat.fasta[[i]][tidyName(fasta.tube.codes.extracted) %in% tidyName(seqs.to.use) , ]
	  } # close i	
	metadata <- metadata[metadata[[col.owner]] %in% source.labs, ]
	}
	
  for(i in names(dat.fasta)) {
    extracted.spm.codes <- dna.to.spm(sapply(row.names(dat.fasta[[i]]), function(x) paste(tail(strsplit(x, '_')[[1]], tail.to)[patt], collapse = "_")), dat.extractions)
	errorLog <- character(0)
	if(exclude.permissionNotGranted) {
	  reject.rows <- which(!as.logical(dat.specimens[match(extracted.spm.codes, dat.specimens[[col.spm]]), col.permission]))
	  errorLog <- c(errorLog, 'SPECIMENS FLAGGED AS NOT TO BE SHARED', paste(row.names(dat.fasta[[i]])[reject.rows], '[from fasta file] --', dat.specimens[match(extracted.spm.codes, dat.specimens[[col.spm]]), 'seqName'][reject.rows], '[as relabelled]'), '', '')
	  if(length(reject.rows) > 0) {
	    dat.fasta[[i]] <- dat.fasta[[i]][-reject.rows, ]
		extracted.spm.codes <- extracted.spm.codes[-reject.rows]
		} # end if length
	  } # end if exclude
	if(change.tip.labels) {
      new.row.names <- dat.specimens[match(extracted.spm.codes, dat.specimens[[col.spm]]), 'seqName']
	  errorLog <- c(errorLog, 'FASTA.LABEL.NO.SPECIMEN.MATCH', row.names(dat.fasta[[i]])[is.na(new.row.names)], '', '')
	  new.row.names[is.na(new.row.names)] <- paste(row.names(dat.fasta[[i]])[is.na(new.row.names)], "FASTA.LABEL.NO.SPECIMEN.MATCH", sep = tip.delim) # for any failed names, adds in original name
	  errorLog <- c(errorLog, 'FASTA.LABEL.NO.RECOGNIZED.TAXON', row.names(dat.fasta[[i]])[substr(new.row.names,1,1) == tip.delim])
	  new.row.names[substr(new.row.names,1,1) == tip.delim] <- paste(row.names(dat.fasta[[i]])[substr(new.row.names,1,1) == tip.delim], ".FASTA.LABEL.NO.RECOGNIZED.TAXON", new.row.names[substr(new.row.names,1,1) == tip.delim], sep = "") # when there is no species, swap in fasta label
	  write.csv(cbind(oldName = row.names(dat.fasta[[i]]), newName = new.row.names), paste('nameChanges.', i, '.csv', sep = ''))
	  row.names(dat.fasta[[i]]) <- new.row.names
	  writeLines(errorLog, con = paste('errorLog.', i,'.log', sep = ''))
	  }
	}
  dat.fasta$concat.gappy <- do.call(cbind.DNAbin, args = c(dat.fasta, fill.with.gaps = TRUE))
  dat.fasta$concat.noGaps <- do.call(cbind.DNAbin, args = c(dat.fasta, fill.with.gaps = FALSE))
  if(!is.na(dna.out[1])) {
    for(i in names(dat.fasta)) {
	  if(dna.out[1] == 'fasta') write.dna(dat.fasta[[i]], paste(i, '.GCG.export.fas', sep = ''), format = 'fasta')
	  if(dna.out[1] == 'phylip') write.dna(dat.fasta[[i]], paste(i, '.GCG.export.phy', sep = ''), format = 'sequential', nbcol = -1, colsep = '')
	  } # close i
	} # close if
  out <- list(seqs = dat.fasta, dat.specimens = dat.specimens, dat.extractions = dat.extractions, sequence.owners = sequence.owners[-1])
  class(out) <- "carex.data"
  out
  }