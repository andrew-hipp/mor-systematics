## tools for extracting datasets
## A Hipp, 2013-11-05

require(ape)

make.all.cariceae.dna <- function(base.dir = getwd(), rax.string = 'raxmlHPC-PTHREADS -f a -x 12345 -m GTRCAT -# 200 -T 4') {
  dat <- read.cariceae.data(NULL, NULL, NULL) # right now, just gets all data iteratively
  phy.files.string <- directories.string <- outfiles.string <- character(0)
  for(i in dat$sequence.owners) {
    seq.count <- numeric(0)
	working.dir <- paste(base.dir, '/', i, format(Sys.time(), ".%Y-%m-%d_%H.%M"), sep = '')
	dir.create(working.dir)
	seq.dat <- read.cariceae.data(dat, source.labs = i)
	for(j in seq(length(seq.dat$seqs))) {
	  fasta.file.out <- paste(working.dir, '/', names(seq.dat$seqs)[j], '.fas', sep = '')
	  phy.file.out <- paste(working.dir, '/', names(seq.dat$seqs)[j], '.phy', sep = '')
	  a = try(write.dna(seq.dat$seqs[[j]], fasta.file.out, format = 'fasta'))
      if(class(a) == 'try-error') message(paste('error in writing dataset', j, 'on data provider', i)) 
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
		write.dna(do.call(cbind.DNAbin, seq.dat$seqs[seq.count]), phy.file.out, format = 'sequential', nbcol = -1, colsep = '')
		phy.files.string <- c(phy.files.string, phy.file.out)
		outfiles.string <- c(outfiles.string, 'combinedAnalysis')
		directories.string <- c(directories.string, working.dir)
	    } # close if
	  write.csv(seq.dat$dat, paste(working.dir, '/', i, '.metadata.csv', sep = ''))
	} # close i
  ## now write the raxml batch file
  writeLines(paste(rax.string, ' -s "', phy.files.string, '" -w "', directories.string, '/" -n analysis.', outfiles.string, seq(length(rax.string)), '.out', sep = ''), format(Sys.time(), "analysis.%Y-%m-%d_%H.%M.bat"))
  return(0)
  }

read.cariceae.data <- function(read.dat.obj = NULL, fasta = NULL, metadata = NULL, source.labs = 'ALL_SEQUENCES') {
  if(!is.null(read.dat.obj)) {
    fasta <- read.dat.obj$seqs
	metadata <- read.dat.obj$dat
	}
  if(is.null(fasta)) {
    fasta.file.names <- choose.files(caption = 'Select fasta files', multi= T)
	fasta = lapply(fasta.file.names, read.dna, format= 'fasta')
	names(fasta) <- sapply(strsplit(fasta.file.names, '\\', fixed = T), function(x) return(tail(x,1)))}
  if(is.null(metadata)) metadata <- read.delim(choose.files(caption = 'Select extractions metadata', multi = F), as.is = TRUE)
  metadata$seqName <- paste(metadata$TAXON, metadata$DNA_TUBE_LABELS, sep = '_')
  sequence.owners <- c('ALL_SEQUENCES', as.character(sort(unique(metadata$SOURCE.LAB..owner.of.DNA.))))
  sequence.owners <- sequence.owners[sequence.owners != ''] # get rid of blanks
  if(!source.labs %in% sequence.owners) source.labs <- select.list(sequence.owners, multi = TRUE, title = 'Select source lab(s)')
  if(source.labs != 'ALL_SEQUENCES') {
    # seqs.to.use <- unique(unlist(sapply(source.labs, grep, x = metadata$SOURCE.LAB..owner.of.DNA., value = TRUE)))
	seqs.to.use <- metadata$seqName[metadata$SOURCE.LAB..owner.of.DNA. %in% source.labs]
    for(i in names(fasta)) {
	  ## fasta[[i]] <- fasta[[i]][tidyName(row.names(fasta[[i]])) %in% tidyName(seqs.to.use), ]
	  fasta[[i]] <- fasta[[i]][unique(unlist(sapply(tidyName(seqs.to.use), grep, x = tidyName(row.names(fasta[[i]]))))), ] 
      } # close i	
	metadata <- metadata[metadata$SOURCE.LAB..owner.of.DNA. %in% source.labs, ]
	}
  out <- list(seqs = fasta, dat = metadata, sequence.owners = sequence.owners[-1])
  out
  }