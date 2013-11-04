## tools for extracting datasets
## A Hipp, 2013-11-05

require(ape)

read.cariceae.data <- function(read.dat.obj = NULL, fasta = NULL, metadata = NULL) {
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
  source.labs <- select.list(c('ALL_SEQUENCES', as.character(sort(unique(metadata$SOURCE.LAB..owner.of.DNA.)))), multi = TRUE, title = 'Select source lab(s)')
  if(source.labs != 'ALL_SEQUENCES') {
    seqs.to.use <- metadata$seqName[metadata$SOURCE.LAB..owner.of.DNA. %in% source.labs]
    for(i in names(fasta)) fasta[[i]] <- fasta[[i]][tidyName(row.names(fasta[[i]])) %in% tidyName(seqs.to.use), ]
	metadata <- metadata[metadata$SOURCE.LAB..owner.of.DNA. %in% source.labs, ]
	}
  out <- list(seqs = fasta, dat = metadata)
  out
  }