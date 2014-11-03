fasta.name.scrubber <- function(x, delim = '_', elements = 3, ...) {
## just scrubs names down to the extraction code, peels off scientific names
## x is a character vector

  x.list <- strsplit(x, delim, ...)
  x.out <- sapply(x.list, function(y) {paste(y[(length(y) - 2):length(y)], collapse = delim)}, USE.NAMES = FALSE)
  x.out
  }

clean.all.names <- function(directory = './SEQS/', pattern = '.fas', ...) {
  newdir <- paste(directory, 'new.seqs.', format(Sys.time(), "%Y-%m-%d"), sep = '')
  dir.create(newdir)
  all.dat <- lapply(dir(directory, patt = pattern, full = T), read.dna, format = 'fasta')
  names(all.dat) <- gsub(pattern, '', dir(directory, patt = pattern))
  for(i in names(all.dat)) {
    out.table <- cbind(old.name = row.names(all.dat[[i]]), 
	                   new.name = fasta.name.scrubber(row.names(all.dat[[i]]), ...))
	row.names(all.dat[[i]]) <- out.table[, 'new.name']
	write.csv(out.table, paste(newdir, '/', i, '.nameTable.csv', sep = ''))
	write.dna(all.dat[[i]], paste(newdir, '/', i, '.renamed.fas', sep = ''), format = 'fasta')
	}
  return(all.dat)
  }

