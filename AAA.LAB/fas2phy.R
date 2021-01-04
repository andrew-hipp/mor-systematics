## batch process a bunch of fasta files into phylip, make analysis batch file
## A Hipp, 2013-06-18

get.files <- function(working.dir = choose.dir(), get.patt = '.fas', ignore.patt = NA, ...) {
  out <- list(working.dir = working.dir, files = dir(working.dir,patt = get.patt, full = TRUE))
  if(!is.na(ignore.patt[1])) out$files <- out$files[-grep(paste(ignore.patt, collapse = '|'), out$files, ...)]
  return(out)
  }

fas2phy <- function(fileList = get.files(), suffix = '.phy', analysis.subdir = paste('raxml.', format(Sys.time(), "%Y-%m-%d"), sep = ''), rax.batch.filename = paste('raxml.', format(Sys.time(), "%Y-%m-%d.bat"), sep = '')) {
  require(ape)
  original.dir <- getwd()
  files <- fileList$files
  setwd(fileList$working.dir)
  if(!analysis.subdir %in% dir()) dir.create(analysis.subdir)
  lapply(files, function(x) write.dna(read.dna(x, format = 'fasta'), paste(analysis.subdir, '/', x, suffix, sep = ''), 
                                      format = 'sequential', nbcol = -1, colsep = ''))
  setwd(analysis.subdir)
  subsubdirs <- sapply(files, function(x) strsplit(x, pattern, fixed = T)[[1]][1])
  lapply(subsubdirs, dir.create)
  rax.out <- paste('raxmlHPC-PTHREADS -f a -x 12345 -m GTRCAT -# 200 -T 2 -s', paste(files, suffix, sep = ''), '-w', subsubdirs, '-n', paste(files, '.raxd.tre', sep = ''))
  writeLines(rax.out, rax.batch.filename)
  setwd(original.dir)
  return('done!')
  }

concat.dna <- function(file.list, report = TRUE, delim = '.', ...) {
## dat is a list of DNA objects
  require(ape)
  dat <- lapply(file.list, read.dna, format = 'fasta')
  names(dat) = sapply(file.list, function(x) strsplit(x, delim, fixed = TRUE)[[1]][1])
  if(report) dna.summary <- dna.summary.report(structure(lapply(dat, row.names), names = names(dat)))
  out <- do.call('cbind.DNAbin', c(dat, fill.with.gaps = TRUE))
  if(report) out <- list(dna = out, dna.summary = dna.summary)
  return(out)
  }

dna.summary.report <- function(names.list) {
## takes a list of names and makes a matrix for who has what; sums by columns and rows
  all.names <- sort(unique(unlist(names.list)))
  out <- do.call('cbind', lapply(names.list, function(x) all.names %in% x))
  out <- cbind(out, totals = apply(out, 1, sum))
  out <- rbind(out, totals = apply(out, 2, sum))
  dimnames(out) <- list(c(all.names, 'total'), c(names(names.list), 'total'))
  return(out)
  }