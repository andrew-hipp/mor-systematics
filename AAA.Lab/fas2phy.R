## batch process a bunch of fasta files into phylip, make analysis batch file
## A Hipp, 2013-06-18

fas2phy <- function(working.dir = choose.dir(), pattern = '.fas', ignore = 'none', suffix = '.phy', analysis.subdir = paste('raxml.', format(Sys.time(), "%Y-%m-%d"), sep = ''), rax.batch.filename = paste('raxml.', format(Sys.time(), "%Y-%m-%d.bat"), sep = '')) {
  require(ape)
  original.dir <- getwd()
  setwd(working.dir)
  files <- dir(patt = pattern)
  if(ignore != 'none') files <- files[!files %in% ignore]
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

concat.DNA <- function(dat) {
## dat is a list of DNA objects
}