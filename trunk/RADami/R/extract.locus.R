extract.locus <- function(pyDat, locName, dat.format = c('text', 'fasta', 'DNAStringSet')) {
  seqs <- pyDat$seqs[pyDat$locus.index == locName]
  if(dat.format[1] == 'text') names(seqs) <- pyDat$tips[pyDat$locus.index == locName]
  as.matrix(seqs)
}