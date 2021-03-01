## bp table for a set of sequences
require(magrittr)
require(seqinr)
require(parallel)
bpScores <- function(
                      files, seqNames = NA, format = 'fasta',
                      sortBy = c(amb(), '-', '?'),
                      ncores = 1, simplify = TRUE
                    )
{
  if(format != 'fasta') stop('only fasta supported now')
  seqs <- lapply(files, read.fasta)
  if(!identical(seqNames, NA)) names(seqs) <- seqNames
  seqSum <- mclapply(seqs, function(x) {
    temp <- x %>% unlist %>% table %>% '['(sortBy)
    names(temp) <- sortBy
    temp[is.na(temp)] <- 0
    temp
  }, mc.cores = ncores)
  if(simplify) seqSum <- simplify2array(seqSum)
  seqSum
}

## example
bpSums <- function(scoreMat){}
