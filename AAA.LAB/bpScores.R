## bp table for a set of sequences
require(magrittr)
require(seqinr)

bpScores <- function(
                      files, seqNames = NA, format = 'fasta',
                      sortBy = c(amb(), '-', '?')
                    )
{
  if(format != 'fasta') stop('only fasta supported now')
  seqs <- lapply(files, read.fasta)
  if(!identical(seqNames, NA)) names(seqs) <- seqNames
  seqSum <- sapply(seqs, function(x) {
    temp <- x %>% unlist %>% table %>% '['(sortBy)
    names(temp) <- sortBy
    })
  seqSum
}

## example
bpSums <- function(scoreMat){}
