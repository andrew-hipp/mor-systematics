## bp table for a set of sequences
require(magrittr)
require(seqinr)

bpScores <- function(
                      files, seqNames = NA,
                      format = 'fasta', ambProp = TRUE,
                      sortBy = c(amb(), '-', '?')
                    )
{
  if(format != 'fasta') stop('only fasta supported now')
  seqs <- lapply(files, read.fasta)
  if(!identical(seqNames, NA)) names(seqs) <- seqNames
  seqSum <- sapply(seqs, function(x) x %>% unlist %>% table %>%'['(sortBy)
  names(seqSum) <- sortBy
  seqSum
}

## example
