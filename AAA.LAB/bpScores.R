## bp table for a set of sequences
## tailored to Quercus hybseq, 2021-03-01, ahipp@mortonarb.org

require(magrittr)
require(seqinr)
require(parallel)

bpScores <- function(
                      files, seqNames = NA, format = 'fasta',
                      sortBy = c(amb(), '-', '?'),
                      ncores = 1, simplify = TRUE,
                      seqsAsRows = TRUE,
                      addMatrixStats = TRUE
                    )
{
  if(format != 'fasta') stop('only fasta supported now')
  message(paste('reading', length(files), 'sequence files...'))
  seqs <- lapply(files, read.fasta)
  if(!identical(seqNames, NA)) names(seqs) <- seqNames
  message('doing seqSums...')
  seqSums <- mclapply(seqs, function(x, stats = addMatrixStats) {
    temp <- x %>% unlist %>% table %>% '['(sortBy)
    names(temp) <- sortBy
    temp[is.na(temp)] <- 0
    if(stats) temp <- c(temp, width = length(x[[1]]), nSeq = length(x))
    temp
  }, mc.cores = ncores)
  if(simplify) seqSums <- simplify2array(seqSums)
  if(seqsAsRows) seqSum <- t(seqSums)
  seqSums
}

bpSums <- function(seqSums, classes = list(
    nucs = c("a", "c", "g", "t"),
    ambs = c("u", "r", "y", "m", "k", "s", "w", "b", "d", "h", "v"),
    indet = c("n", "-", "?")),
  props = TRUE, roundTo = 4, addMatrixStats = TRUE
) {
  testL <- ifelse(addMatrixStats,
                  2 + classes %>% unlist %>% length,
                  classes %>% unlist %>% length)
  if(dim(seqSums)[2] != testL) seqSums <- t(seqSums)
  if(dim(seqSums)[2] != testL) stop('array dims wrong')
  if(dim(seqSums)[1] == dim(seqSums)[2])
    warning('square array: make sure it is oriented correctly')
  out <- matrix(NA, dim(seqSums)[1], length(classes),
                dimnames = list(row.names(seqSums), names(classes)))
  for(i in names(classes)) out[, i] <- apply(seqSums[, classes[[i]]], 1, sum)
  if(props) out[, names(classes)] <-
    apply(out[, names(classes)], 1, function(x) round(x / sum(x), roundTo)) %>%
      t
  if(addMatrixStats) out <- cbind(out, seqSums[, c('width', 'nSeq')])
  out <- data.frame(out)
  out
}
