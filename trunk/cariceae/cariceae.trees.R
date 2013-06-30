## functions for parsing and interpreting Cariceae trees
## Andrew Hipp, 2013-06-19

## variables
top12 <- c('ITS', 'ETS','trnLF','rbcL','matK','18S','trnK','rpoC1','rps16','psbA','atpF','rpoB')
top5 <- c('ITS', 'ETS','trnLF','rbcL','matK')
dregs <- c('trnK','rpoC1','rps16','psbA','atpF','rpoB')

## tree traversal and counting
tip.spp <- function(tree, delim = '_') {
## finds spp at tips, assuming first two elements are the genus and sp epithet
  sapply(tree$tip.label, function(x) paste(strsplit(x, '_', fixed = T)[[1]][1:2], collapse = ' '))
  }

spp.ci <- function(tree, ...) {
## finds CI of all species
  require(phangorn)