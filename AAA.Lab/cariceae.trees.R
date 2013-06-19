## functions for parsing and interpreting Cariceae trees
## Andrew Hipp, 2013-06-19

tip.spp <- function(tree, delim = '_') {
## finds spp at tips, assuming first two elements are the genus and sp epithet
  sapply(tree$tip.label, function(x) paste(strsplit(x, '_', fixed = T)[[1]][1:2], collapse = ' '))
  }

spp.ci <- function(tree, ...) {
## finds CI of all species
  require(phangorn)