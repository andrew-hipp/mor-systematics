# helper functions for dealing with short-read data
# Andrew Hipp, December 2012

str2mat <- function(seqs) t(sapply(seqs, function(x) strsplit(x, "")[[1]])) # turns a vector of sequences into a matrix, one nucleotide per cell

