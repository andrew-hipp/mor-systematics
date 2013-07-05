## functions for parsing and interpreting Cariceae trees
## Andrew Hipp, 2013-06-19

## variables
top12 <- c('ITS', 'ETS','trnLF','rbcL','matK','18S','trnK','rpoC1','rps16','psbA','atpF','rpoB')
top5 <- c('ITS', 'ETS','trnLF','rbcL','matK')
dregs <- c('trnK','rpoC1','rps16','psbA','atpF','rpoB')

summarize.gene.coverage <- function(datMat = cariceae.concat.4genes.summary, extra.left = 2, ...) {
  genes.no = dim(datMat)[2]
  inds.no = dim(datMat)[1]
  par(mar = c(5, 4 + extra.left, 4, 2) + 0.1)
  plot(0, type = "n", xlim = c(1,inds.no), ylim = c(0, genes.no+1), xlab = 'Individual sequenced [number]', ylab = '', axes = FALSE)
  axis(2, at = seq(1:genes.no), labels = names(datMat), tick = FALSE, las = 1)
  axis(1)
  for(i in 1:genes.no) points(which(datMat[1:inds.no, i] == 1), rep(i, sum(datMat[1:inds.no, i])), col = i, ...)
  return('done')
  }
  
## tree traversal and counting
tip.spp <- function(tree, delim = '[_.]', ...) {
## finds spp at tips, assuming first two elements are the genus and sp epithet
  sapply(tree$tip.label, function(x) paste(strsplit(x, delim, ...)[[1]][1:2], collapse = ' '))
  }

spp.ci <- function(tree, ...) {
## finds CI of all species
  require(phangorn)
  ## NOT DONE
  }
  
monophyletic.spp <- function(tree, ...) {
  require(ape)
  all.spp <- tip.spp(tree, ...)
  unique.spp <- sort(unique(all.spp))
  out <- list(numSpp = length(unique.spp), spp.summary = cbind(count = sapply(unique.spp, function(x) sum(all.spp == x)), monophyletic = sapply(unique.spp, function(x) is.monophyletic(tree, names(all.spp)[all.spp == x]))))
  return(out)
  }