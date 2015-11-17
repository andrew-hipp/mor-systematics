## functions for parsing and interpreting Cariceae trees
## Andrew Hipp, 2013-06-19
## additions 2015-03-03 to deal with sections

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
label.elements <- function(tree, delim = '|', returnNum = 1, returnDelim = ' ', ...) {
## finds any label at the tips; default assumes pipe delimitation, and the element of interest is b/f the first pipe
  out <- sapply(tree$tip.label, function(x) paste(strsplit(x, delim, ...)[[1]][returnNum], collapse = returnDelim))
  out
  }

tip.select <- function(tr, ...) {
## this function returns the number of nodes that would need to be collapsed or traversed or something to make a set of tips monophyletic; these can be defined by species, sections, or any other attribute; could also return CI
}

tips.ci <- function(tr, tipsBoolean, ...) {
## finds CI of any boolean tip vector
  require(phangorn)
  tipsBoolean <- as.integer(tipsBoolean)
  tips <- phyDat(matrix(tipsBoolean, length(tipsBoolean), 1, dimnames = list(tr$tip.label, NULL)), type = 'USER', levels = 0:1)
  out <- CI(tr, tips)
  out
  }

tips.expected <- function(tr, tips, value = FALSE, ...) {
## returns the number of tips descended from the mrca of a vector of tips
## called "expected" because this is the number of tips expected if the tips provided were a clade
  require(phangorn)
  if(length(tips) == 1) return(1)
  tips.mrca <- getMRCA(tr, tips)
  tips.descendants <- Descendants(tr, tips.mrca, ...)[[1]]
  if (value) out <- tips.descendants
  else out <- length(tips.descendants)
  out 
  }
  
summary.by.elements <- function(tr, ...) {
## this is the taxonomic disparity function
  all.elements <- label.elements(tr, ...)
  unique.elements <- sort(unique(all.elements))
  out <- cbind(count = sapply(unique.elements, function(x) sum(all.elements == x)),
               expected = sapply(unique.elements, function(x) tips.expected(tr, names(all.elements)[all.elements == x]))
			   )
  out <- cbind(out, disparity = out[, 'expected'] - out[, 'count'])
  out
  }
  
monophyletic.spp <- function(tree, ...) {
  require(ape)
  all.spp <- label.elements(tree, ...)
  unique.spp <- sort(unique(all.spp))
  out <- list(numSpp = length(unique.spp), 
              spp.summary = cbind(count = sapply(unique.spp, function(x) sum(all.spp == x)), 
			                      monophyletic = sapply(unique.spp, function(x) is.monophyletic(tree, names(all.spp)[all.spp == x])),
								  ci = sapply(unique.spp, function(x) tips.ci(tree, all.spp == x))
								  ))
  return(out)
  }

add.data.to.tips <- function(tr, datMat, delim = '[_|]', returnNum = 1:2, returnDelim = "_", addCols = c('GROUP'), addDelim = "|", ...) {
## adds label to the end of the tips using a standard formula
## returns the relabelled tree and a matrix indicating what tip got what section
## arguments:
##   tr = tree
##   datMat = matrix with the data you want added as columns, row.names matching tips of the tree after scrubbing through label.elements
  tips.to.match <- label.elements(tr, delim, returnNum, returnDelim, ...)
  addVect <- datMat[match(tips.to.match, row.names(datMat)), addCols]
  if(!is.null(dim(addVect))) addVect <- apply(addVect, 1, paste, collapse = addDelim)
  oldNames <- tr$tip.label
  tr$tip.label <- paste(oldNames, addVect, sep = addDelim)
  out = list(tr.relabelled = tr, labelMat = cbind(oldLabel = oldNames, tipMatched = tips.to.match, newElement = addVect, newLabel = tr$tip.label))
  return(out)
}

color.tips.by.element <- function(tr, element = 6, delim = "|", fixed = TRUE, whiteOut = 'NA', addLegend = T, colorIt = FALSE, byLabels = TRUE, tip.cex = 0.1, dot.pch = 1, blank.tips = T, ...) {
  vectorToColorBy <- label.elements(tr, delim, returnNum = element, fixed = fixed)
  #tr$tip.label <- vectorToColorBy
  colors = colors()[as.factor(vectorToColorBy)]
  colors[vectorToColorBy %in% whiteOut] <- 'white'
  par(mar = c(5,10,4,2))
  if(blank.tips) tr$tip.label <- sapply(tr$tip.label, function(x) "")
  a = plot(tr, ...)
  if(colorIt) tiplabels(col = colors, pch = dot.pch, cex = tip.cex)
  if(byLabels) tiplabels(vectorToColorBy, cex = tip.cex, align.tip.label = T)
  if(addLegend) legend(a$x.lim[1] - abs(diff(a$x.lim) / 4), a$y.lim[2], legend = unique(vectorToColorBy), pch = dot.pch, cex = 1, col = unique(colors), bty = 'n')
  }



