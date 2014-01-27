locus.picker <-
function(pyDat, minThreshold = 3, inds = row.names(pyDat$radSummary$inds.mat)) {
  out = names(which(colSums(pyDat$radSummary$inds.mat[inds, ]) >= minThreshold))
  }
