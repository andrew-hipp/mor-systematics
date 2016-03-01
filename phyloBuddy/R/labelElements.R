label.elements <- function(x, delim = '|', returnNum = 1, returnDelim = ' ', ...) {
## finds any label at the tips; default assumes pipe delimitation, and the element of interest is b/f the first pipe
  if('phylo' %in% class(x)) labelVector <- x$tip.label
  else labelVector <- x
  out <- sapply(labelVector, function(x) paste(strsplit(x, delim, ...)[[1]][returnNum], collapse = returnDelim))
  out
  }
