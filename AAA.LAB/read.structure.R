## read a structure file

read.structure <- function(x) {
  a <- readLines(x)
  startLine <- grep('Inferred clusters', a) + 1
  endLine <- which(a == '')[which(a == '') > startLine][1] - 1
  writeLines(a[startLine:endLine], 't.tempqperopweru.txt')
  a <- read.table('t.tempqperopweru.txt', as.is = T)
  names(a)[[2]] <- 'sample'
  a[c(1,3,4)] <- NULL
  names(a)[2:dim(a)[2]] <- paste('cluster',1:(dim(a)[2] - 1, sep = '')
  class(a) <- 'strObj'
  return(a)
}
