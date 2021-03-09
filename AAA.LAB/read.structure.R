## read a structure file

read.structure <- function(x) {
  a <- readLines(x)
  startLine <- grep('Inferred clusters', a) + 1
  endLine <- which(a == '')[which(a == '') > startLine][1] - 1
  writeLines(a[startLine:endLine], 't.tempqperopweru.txt')
  a <- read.table('t.tempqperopweru.txt', row.names = 2, as.is = T)
  a[1:3] <- NULL
  return(a)
}
