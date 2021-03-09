## read a structure file

read.structure <- function(x) {
  x.dat <- readLines(x)
  startLine <- grep('Inferred clusters', x.dat) + 1
  endLine <- which(x.dat == '')[which(x.dat == '') > startLine][1]
  c(startLine, endLine)
}
