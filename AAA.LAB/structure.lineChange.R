## change from one-line to two-line STRUCTURE format
## currently v. simple, no pops or extra columns

str1to2 <- function(x, write.str = NULL, ...) {
  if(!class(x) %in% c('matrix', 'data.frame')) {
    infile <- as.matrix(read.table(x, row.names = 1, as.is = TRUE))
  } else infile <- as.matrix(x)
  line1 <- seq(from = 1, to = (dim(infile)[2] - 1), by = 2)
  line2 <- seq(from = 2, to = dim(infile)[2], by = 2)
  outfile <- matrix(NA, dim(infile)[1] * 2, dim(infile)[2] / 2)
  for(i in 1:dim(infile)[1]) {
    outfile[(i * 2 - 1), ] <- infile[i, line1]
    outfile[(i * 2), ] <- infile[i, line2]
  }
  newLabels <- unlist(strsplit(paste(row.names(infile), row.names(infile)), ' '))
  outfile <- cbind(newLabels, outfile)
  if(!identical(write.str, NULL)) {
    write.table(outfile, write.str, quote = F, row.names = F, col.names = F)
  }
  return(outfile)
}
