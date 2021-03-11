## read a structure file

read.structure <- function(x, fileType = c('structure', 'faststructure'),
                            sampleNames = NULL) {
  if(fileType == 'structure') {
    a <- readLines(x)
    startLine <- grep('Inferred clusters', a) + 1
    endLine <- which(a == '')[which(a == '') > startLine][1] - 1
    writeLines(a[startLine:endLine], 't.temp.read.structure.swapFile.txt')
    a <- read.table('t.temp.read.structure.swapFile.txt', as.is = T)
    names(a)[[2]] <- 'sample'
    a[c(1,3,4)] <- NULL
  }
  if(fileType == 'faststructure') {
    a <- read.table(x)
    if(!identical(sampleNames, NULL)) a <- cbind(sampleNames, a)
  }
  names(a)[2:dim(a)[2]] <- paste('cluster',1:(dim(a)[2] - 1), sep = '')
  class(a) <- c('strObj', 'data.frame')
  return(a)
}
