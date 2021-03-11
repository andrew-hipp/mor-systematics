## name clusters in structure output
## args:
##  x = output from read.structure
##  nameList = list of indiduals by sp

str.rename <- function(x, nameList, doSort = TRUE) {
  for(clustCol in 2:length(x)) {
    temp = sapply(spList, function(y) mean(x[y, clustCol]))
    names(x)[clustCol] <- names(nameList)[which(temp == max(temp))[1]]
  }
  names(x) <- make.unique(names(x))
  if(doSort) x <- x[order(names(x))]
  return(x)
}
