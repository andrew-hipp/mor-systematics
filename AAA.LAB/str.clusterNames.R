## name clusters in structure output
## args:
##  x = output from read.structure
##  nameList = list of indiduals by sp

str.rename <- function(x, nameList) {
  for(clustCol in 2:length(x)) {
    temp = sapply(spList, function(x) mean(fsIn[[1]][x, clustCol]))
    names(x)[clustCol] <- names(nameList)[which(temp == max(temp))]
  }
  names(x) <- make.unique(names(x))
  return(x)
}
