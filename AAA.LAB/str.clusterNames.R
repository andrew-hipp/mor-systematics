## name clusters in structure output
## args:
##  x = output from read.structure
##  nameList = list of indiduals by sp

str.rename <- function(x, nameList) {
  for(clustCol in 2:length(x)) {
    temp = sapply(spList, function(x) sum(fsIn[[1]][x, 2]))
    names(x)[clustCol] <- names(nameList)[which(temp == max(temp))]
  }
  return(x)
}
