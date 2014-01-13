rad2phy <- function(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], outfile = 'pyMat.out.phy', padding = 50) {
## makes a phylip-style data matrix from rad.mat output, limiting by individuals and loci
  if(class(pyDat) != "rad.mat") warning("I'm expecting output from rad.mat")
  outfile = file(outfile, "wt")
  open(outfile)
  cat(paste(length(inds), sum(sapply(pyDat[inds[1], loci], nchar)), "\n"), file = outfile) #header: number of individuals, number of bases
  for(i in inds) {
    message(paste("Writing DNA line for individual", i))
	cat(i, file = outfile)
	cat(paste(rep(" ", padding - nchar(i)), collapse = ""), file = outfile)
	cat(paste(pyDat[i, loci], collapse = ""), file = outfile)
	cat("\n", file = outfile) # endline
	}
  close(outfile)
  }
 
write.DNAStringSet <- function(x, format= 'phylip', padding = 30, filename = 'DNAStringSetOut.phy') {
  # writes a sequence matrix to phylip format
  x.width <- width(x)[1]
  x <- as.character(x)
  for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
  writeLines(c(paste(length(x), x.width), x), filename)
  return(0)
  }
  
