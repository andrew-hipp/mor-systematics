## for linking RADs to Durand linkage map

loadDat <- FALSE

if(loadDat) durand <- read.delim('durand-2010')

seqsFromDurand <- function(dat = durand[durand$linkage.group.assignement. != "NO", ]) {
  # 1. make a table with clone as the row.name, then SSR.ID as the field
  # 2. search all the clones and export those as a fasta file
  # 3. blast against the clones fasta file, export b6 file
  # 4. read b6 file back in, get the linkage group assignment for each marker
  # 5. hist of linkage groups by individual
  # 6. filter phylogeny matrix by these loci, si possible.
  
  # 1
  cloneMat = matrix(NA,0,1,dimnames = list(NULL, "SSR.ID"))
  for(i in 1:dim(dat)[1]) {
    clones = sapply(strsplit(unlist(strsplit(as.character(dat$EST.ID[i]), ",")), ".", T), function(x) gsub(" ", "", x[1]))
    cloneMat = rbind(cloneMat, matrix(dat$SSR.ID[i], length(clones), 1, dimnames = list(clones, "SSR.ID")))
	} # not robust... will fail if any clone names are non-unique
	
  # 2
  fastaFile <- entrez(row.names(cloneMat))
  write.dna(fastaFile, format = "fasta", file = "trial.fas")
  }

linkBack <- function(b6.files = durand.b6.files, fileNames = durand.b6.file.names, refDat = durand[durand$linkage.group.assignement. != "NO", ], threshold = -14) {
  inDat <- lapply(b6.files, read.delim, header = F, colClasses = "character")
  names(inDat) <- fileNames
  outDat <- sapply(inDat, function(x) apply(x, 1, function(y) refDat[grep(y[[2]], refDat[, 'EST.ID']), 'linkage.group.assignement.'][1]))
  for(i in names(inDat)) outDat[[i]] <- cbind(inDat[[i]], outDat[[i]])
  if(!is.na(threshold)) outDat <- list(durand.mats = outDat, linkage.mats = linkMats(outDat, threshold))
  return(outDat)
  }

linkMats <- function(durand.dat, threshold, lgNo = 12) {
  #takes linkBack$durand.mats and make matrix of linkage groups for each one
  if(class(durand.dat) != "list") durand.dat <- list(durand.dat)
  outMat <- matrix(NA, nrow = length(durand.dat), ncol = lgNo, dimnames = list(names(durand.dat), paste("LG", seq(lgNo), sep = "")))
  for(i in 1:length(durand.dat)) {
    dat <- split(durand.dat[[i]], durand.dat[[i]][[1]]) # splits into matrices with all rows of a given name
	dat <- dat[log10(sapply(dat, function(x) min(as.numeric(x[, 11])))) <= threshold] # filters out only those matrices with a fit better than or equal to threshold
	lg.vector <- sapply(dat, function(x) as.integer(as.character(x[x[, 12] == min(x[, 12]), 13][1])))
	outMat[i, ] <- sapply(1:12, function(x) sum(lg.vector == x))
	}
  return(outMat)
  }
	