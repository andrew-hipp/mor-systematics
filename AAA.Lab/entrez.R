## Functions for calling down data from Entrez
## Andrew Hipp, March 2011

## http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucest&usehistory=y&term=LG0AAA6YD14RM1
## then:
## http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucest&query_key=1&WebEnv=NCID_1_80041556_130.14.18.47_9001_1302193105_1430447178&retmode=fasta
## Longer example:
## http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucest&usehistory=y&term=LG0AAD22YM14RM1,LG0AAB7YM23RM1,LG0AAB2YK05RM1,LG0AAB9YK17RM1,LG0AAB6YF09RM1,LG0AAD34YN16RM1,LG0AAB17YE03RM1,LG0AAD36YO24RM1,LG0AAB7YM23RM1,LG0AAD31YB09RM1,LG0AAA8YJ20RM1,LG0AAB21YG10RM1,LG0AAA9YB15RM1,LG0AAD31YJ06RM1,LG0AAB19YG07RM1,LG0AAC6YD22RM1,LG0AAB12YC01RM1
## then:
##http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucest&query_key=1&WebEnv=NCID_1_80231330_130.14.18.47_9001_1302197555_145581083&retmode=fasta
## but this returns nothing... odd, b/c the esearch query pulled down info

loadDat <- FALSE

if(loadDat) durand <- read.delim('durand-2010')
# note: there are two rows duplicated in the durand table, POR047 and VIT064; neither is assigned to a linkage group

entrez <- function(queryVector) {
  require(XML)
  require(ape)

  urlBase <- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
  fetch.url <- '/efetch.fcgi?'
  search.url <- '/esearch.fcgi?'
  database <- 'nucest' # we are pulling data from the nuclear EST database of NCBI
  retmode <- 'fasta' # because we want our data back as fasta files
  fastaOut <- ncbiOut <- vector('list', 0)
  for(queryString in queryVector) {
    message(paste('Working on query', queryString))
    searchDat <- xmlTreeParse(paste(urlBase, search.url, 'db=', database, '&usehistory=y&term=', queryString, sep = ""))
    if(try(xmlValue(searchDat$doc$children$eSearchResult$children$Count) != "1")) { # turned this into a try b/c I got an error that I couldn't replicate
      message(paste(xmlValue(searchDat$doc$children$eSearchResult$children$Count), 'values returned for', queryString, '-- skipping on to the next!'))
      next
      }
    WebEnv <- xmlValue(searchDat$doc$children$eSearchResult$children$WebEnv)
    query_key <- xmlValue(searchDat$doc$children$eSearchResult$children$QueryKey)
    fetchDat <- read.dna(paste(urlBase, fetch.url, 'db=', database, '&query_key=', query_key, '&WebEnv=', WebEnv, '&retmode=', retmode, sep = ""), format = 'fasta')
    if(identical(fastaOut, NULL)) {
      fastaOut <- fetchDat
      names(fastaOut)[1] <- queryString
      } # close if
    else {
      fastaOut <- c(fastaOut, list(fetchDat))
      names(fastaOut)[length(fastaOut)] <- queryString
      } # close else
    } # close queryString loop
  return(fastaOut)
  }  
## this returns a DNAbin object, and these can be rbind'd together

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
	
	