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


	