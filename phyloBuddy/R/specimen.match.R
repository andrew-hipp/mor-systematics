## specimens based on NCBI parsed output
## export specimen-linked datasets
## A Hipp, 2013-02-05
## version 2: 2013-06-12
## v3, 2015-03-26: gets rid of tidyName, makes this work directly on the entire parsed table from parse.INSDSeq

parse.specimen <- function(ncbiDat) {
## try to parse out the pieces of the specimen field
## xmlDat is exported from parse.INSDSeq
  require(gdata)
  fields <- c('All collectors', 'Primary collector last name', 'Collector number', 'Collection', 'Unedited text')
  obj <- ncbiDat[, 'specimen_voucher']
  ac <- pcln <- cn <- coll <- character(length(obj))
  obj.split <- strsplit(obj, " ")
  
  ## parse out collectors
  ac <- gsub("s.n.", "", trim(sapply(strsplit(obj, '[0123456789]'), function(x) x[1])), fixed = T)
  ac[grep(":", ac)] <- sapply(strsplit(ac[grep(":", ac)], ":"), function(x) x[2])
  
  ## swap in "collected_by" field when that is relevant
  use.collected.by <- which(apply(cbind(ncbiDat[, 'collected_by'], ac), 1, function(x) nchar(x[1]) > nchar(x[2])))
  ac[use.collected.by] <- ncbiDat[use.collected.by, 'collected_by']
  
  ## if there are no good collectors, use authors
  use.authors <- which(is.na(ac))
  ac[use.authors] <- paste('AUTHOR', ncbiDat[use.authors, 'authors'], sep = ":")
  
  ## get last name of primary collectors: Seems to be easiest just to take first word > 1 character long, but some editing is needed
  
  ac.splitted <- strsplit(ac, "[. ,]")
  for(i in which(!is.na(ac))) {
    ac.temp <- ac.splitted[[i]][nchar(ac.splitted[[i]]) > 1]
	pcln[i] <- ac.temp[1]
	}
  
  ## for collector number, assume that > 1/2 of characters are numbers
  countInstances <- function(x, pattern = as.character(0:9), proportion = TRUE) {
    x.split <- strsplit(x, NULL)[[1]]
	x.sum <- sum(x.split %in% pattern)
    if(proportion) out <- x.sum / nchar(x)
	else out <- x.sum
    return(out)
  }
  for(i in 1:length(obj)) {
 	if(obj[i] == "" | is.na(obj[i])) next
	temp <- sapply(obj.split[[i]], countInstances)
	cn[i] <- obj.split[[i]][which(temp == max(temp))][1]
	}
  ## make anything with "s.n." into a collector number of "s.n."
  sn <- unique(c(grep('s.n.', as.character(obj), fixed = T), grep('s. n.', as.character(obj), fixed = T)))
  cn[sn] <- 's.n.'

  ## for collection, assume either a "(" or all caps followed by ":"; exclude any with "s." in them
  colonDelimits <- grep(":", obj)
  coll[-colonDelimits] <- sapply(obj.split[-colonDelimits], function(x) {
                                   out <- gsub("(", "", gsub(")", "", grep("(", x, fixed = T, value = T), fixed = T), fixed = T)
								   if(length(out) == 0) return("")
								   else return(out[length(out)])
								   })
  coll[colonDelimits] <- sapply(obj.split[colonDelimits], function(x) strsplit(x, ":")[[1]][1])
  coll <- unlist(coll)
  
  out <- data.frame(ac = ac, pcln = pcln, cn = cn, coll = coll, dat = obj)
  names(out) <- fields
  out}
