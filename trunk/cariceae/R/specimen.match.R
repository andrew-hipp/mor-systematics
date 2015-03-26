## specimens based on NCBI parsed output
## export specimen-linked datasets
## A Hipp, 2013-02-05
## version 2: 2013-06-12

## handy functions
nospace <- function(x) gsub(" ", "", x)
tidyNames <- function(x) gsub("[# -.,]", "_", x)

clean.specimen <- function(x) {
## cleans specimen names for known genbank oddities
   x <- gsub('                     /specimen_voucher=', '', gsub('\"', '', x, perl = T), fixed = T)
   x
   }
 
parse.specimen <- function(obj) {
## try to parse out the pieces of the specimen field
  require(gdata)
  fields <- c('All collectors', 'Primary collector last name', 'Collector number', 'Collection', 'Unedited text')
  ac <- pcln <- cn <- coll <- character(length(obj))
  obj.split <- strsplit(obj, " ")
  
  ## get collectors
  ac <- gsub("s.n.", "", trim(sapply(strsplit(obj, '[0123456789]'), function(x) x[1])), fixed = T)
  ac[grep(":", ac)] <- sapply(strsplit(ac[grep(":", ac)], ":"), function(x) x[2])
  
  ## would it make sense to take the pcln as the longest element in ac, after strspliting by " "? no... second collectors...
  
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
  browser()
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