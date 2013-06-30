## read-write-check functions
write.wcm <- function(wcm, fileBase, timeStamp = "generate", fileEnd = '.txt') {
## write wcm tables in the format that we want them in
  if(timeStamp == "generate") fileName = paste(fileBase, "_", format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), fileEnd, sep = "")
  else fileName = paste(fileBase, "_", timeStamp, fileEnd, sep = "")
  write.table(wcm, fileName, fileEncoding = "ANSI", quote = F, row.names = F, sep = "\t") # changed to ANSI in May 2013, b/c no longer working with WCM export directly
  }
 
read.wcm <- function(fileName = file.choose()) {
## read wcm tables in the format that I want them in; default to interactive file selection
  out <- read.delim(fileName, encoding = "ANSI", as.is = T) # changed to ANSI in May 2013, b/c no longer working with WCM export directly
  return(out)
  }

check.wcm <- function(WCM) {
## checks WCM for discrepancies between parents and names
  nameCols <- c(1,3,2,5,4)
  WCM.names <- gsub(" ", "", as.character(apply(WCM[, nameCols], 1, paste, collapse = "")))
  WCM.parents <- gsub(" ", "", WCM$parent_name)
  missingParents <- which(!WCM.parents %in% WCM.names)
  return(missingParents)
  }

