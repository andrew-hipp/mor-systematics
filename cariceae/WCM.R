### WCM manipulation and display scripts
### ahipp@mortonarb.org
### 2012-06-13: created
### 2013-05-13: updating to update the checklist in a more streamlined manner

read.all.collector.sheets <- function(datDir, special.files = c('CHECKLIST.csv', 'GEO.csv', 'PermittedValues.csv', 'VersionNotes.csv'), rm.csv = T) {
## 13 May 2013 -- assumes all files have been exported using as CSV from standard collecting workbook
  files <- dir(datDir, full = T)
  names(files) <- dir(datDir, full = F)
  out <- files[!names(files) %in% special.files]
  out <- lapply(out, read.csv, as.is = TRUE)
  if(rm.csv) names(out) <- gsub('.csv', '', names(out), fixed = T)
  return(out)
  }
  
taxonByLab <- function(x = 'EXPORT.2013-05-13') {
## This function appends to the WCM matrix a column for each lab, indicating whether that lab has reported any material or sequences
## Arguments:
##  x = data directory to look in for all data collecting sheets and checklist
  collSheets <- read.all.collector.sheets(x)
  wcm <- read.csv(dir(x, full = T, patt = 'CHECKLIST'))
  headers <- names(collSheets)
  for(i in headers) {
    message('doing', i)
	# browser()
	wcm[[i]] <- (wcm$WCM.ID %in% collSheets[[i]]$ID) & (wcm$WCM.ID != "")
	}
  return(wcm)
  }

makeOneDatasheet <- function(datDir = 'EXPORT.2013-05-13') {
## ultimately, all extractions data should be in a single table
## this does that temporarily for our datasheets


voucherBySpecies <- function(datDir = 'EXPORT.2013-05-13') {
## this function:
##	(1) creates a matrix with taxonomic terms, parents, term authors, WCM id, and geography;
##	(2) expands that matrix with a row for each voucher, with Rank = 'VOUCHER', placing the taxon for that voucher as the parent of that voucher; note that this will have the odd side-effect of making the voucher of a species taxonomically equivalent to infraspeces of that species. Perhaps set off with a big character? e.g., -V- or ***; voucher field includes lab ownership / contact for extraction
##	(3) Adds sequences as children of each voucher, with Rank = 'SEQUENCE'; sequence fields include region name and genbank number
##	(4) optionally:
##		(a) calls traverseChecklist to create a hierarchical checklist (text file, with indents)
##		(b) calls assignSectionsAndClades to assign every taxon to a section and higher clade (spreadsheet)
##		(c) highlights taxa for which no samples are available
##		(d) aggregates vouchers for synonyms up to accepted names, holding the synonym under 'Deposited as'

## 1
  

  return(out)
  }
  
traverseChecklist <- function(wcm = wcv.matched.v6.8$WCM.out, topParent = "Cyperaceae", include = c("all", "accepted"), tabs = "", tabChar = "\t", outfileName = "wcm.out.txt", includeGeog = F, includeDNA = F, wikiStyle = T) {
  children <- which(gsub(" ", "", tolower(wcm$parent_name)) == gsub(" ", "", tolower(topParent)))
  if(length(children) == 0) return(0)
  outFile = file(outfileName, open = "a", encoding = "UTF-8")
  nameCols <- c(1,3,2,5,4) # orders name columns in the right order
  # childrenFormatted <- paste(tabs, as.character(gsub("  {2, }", " ", apply(wcm[children, nameCols], 1, paste, collapse = " "))), sep = "")
  childrenFormatted <- paste(as.character(gsub("  {2, }", " ", apply(wcm[children, nameCols], 1, paste, collapse = " "))), sep = "")
  if(includeGeog) childrenAuxData <- as.character(paste('[', apply(wcm[children, c('rank_name', 'usage', 'geo_cariceae_wcm')], 1, paste, collapse = ", "), ']', sep = ''))
  else childrenAuxData <- as.character(paste('[', apply(wcm[children, c('rank_name', 'usage')], 1, paste, collapse = ", "), ']', sep = ''))
  if(includeDNA) dnaIndicator <- paste(ifelse(wcm[children, 'onHand'], "[%] ", ""), ifelse(wcm[children, 'toGet'], "[@] ", ""), sep = "")

  sortOrder <- order(childrenFormatted)
  rank_name <- wcm[children, c('rank_name')][sortOrder]
  childrenSorted <- childrenFormatted[sortOrder]
  auxDataSorted <- childrenAuxData[sortOrder]
  if(includeDNA) dnaSorted <- dnaIndicator[sortOrder]
  for(i in 1:length(childrenSorted)) {
	if(wikiStyle) {  
	  if(rank_name[i] == "SPECIES") tabsPrint <- c("#","") else tabsPrint = c(tabs, tabs)	#this is just a work-around for wiki formatting
	  out <- ifelse(includeDNA, paste(tabsPrint[1], "'''", dnaSorted[i], childrenSorted[i], "''' ", auxDataSorted[i], tabsPrint[2], sep = ""), paste(tabsPrint[1], "'''", childrenSorted[i], "''' ", auxDataSorted[i], tabsPrint[2], sep = ""))
	  }
	else out <- ifelse(includeDNA, paste(tabs, dnaSorted[i], childrenSorted[i], " ", auxDataSorted[i], sep = ""), paste(tabs, childrenSorted[i], " ", auxDataSorted[i], sep = ""))
	if(outfileName == "screen") cat(c(out, "\n")) # to screen
	else cat(c(out, "\n"), file = outFile) # to file
	iAsParent <- gsub("[ \t]", "", childrenSorted[i])
	# print(iAsParent)
	nextTabs <- paste(tabs, tabChar, sep = "")
	traverseChecklist(wcm, iAsParent, tabs = nextTabs, tabChar = tabChar, outfileName = outfileName, includeGeog = includeGeog, includeDNA = includeDNA, wikiStyle = wikiStyle)
	}
  close(outFile)
  }

summaryByRegion <- function(wcm = scan(file.choose(), what = character(), sep = "\n", fileEncoding = "UTF-8"), tdwg = read.delim('TDWG.regions.level2', as.is = TRUE), sampleChar = "%", toGetChar = "@", exportByRegion = F, excludeHybrids = TRUE) {
  ## get wcm by reading in a summary file using wcm <- scan(file.choose(), what = character(), sep = "\n", fileEncoding = "UTF-8")
  out <- matrix(NA, nrow = dim(tdwg)[1] + 1, ncol = 5, dimnames = list(c(tdwg$label, "GLOBAL"), c('Total species', 'Species sampled', 'Species sampled or easily obtained', 'Species not accounted for', 'Percent')))
  wcm.spp <- wcm[grep("SPECIES", wcm, fixed = TRUE)]
  if(excludeHybrids) wcm.spp <- wcm.spp[-grep("×", wcm.spp, fixed = TRUE)]
  sppNum <- length(wcm.spp)
  sampled <- seq(sppNum) %in% grep(sampleChar, wcm.spp, fixed = T)
  toSample <- seq(sppNum) %in% grep(toGetChar, wcm.spp, fixed = T)
  sampledOrObtained <- sampled | toSample
  for(i in 1:dim(tdwg)[1]) {
    temp <- seq(sppNum) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm.spp, fixed = F, perl = F)
	out[i, ] <- c(sum(temp), sum(temp & sampled), sum(temp & sampledOrObtained), sum(temp) - sum(temp & sampled), paste(formatC((sum(temp & sampled) / sum(temp)) * 100, 1, format = "f"), "%", sep = ""))
    if(exportByRegion) {
	  outfile = file(tdwg$label[i], open = "w", encoding = "UTF-8")
	  # writeLines(wcm[(seq(length(wcm)) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm, fixed = F, perl = F)) | (!seq(length(wcm)) %in% grep("SPECIES", wcm, fixed = F, perl = F))], con = outfile)
	  writeLines(sort(wcm[seq(length(wcm)) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm, fixed = F, perl = F)]), con = outfile)
	  close(outfile)
	  }
	}
  out["GLOBAL", ] <- c(sppNum, sum(sampled), sum(sampledOrObtained), sppNum - sum(sampled), paste(formatC((sum(sampled) / sppNum) * 100, 1, format = "f"), "%", sep = ""))
  return(out) 
  }