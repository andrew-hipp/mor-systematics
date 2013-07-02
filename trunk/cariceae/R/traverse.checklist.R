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
