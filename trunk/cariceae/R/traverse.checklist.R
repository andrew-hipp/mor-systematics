traverse.checklist <- function(spTaxonomy, topParent = "Cyperaceae", lowestRank = NA, tabs = "", tabChar = "\t", outfileName = "spTaxonomy.out.txt", includeGeog = F, wikiStyle = T, file.encoding = "UTF-8") {
##  Currently, you need to send it an accepted-names only list if you want an accepted-names only result
##  Arguments:
##    spTaxonomy = a Scratchpads checklist
##    topParent = node at which to begin building checklist
##    lowestRank = vector of ranks at which to stop (e.g., "SEQUENCE", c("Subspecies", "Variety"), or "Species"); use NA to go to lowest rank
  
  if(spTaxonomy$Rank[spTaxonomy$Term.name == topParent] %in% lowestRank) return(0) # doesn't recurse if we've hit the lowest rank we care about
  children <- which(gsub(" ", "", tolower(spTaxonomy$Parent.Term.Name)) == gsub(" ", "", tolower(topParent)))
  if(length(children) == 0) return(0) # doesn't recurse if there are no children
  outFile = file(outfileName, open = "a", encoding = file.encoding)
  nameCols <- c('Term.name','Authors') # orders name columns in the right order
  childrenFormatted <- paste(as.character(gsub("  {2, }", " ", apply(spTaxonomy[children, nameCols], 1, paste, collapse = " "))), sep = "")
  if(includeGeog) childrenAuxData <- as.character(paste('[', apply(spTaxonomy[children, c('Rank', 'Usage', 'Geography')], 1, paste, collapse = ", "), ']', sep = ''))
  else childrenAuxData <- as.character(paste('[', apply(spTaxonomy[children, c('Rank', 'Usage')], 1, paste, collapse = ", "), ']', sep = ''))
  sortOrder <- order(childrenFormatted)
  rank_name <- spTaxonomy[children, c('Rank')][sortOrder]
  childrenSorted <- childrenFormatted[sortOrder]
  auxDataSorted <- childrenAuxData[sortOrder]
  for(i in 1:length(childrenSorted)) {
	if(wikiStyle) {  
	  if(rank_name[i] == "Species") tabsPrint <- c("#","") else tabsPrint = c(tabs, tabs)	#this is just a work-around for wiki formatting
	  out <- paste(tabsPrint[1], "'''", childrenSorted[i], "''' ", auxDataSorted[i], tabsPrint[2], sep = "")
	  }
	else out <- paste(tabs, childrenSorted[i], " ", auxDataSorted[i], sep = "")
	if(outfileName == "screen") cat(c(out, "\n")) # to screen
	else cat(c(out, "\n"), file = outFile) # to file
	iAsParent <- spTaxonomy$Term.name[children[sortOrder]][i]
	# print(iAsParent) # for debugging only
	nextTabs <- paste(tabs, tabChar, sep = "")
	traverse.checklist(spTaxonomy, iAsParent, tabs = nextTabs, tabChar = tabChar, outfileName = outfileName, includeGeog = includeGeog, wikiStyle = wikiStyle)
	}
  close(outFile)
  }
