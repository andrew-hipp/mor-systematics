traverse.checklist <- function(spTaxonomy, topParent = "Cyperaceae", lowestRank = NA, tabs = "", tabChar = "\t", DNA.flag = 'DNA VOUCHER', dna.missing.char = c('***','*','*'), show = c('all', 'missing', 'dna'), outfileName = "spTaxonomy.out.txt", includeGeog = T, wikiStyle = F, file.encoding = "UTF-8") {
##  Currently, you need to send it an accepted-names only list if you want an accepted-names only result
##  Arguments:
##    spTaxonomy = a Scratchpads checklist
##    topParent = node at which to begin building checklist
##    lowestRank = vector of ranks at which to stop (e.g., "SEQUENCE", c("Subspecies", "Variety"), or "Species"); use NA to go to lowest rank
##    DNA.flag = character string used to flag a DNA extraction. If !is.na(DNA.flag), then species and infraspecies are flagged with dna.missing.char if they have no DNA voucher
##    dna.missing.character = character used if DNA is missing from species[1], subspecies[2], or variety[3]   
  names(dna.missing.char) <- c('Species', 'Subspecies', 'Variety')
  topParentRank <- spTaxonomy$Rank[spTaxonomy$Term.name == topParent]
  if(length(topParentRank) > 1) {
    message(paste('There are', length(topParentRank), 'instances of taxon', topParent, '-- only using the first to determine rank of the parent'))
	topParentRank <- topParentRank[1]
	}
  if(topParentRank %in% lowestRank) return(0) # doesn't recurse if we've hit the lowest rank we care about
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
	  iAsParent <- spTaxonomy$Term.name[children[sortOrder]][i]
      iRank <- spTaxonomy$Rank[children[sortOrder]][i]
	  if(!is.na(DNA.flag)) {
        children.temp <- which(gsub(" ", "", tolower(spTaxonomy$Parent.Term.Name)) == gsub(" ", "", tolower(iAsParent)))
        ## following line to try to aggregate DNA samples of first-level infraspecies up to species
		children.of.children.temp <- which(gsub(" ", "", tolower(spTaxonomy$Parent.Term.Name)) %in% gsub(" ", "", tolower(spTaxonomy$Term.name[children.temp])))
 	    if(iRank %in% c('Species', 'Subspecies', 'Variety')) {
	      if(length(grep(DNA.flag, spTaxonomy$Term.name[c(children.temp, children.of.children.temp)], fixed = TRUE)) == 0) {
		    if(show[1] %in% c('all', 'missing')) childrenSorted[i] <- paste(dna.missing.char[as.character(iRank)], childrenSorted[i])
			else return(0)
			}
		  }}
      if(wikiStyle) {  
	    if(rank_name[i] == "Species") tabsPrint <- c("#","") else tabsPrint = c(tabs, tabs)	#this is just a work-around for wiki formatting
	    out <- paste(tabsPrint[1], "'''", childrenSorted[i], "''' ", auxDataSorted[i], tabsPrint[2], sep = "")
	    }
	  else out <- paste(tabs, childrenSorted[i], " ", auxDataSorted[i], sep = "")
	  if(outfileName == "screen") cat(c(out, "\n")) # to screen
	  else cat(c(out, "\n"), file = outFile) # to file
	  # print(iAsParent) # for debugging only
	  nextTabs <- paste(tabs, tabChar, sep = "")
	  traverse.checklist(spTaxonomy, iAsParent, tabs = nextTabs, tabChar = tabChar, outfileName = outfileName, includeGeog = includeGeog, wikiStyle = wikiStyle)
	  }
    close(outFile)
    return(outfileName)
    }
