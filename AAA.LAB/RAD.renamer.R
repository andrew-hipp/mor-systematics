## rename RAD-seq trees after you run them

library(magrittr)
library(ape)

dat.metaFile = 'https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/AAA.DATA/RADSpecimensFloragenexLiveCopy_2020-11-16.csv'
matchCol = 'Extraction_Tube_no'
labelCols = c(
  "Species",
  "Inst.",
  "Extraction_Tube_no"
  )
delim = '|'
  # "Collector_no", "Wild_Garden_Collected",
  # "Original_Source_Information", "Vouchers", "MOR_Herbarium_Accession_Number",
  # "Notes", "Voucher_Image", "Photo_1", "Photo_2", "City", "County",
  # "State_Province", "Country", "Lat", "Long", "Georeference_Source",
  # "Cronn_CP", "AFLP", "Section", "Floragenex.Library.Name", "SEQUENCED.ON",
  # "SEQUENCING.DETAILS")

if(!exists('treesToDo')) {
  stop('Assign the tree(s) you want to relabel as elements of a list in treesToDo and rerun')
  }

dat.meta <- read.csv(dat.metaFile)
treesOut <- treesToDo
for(i in 1:length(treesOut)) {
  matchTips <- match(toupper(treesOut[[i]]$tip.label), toupper(dat.meta[[matchCol]]))
  if(matchTips %>% is.na %>% any) warning(paste('Missing', sum(is.na(matchTips)), 'tips in metadata'))
  treesOut[[i]]$tip.label <-
    apply(dat.meta[matchTips, labelCols],
          1,
          paste, collapse = delim
        )
  rm(matchTips)
}
