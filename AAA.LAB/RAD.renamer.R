## rename RAD-seq trees after you run them

require(magrittr)
require(ape)

dat.acer <-
  read.xlsx('T:/Systematics/DNA_DATABASES-LIVE/AAA.ACER_Specimen_Database.xlsx') %>% try
dat.tilia <-
  read.xlsx('T:/Systematics/DNA_DATABASES-LIVE/AAA.Tilia_Specimen_Database.xlsx') %>% try
dat.malus <-
  read.xlsx('T:/Systematics/DNA_DATABASES-LIVE/AAA.Malus_Specimen_Database.xlsx') %>% try

morTipLabels <- function(
    treesToDo,
    dat.meta,
    labelCols = c(
      "TAXA-Current_determination",
      "State", "Country","Collector_no",
      'RAD-SEQ-FLORAGENEX',
      'SPMCODE'),
    matchCol = 'RAD-SEQ-FLORAGENEX',
    grepCols = TRUE,
    tipRemovals = '',
    delim = '|',
    outfile = NA)
    {
  if('phylo' %in% class(treesToDo)) treesToDo <- list(treesToDo)
  if(length(intersect(class(treesToDo), c('phylo', 'list'))) < 1) {
    stop('treesToDo should be either a phylo object or a list of phylo objects')
    } # close if

  if(grepCols)
    labelCols <- sapply(labelCols, function(x) grep(x, names(dat.meta), value = T)[1])
  treesOut <- treesToDo
  for(i in 1:length(treesOut)) {
    matchTips <- match(tidyName(treesOut[[i]]$tip.label), tidyName(dat.meta[[matchCol]]))
    if(matchTips %>% is.na %>% any) warning(paste('Missing', sum(is.na(matchTips)), 'tips in metadata'))
    treesOut[[i]]$tip.label <-
      apply(dat.meta[matchTips, labelCols],
            1,
            paste, collapse = delim
          )
    rm(matchTips)
  } # close i
} # close function
