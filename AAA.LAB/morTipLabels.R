## rename RAD-seq trees after you run them

require(magrittr)
require(ape)
require(openxlsx)
if(!exists('tidyName')) source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/AAA.LAB/tidyName.R')

if(!exists('dat.acer')) dat.acer <-
  read.xlsx('T:/Systematics/DNA_DATABASES-LIVE/AAA.ACER_Specimen_Database.xlsx') %>% try
if(!exists('dat.tilia')) dat.tilia <-
  read.xlsx('T:/Systematics/DNA_DATABASES-LIVE/AAA.Tilia_Specimen_Database.xlsx') %>% try
if(!exists('dat.malus')) dat.malus <-
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
    tipRemovals = '.fq.gzbarcodeStripped',
    outgroup = NA,
    outgroupGrep = TRUE,
    delim = '|',
    ladder = TRUE,
    isNameVect = FALSE,
    outfile = NA,
    version = 1,
    pdfTree = TRUE,
    tip.cex = 0.6,
    nodeLabels = TRUE,
    pdfW = 10, pdfH = 20)
    {
  if(isNameVect) {
    treesToDo <- list(tip.label = treesToDo)
    treesToDo <- list(treesToDo)
  } else {
    if('phylo' %in% class(treesToDo)) treesToDo <- list(treesToDo)
    if(length(intersect(class(treesToDo), c('phylo', 'list'))) < 1) {
      stop('treesToDo should be either a phylo object or a list of phylo objects')
      } # close if
    } # close else

  if(grepCols)
    labelCols <- sapply(labelCols, function(x) grep(x, names(dat.meta), value = T)[1])
  treesOut <- treesToDo
  for(i in 1:length(treesOut)) {
    treesOut[[i]]$tip.label <-
      gsub(tipRemovals, '', treesOut[[i]]$tip.label)
    matchTips <- match(tidyName(treesOut[[i]]$tip.label), tidyName(dat.meta[[matchCol]]))
    if(matchTips %>% is.na %>% any) warning(paste('Missing', sum(is.na(matchTips)), 'tips in metadata'))
    treesOut[[i]]$tip.label <-
      apply(dat.meta[matchTips, labelCols],
            1,
            paste, collapse = delim
          )
    rm(matchTips)
    if(!isNameVect) {
      if(ladder) treesOut[[i]] <- ladderize(treesOut[[i]])
      if(outgroupGrep & !is.na(outgroup[1])) {
        outgroup = grep(outgroup, treesOut[[1]]$tip.label, value = T) %>%
          as.character
        } # close outgroup
        if(!is.na(outgroup[1])) {
          treesOut[[i]] <- try(root(treesOut[[i]], outgroup))
          message('Rooting by this outgroup:')
          message(outgroup)}
        }

  } # close i
  if(!is.na(outfile) & !isNameVect) {
    for(i in 1:length(treesOut)) {
      write.tree(treesOut[[i]], paste(outfile, '.', i, '_v', version, '.tre', sep = ''))
      if(pdfTree) {
        pdf(paste(outfile, '.', i, '_v', version, '.pdf', sep = ''), pdfW, pdfH)
        plot(treesOut[[i]], cex = tip.cex, show.node.label = nodeLabels)
        dev.off()
      }
    } # close write.tree
  } # close if(!is.na(outfile))
  return(treesOut)
} # close function
