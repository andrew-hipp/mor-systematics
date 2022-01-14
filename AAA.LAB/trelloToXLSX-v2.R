require(trelloR)
require(magrittr)
require(openxlsx)

readData <- TRUE
writeAll <- TRUE
writeMasterList <- TRUE
addMonths  <- TRUE
includeTasks <- TRUE

#header <- c('2020 goals', 'Herbarium and Systematics', '--------------------')
header <- paste('Goal', 'Project', 'Objective', 'Task', 'DONE',
            'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC',
            sep = '\t')
if(readData) {
syst <- get_id_board('https://trello.com/b/xkR9QfSm/systematics-2020')

tr <- list(
    lists = get_board_lists(syst) %>% as.data.frame,
    cards = get_board_cards(syst) %>% as.data.frame)
## for check items, I think you need to grab the card ID AND the checkItems ID
## working example: https://api.trello.com/1/cards/5defca70a9ad294304bf47ee/checkItem/5defcc6bf697458a1e541fcf
row.names(tr$lists) <- tr$lists$id
row.names(tr$cards) <- tr$cards$id
tr$checklists <- lapply(tr$cards$id, function(x) get_card_checklists(x) %>% as.data.frame) # checkItems is currently the comma-delimited vector of things to do
#  names(tr$checklists) <- tr$cards$id
    tr$checklists <- do.call('rbind', tr$checklists)
tr$checkItems <- lapply(tr$checklists$id, function(x) {
    trello_get(id = x, parent = 'checklists', child='checkItems') %>%
      as.data.frame
    }
  )
#  names(tr$checkItems) <- do.call(rbind, tr$checklists)$id
  tr$checkItems <- tr$checkItems[sapply(tr$checkItems, dim)[1, ] > 0]
  tr$checkItems <- lapply(tr$checkItems, function(x) {x$nameData <- NULL; x})
  tr$checkItems <- do.call('rbind', tr$checkItems)
row.names(tr$checklists) <- tr$checklists$id
row.names(tr$checkItems) <- tr$checkItems$id
} # close if readData

if(writeAll)
    lapply(names(tr), function(i)
        write.xlsx(tr[[i]], paste(i, '.xlsx', sep = ''))
        )


if(writeMasterList) {
  master <- header
  for(list in tr$lists$id) {
    master <- c(master, paste('GOAL:', tr$lists[list, 'name']))
    for(card in tr$cards$id[which(tr$cards$idList == list)]) {
      master <- c(master, paste('\tPROJECT: ', tr$cards[card, 'name'], sep = ''))
      for(checklist in tr$checklists$id[which(tr$checklists$idCard == card)]) {
          master <- c(master, paste('\t\tOBJECTIVE: ', tr$checklists[checklist, 'name'], sep = ''))
          for(checkItem in tr$checkItems$id[which(tr$checkItems$idChecklist == checklist)]){
              if(includeTasks) master <- c(master, paste('\t\t\tTASK: ', tr$checkItems[checkItem, 'name'], sep = ''))
          } # close checkItem
      } # close for checklist
    } # close for card
  } # close for list
writeLines(master, 'masterList.tabDelim.tsv')
if(addMonths) {
  #master.df <- read.delim('masterList.tabDelim.tsv')
  master.df <- read.delim(text = master)
  master.df[is.na(master.df)] <- ''
  ml <- grep('|', master, fixed = T) - 1 # subtract 1 b/c master has header line, master.df has header
  m <- apply(master.df[ml,], 1, paste, collapse = '') %>% as.character
  m.split <- strsplit(m, '|', fixed = T)
  m.months <- sapply(m.split, '[', 2)
  m.months <- lapply(m.months, function(x) {
    temp <- strsplit(x, ' ', fixed = TRUE)[[1]]
    temp <- temp[temp != '']
  })
  #for(i in seq(length(ml))) master.df[ml[i], m.months[[i]]] <- m.months[[i]]
  for(i in seq(length(ml))) master.df[ml[i], m.months[[i]]] <- rep('XX', length(m.months[[i]]))
  write.table(master.df,
    paste('masterList.',
          ifelse(includeTasks, 'tasks.', ''),
          ifelse(addMonths, 'months.', ''),
          'tsv', sep = ''),
    row.names = FALSE, sep = '\t', quote = F)
  writeLines(grep('|', master, fixed= T, value = T, invert = T), 'missingMonths.txt')
}

} # close if write masterList
