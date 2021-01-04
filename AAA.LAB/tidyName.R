tidyName <- function(x, fixes = c('_', '.', ' ')) {
  x = tolower(x)
  for(i in fixes) x <- (gsub(i, "", x, fixed = T))
  x
  }