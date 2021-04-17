quickPassword <- function(n = 16) {
  require(magrittr)
  paste(c(letters,LETTERS, 0:9, rep(c('#', '$', '%', '_'), 2))) %>% 
    sample(n) %>% 
    paste(collapse = '')
  }
