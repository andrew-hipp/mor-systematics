## convert raxml-style long-format pairwise distance table into a square distance matrix

rax2dist <- function(file, diag.fill = 0, write = TRUE) {
  dat <- read.table(file, as.is = T)
  allTips <- sort(unique(unlist(dat[1:2])))
  out <- matrix(NA, length(allTips), length(allTips),
                dimnames = list(allTips, allTips))
  for(i in 1:dim(dat)[1]) out[dat[i, 'V1'], dat[i, 'V2']] <-
                          out[dat[i, 'V2'], dat[i, 'V1']] <-
                          dat[i, 'V3']
  diag(out) <- diag.fill
  if(write) write.csv(out, paste(file, '.dist.csv', sep = ''))
  out
  }
