make.binary <- function(dat, binChars = strsplit('ACDEFGHIKLMNPQRSTVWY','')[[1]], extras = c('-')) {
  ## dat is a character matrix
  dat <- dat[, apply(dat, 2, function(x) all(x %in% c(binChars, extras)))] # keeps only columns with the binChars and extras
  charCount <- apply(dat, 2, function(x) sum(unique(x) %in% binChars)) # counts unique binChars
  dat <- dat[, charCount == 2]
  return(dat)
  }
  