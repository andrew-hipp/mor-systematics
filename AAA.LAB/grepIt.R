## grep for a vector
grepIt <- function(x, pattern, use = 1, as.char = T) {
  out <- sapply(x, function(i) grep(i, pattern, value = T)[use])
  if(as.char) out <- as.character(out)
  return(out)
}
