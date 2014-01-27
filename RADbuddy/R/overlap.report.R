overlap.report <-
function(dat, repPattern = "_re", origPattern = "_h") {
## reports on how replicate individuals fall out in the loci
  if(class(dat) == 'summary.pyRAD.loci') dat <- dat$inds.mat
  reps <- grep(repPattern, row.names(dat), fixed = TRUE, value = TRUE)
  orig <- gsub(repPattern, origPattern, reps, fixed = TRUE)
  outCols <- c("Original", "Replicate", "Intersection", "Union", "Not intersection", "Overlap proportion", "Original loci replicated", "Increase, original to replicate")
  out <- matrix(NA, nrow = length(orig), ncol = length(outCols), dimnames = list(orig, outCols))
  for(i in 1:length(orig)) {
    message(paste("Doing names", orig[i]))
	out[orig[i], "Original"] <- sum(dat[orig[i],])
	out[orig[i], "Replicate"] <- sum(dat[reps[i],])
	out[orig[i], "Intersection"] <- sum(colSums(dat[c(orig[i], reps[i]), ]) == 2, na.rm = T)
	out[orig[i], "Union"] <- sum(colSums(dat[c(orig[i], reps[i]), ]) %in% 1:2, na.rm = T)
	}
  out[, "Not intersection"] <- out[, "Union"] - out[, "Intersection"]
  out[, "Overlap proportion"] <- round(out[, "Intersection"] / out[, "Union"], 3)
  out[, "Original loci replicated"] <- round(out[, "Intersection"] / out[, "Original"], 3)
  out[, "Increase, original to replicate"] <- round(out[, "Replicate"] / out[, "Original"], 3)
  return(out)
  }
