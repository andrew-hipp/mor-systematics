LG.plot <- function(lociBlast, markerPositions, max.evalue = NA, min.alignL = NA, lg.name = 'LG', pos.name = 'consensus', lg = NA, tickBounds = c(-0.1, 0.1), label = TRUE, tick.cex = 0.4, markDupes = c('left', 'right', 'n'), markOverlaps = TRUE, totals = TRUE, extraValues = NA, ...) {
## plots the linkage group locations of loci, given a b6 output of loci blasted to mapped markers, and map positions of the markers
## doesn't really belong in pyRAD, but it's the most convenient spot for it to live right now
## ARGUMENTS:
##  lociBlast - BLASTN of loci against mapped markers, b6 format (http://drive5.com/usearch/manual/blast6out.html)
##  markerPositions - a data.frame with map positions for the targets of the loci blast, with row.names being the marker names
##  lg.name - column name for the linkage groups indicated in markerPositions
##  pos.name - column name for the map position indicated in markerPostions
##  lg - vector of linkage group names to map

  names(lociBlast) <- c("query", "target", "idPercent", "alignL", "mismatchN", "gapN", 
                        "startPosQuery", "endPosQuery", "startPosTarget", "endPosTarget", 
                        "evalue", "bitscore")
  if(!is.na(max.evalue)) lociBlast <- lociBlast[lociBlast$evalue <= max.evalue, ]
  if(!is.na(min.alignL)) lociBlast <- lociBlast[lociBlast$alignL >= min.alignL, ]
  x <- as.data.frame(cbind(lociBlast, markerPositions[lociBlast$query, ]))
  if(is.na(lg[1])) lg <- sort(unique(x[[lg.name]]))
  lg.ranges <- t(sapply(lg, function(z) range(x[[pos.name]][(x[[lg.name]] == z)])))
  plot(1, xlim = c(0, length(lg) + 1), ylim = c(-1, max(lg.ranges) + 20), type = 'n', xaxt = 'n', ...)
  axis(1, at = seq(length(lg)), labels = lg, cex.axis = 0.6)
  for(i in 1:length(lg)) {
	segments(i, 0, i, lg.ranges[i, 2], ...)
	x.temp <- x[x[[lg.name]] == lg[i], ]
	segments(i + tickBounds[1], x.temp[[pos.name]], i+tickBounds[2])
	if(label) text(i + tickBounds[2], x.temp[[pos.name]], x.temp$query, pos = 4, cex = tick.cex)
	if(markDupes[1] != 'n') {
	  xpos = switch(markDupes[1], left = i + tickBounds[1] - 0.05, right = i + tickBounds[2] + 0.05)
	  x.dupes <- names(duplicated.mapped.loci(x))
	  for(j in 1:length(x.dupes)) {
	    ypos <- x.temp[[pos.name]][x.temp$query == x.dupes[j]]
		points(rep(xpos, length(ypos)), ypos, col = j)
        } # close for
	  } # close if
	if(markOverlaps) {
	  # browser()
	  x.overlaps <- x.temp[[pos.name]][which(duplicated(x.temp[[pos.name]]))]
	  for(j in unique(x.overlaps)) text(x = i + tickBounds[1], y = j, labels = sum(x.temp[[pos.name]] == j), cex = tick.cex, pos = 2, offset = 0.1)
	  } # close if
	if(totals) {
	  total.queries <- length(unique(x.temp$query))
	  total.targets <- length(unique(x.temp$target))
	  text(i, max(lg.ranges) + c(20,15,10), c(total.queries, total.targets), cex = 0.7)
	  } # close if
    } # close i
	if(totals) {
	  text(x = (length(lg) + 0.2), y = (max(lg.ranges) + c(20,15)), as.character(c(length(unique(x$query)), length(unique(x$target)))), cex = 0.7, pos = 4)
	  text(0.8, max(lg.ranges) + c(20,15,10), c('RAD loci', 'Contigs'), cex = 0.7, pos = 2)
	  print(length(unique(x$query)))
	  print(length(unique(x$target)))
	  }
  } # done