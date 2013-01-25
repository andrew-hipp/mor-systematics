# These are things I did on 3 feb 2011 that produced a sensible lnL plot (trees-best-worst-by-loci_3feb11.pdf)

rtree.200.likelihoods.best <- apply(rtree.200.likelihoods$locusScores, 2, function(x) which(x == min(x)))
## note, rtree.200.likelihoods is a 'swulLikelihoods' object

rtree.200.likelihoods.best.byLocus <- numeric(201)
for(i in 1:201) {
rtree.200.likelihoods.best.byLocus[i] <- sum(rtree.200.likelihoods.best == i) }
