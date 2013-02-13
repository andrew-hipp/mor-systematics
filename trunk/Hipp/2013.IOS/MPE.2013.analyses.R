## MPE analyses, Andrew Hipp
## 13 February 2013

r2pyrad <- 'http://mor-systematics.googlecode.com/svn/trunk/AAA.Lab/r2pyRAD.R'
source(r2pyrad)


## get data in, summarize
d6m10p2.wRE <- read.pyRAD('../eatonAnalyses/c88_d6m10p2_wRE.readloci.txt')
d6m4p2.wRE <- read.pyRAD('../eatonAnalyses/c88_d6m4p2_wRE.readloci.txt')

overlap.report.d6m10p2.wRE <- overlap.report(d6m10p2.wRE$radSummary$inds.mat)
layout(matrix(1:9,3,3))
apply(overlap.report.d6m10p2.wRE[, c(3,5)], 1, function(x) pie(x, radius = sum(x)/32000, col = c('black', 'white'), labels =""))

overlap.report.d6m4p2.wRE <- overlap.report(d6m4p2.wRE$radSummary$inds.mat)
layout(matrix(1:9,3,3))
apply(overlap.report.d6m4p2.wRE[, c(3,5)], 1, function(x) pie(x, radius = sum(x)/40000, col = c('black', 'white'), labels =""))

## ordinate by shared loci; note that this and following analyses were not shown in the IOS paper
library(vegan)
d6m10p2.wRE.jaccard.MDS.stresses <- sapply(1:10, function(x) metaMDS(d6m10p2.wRE$radSummary$inds.mat[-which(row.names(d6m10p2.wRE$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard', k = x, autotransform = F, wasscores = F, noshare = F, trymax = 50)$stress)
names(d6m10p2.wRE.jaccard.MDS.stresses) <- paste("K=",1:10, sep = "")
plot(1:10, d6m10p2.wRE.jaccard.MDS.stresses, pch = 19); lines(1:10, d6m10p2.wRE.jaccard.MDS.stresses)
d6m10p2.wRE.jaccard.MDS.k2 <- metaMDS(d6m10p2.wRE$radSummary$inds.mat[-which(row.names(d6m10p2.wRE$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard', k = 2, autotransform = F, wasscores = F, noshare = F, trymax = 2000)
d6m10p2.wRE.jaccard.MDS.k3 <- metaMDS(d6m10p2.wRE$radSummary$inds.mat[-which(row.names(d6m10p2.wRE$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard', k = 3, autotransform = F, wasscores = F, noshare = F, trymax = 2000)

d6m4p2.wRE.jaccard.MDS.k2 <- metaMDS(d6m4p2.wRE$radSummary$inds.mat[-which(row.names(d6m4p2.wRE$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard', k = 2, autotransform = F, wasscores = F, noshare = F, trymax = 2000)
d6m4p2.wRE.jaccard.MDS.k3 <- metaMDS(d6m4p2.wRE$radSummary$inds.mat[-which(row.names(d6m4p2.wRE$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard', k = 3, autotransform = F, wasscores = F, noshare = F, trymax = 2000)


layout(matrix(1:4,2,2))
## plot k = 2 solution
cols <- rep('red', dim(d6m10p2.wRE.jaccard.MDS.k2$points)[1])
cols[grep("_re", row.names(d6m10p2.wRE.jaccard.MDS.k2$points), fixed = TRUE)] <- 'black'
plot(d6m10p2.wRE.jaccard.MDS.k2$points, pch = 19, col= cols, cex = 2, main = 'd6m10p2.wRE')
text(d6m10p2.wRE.jaccard.MDS.k2$points, row.names(d6m10p2.wRE.jaccard.MDS.k2$points), cex = 0.7, adj = -0.2)

## plot k = 3 solution
cols <- rep('red', dim(d6m10p2.wRE.jaccard.MDS.k3$points)[1])
cols[grep("_re", row.names(d6m10p2.wRE.jaccard.MDS.k3$points), fixed = TRUE)] <- 'black'
limits  = apply(d6m10p2.wRE.jaccard.MDS.k3$points, 2, range)
scatterplot3d(d6m10p2.wRE.jaccard.MDS.k3$points, pch = 19, color = cols, cex.symbols = 2, type = 'h', xlim = limits[, 1], ylim = limits[, 2], zlim = limits[, 3], highlight.3d = F, main = 'd6m10p2.wRE, K = 3')

## plot k = 2 solution
cols <- rep('red', dim(d6m4p2.wRE.jaccard.MDS.k2$points)[1])
cols[grep("_re", row.names(d6m4p2.wRE.jaccard.MDS.k2$points), fixed = TRUE)] <- 'black'
plot(d6m4p2.wRE.jaccard.MDS.k2$points, pch = 19, col= cols, cex = 2, main = 'd6m4p2.wRE')
text(d6m4p2.wRE.jaccard.MDS.k2$points, row.names(d6m4p2.wRE.jaccard.MDS.k2$points), cex = 0.7, adj = -0.2)

## plot k = 3 solution
cols <- rep('red', dim(d6m40p2.wRE.jaccard.MDS.k3$points)[1])
cols[grep("_re", row.names(d6m4p2.wRE.jaccard.MDS.k3$points), fixed = TRUE)] <- 'black'
limits  = apply(d6m4p2.wRE.jaccard.MDS.k3$points, 2, range)
scatterplot3d(d6m4p2.wRE.jaccard.MDS.k3$points, pch = 19, color = cols, cex.symbols = 2, type = 'h', xlim = limits[, 1], ylim = limits[, 2], zlim = limits[, 3], highlight.3d = F, main = 'd6m4p2.wRE, K = 3')

########################################
## plots below reported on, not included
########################################
d6m10p2.wRE.jaccard.MDS.k2.withOddballs <- metaMDS(d6m10p2.wRE$radSummary$inds.mat, 'jaccard', k = 2, autotransform = F, wasscores = F, noshare = F, trymax = 2000)
d6m10p2.wRE.jaccard.MDS.k3.withOddballs <- metaMDS(d6m10p2.wRE$radSummary$inds.mat, 'jaccard', k = 3, autotransform = F, wasscores = F, noshare = F, trymax = 2000)

layout(matrix(1:2,1,2))
## plot k = 2 with-oddballs solution
cols <- rep('red', dim(d6m10p2.wRE.jaccard.MDS.k2.withOddballs$points)[1])
cols[grep("_re", row.names(d6m10p2.wRE.jaccard.MDS.k2.withOddballs$points), fixed = TRUE)] <- 'black'
plot(d6m10p2.wRE.jaccard.MDS.k2.withOddballs$points, pch = 19, col= cols, cex = 2)
text(d6m10p2.wRE.jaccard.MDS.k2.withOddballs$points, row.names(d6m10p2.wRE.jaccard.MDS.k2.withOddballs$points), cex = 0.7, adj = -0.2)

## plot k = 3 with-oddballs solution
cols <- rep('red', dim(d6m10p2.wRE.jaccard.MDS.k3.withOddballs$points)[1])
cols[grep("_re", row.names(d6m10p2.wRE.jaccard.MDS.k3.withOddballs$points), fixed = TRUE)] <- 'black'
limits  = apply(d6m10p2.wRE.jaccard.MDS.k3.withOddballs$points, 2, range)
scatterplot3d(d6m10p2.wRE.jaccard.MDS.k3.withOddballs$points, pch = 19, color = cols, cex.symbols = 2, type = 'h', xlim = limits[, 1], ylim = limits[, 2], zlim = limits[, 3], highlight.3d = F)


## read in and concatenate blast results to EST databases
blast.results <- lapply(dir('../../../OAK-RAD-ANNOTATIONS/QuercusFeb0513/', full = T, patt = 'blastN'), read.delim, header = F, as.is = T)
names(blast.results)<- (dir('../../../OAK-RAD-ANNOTATIONS/QuercusFeb0513/', full = F, patt = 'blastN')
blast.results.concat <- blast.results[[1]]
for(i in 2:length(blast.results)) blast.results.concat <- rbind(blast.results.concat, blast.results[[i]])

loci.by.threshold <- function(x = blast.results.concat, threshold = -15) {
  out = sort(unique(x[[1]][x[[4]]< (1 * 10^threshold)]))
  return(out)
  }

## Making a data matrix for the EST-linked markers
oaks.d6m4.mat <- rad2mat(oaks.d6m4)
rad2phy(oaks.d6m4.mat, loci = loci.by.threshold(), outfile = 'oaks.dna.d6m4.blast.e-15.phy')
rad2phy(oaks.d6m4.mat, loci = loci.by.threshold(threshold = -25), outfile = 'oaks.dna.d6m4.blast.e-25.phy')

## Subsampling loci to see how phylogeny based on non-EST linked markers compares
## PROBLEMATIC: mean length of the EST-linked markers is longer, so this results in a shorter data matrix
subsampled.loci <- lapply(rep(8833,10), sample, x = dimnames(oaks.d6m4.mat)[[2]][!dimnames(oaks.d6m4.mat)[[2]] %in% blast.results.concat[[1]]])
for(i in 1:10) rad2phy(oaks.d6m4.mat, loci = subsampled.loci[[i]], outfile = paste('oaks.dna.d6m4.blast.subsample.', i,'.phy', sep = ''))

do.EST.phylo <- function(dat = oaks.d6m4.mat, lengths = oaks.d6m4.lengths, blast = blast.results.concat, t.hold = -15, minLength = 85, maxLength = 95, sampleReps = 10) {
  if(is.na(lengths)) lengths <- lengths.report(dat,0) # as written, this will fail if lengths = NA, b/c lengths.report needs a pyRAD object, whereas this function needs a pyRAD.mat object
  sizeFilteredLoci <- names(lengths[lengths >= minLength & lengths <= maxLength])
  estLoci <- loci.by.threshold(blast, t.hold)
  estLoci <- estLoci[which(estLoci %in% sizeFilteredLoci)]
  otherLoci <- sizeFilteredLoci[!sizeFilteredLoci %in% unique(blast[[1]])]
  rad2phy(dat, loci = estLoci, outfile = paste('oak.ests.tHold.', t.hold, '.m', maxLength, '.', minLength, '.', format(Sys.time(), "%Y-%d-%m"), '.phy', sep = ''))
  for(i in 1:sampleReps) {
    message(paste("*** writing rep", i))
	sampledLoci <- sample(otherLoci, length(estLoci))
	rad2phy(dat, loci = sampledLoci, outfile = paste('oak.otherLoci.tHold.', t.hold, '.m', maxLength, '.', minLength, '.rep', i, '.', format(Sys.time(), "%Y-%d-%m"), '.phy', sep = ''))
	}
  return(0)
  }
  
