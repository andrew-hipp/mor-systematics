## IOS analyses, Andrew Hipp
## 25 January 2013

## get data in, summarize
oaks.d6m4 <- read.pyRAD('./pyRAD_2013-01-23/c_hipp_d6m4p2_wRE.readloci.txt')
overlap.report.d6m4 <- overlap.report(oaks.d6m4$radSummary$inds.mat)
layout(matrix(1:9,3,3))
apply(overlap.report.d6m4[, c(3,5)], 1, function(x) pie(x, radius = sum(x)/32000, col = c('black', 'white'), labels =""))

## ordinate by shared loci; note that this and following analyses were not shown in the IOS paper
library(vegan)
oaks.d6m4.jaccard.MDS <- metaMDS(oaks.d6m4$radSummary$inds.mat[-which(row.names(oaks.d6m4$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard')

cols <- rep('red', dim(oaks.d6m4.jaccard.MDS$points)[1])
cols[grep("_re", row.names(oaks.d6m4.jaccard.MDS$points), fixed = TRUE)] <- 'black'
plot(oaks.d6m4.jaccard.MDS$points, pch = 19, col= cols, cex = 2)
text(oaks.d6m4.jaccard.MDS$points, row.names(oaks.d6m4.jaccard.MDS$points), cex = 0.7, adj = -0.2)

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
  
