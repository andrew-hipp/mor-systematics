## MPE analyses, Andrew Hipp
## 13 February 2013

r2pyrad <- 'http://mor-systematics.googlecode.com/svn/trunk/AAA.Lab/r2pyRAD.R'
source(r2pyrad)


## get data in, summarize
d6m10p2.wRE <- read.pyRAD('../eatonAnalyses/c88_d6m10p2_wRE.readloci.txt')
d6m4p2.wRE <- read.pyRAD('../eatonAnalyses/c88_d6m4p2_wRE.readloci.txt')
d6m10p2.wRE.mat <- rad2mat(d6m10p2.wRE)

## get blast results in, concatenate matrices
blasts.2013.02.21 <- dir('./blasts.2013-02-21/', full = T)
blast.results <- lapply(blasts.2013.02.21, read.delim, header = F, as.is = T)
names(blast.results)<- dir('./blasts.2013-02-21/', full = F)
blast.results.concat <- blast.results[[1]]
for(i in 2:length(blast.results)) blast.results.concat <- rbind(blast.results.concat, blast.results[[i]])
blast.results.concat$log.eValue <- log10(as.numeric(blast.results.concat[[4]]))
blast.results.concat <- blast.results.concat[-which(is.na(blast.results.concat$log.eValue)), ]

## do some summaries
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

loci.by.threshold <- function(x = blast.results.concat, threshold = -15) {
  out = sort(unique(x[[1]][x$log.eValue < threshold]))
  return(out)
  }

make.tree.distMat <- function(treelist, ...) {
  require(phangorn)
  outmat.sd <- outmat.bsd <- outmat.pd <- outmat.qpd <- outmat.RF <- matrix(NA, length(treelist), length(treelist), dimnames = list(names(treelist), names(treelist)))
  for(i in 2:length(treelist)) {
    for(j in 1:(i-1)) {
	  temp <- phangorn::treedist(treelist[[i]], treelist[[j]]) ## to differentiate from vegan
	  outmat.sd[i,j] <- temp['symmetric.difference']
	  outmat.bsd[i,j] <- temp['branch.score.difference']
	  outmat.pd[i,j] <- temp['path.difference']
	  outmat.qpd[i,j] <- temp['quadratic.path.difference']
	  outmat.RF[i,j] <- RF.dist(treelist[[i]], treelist[[j]])
	  }
	}
  out <- list(symmetric.difference = outmat.sd, branch.score.difference = outmat.bsd, path.difference = outmat.pd, quadratic.path.difference = outmat.qpd, RF = outmat.RF)
  return(out)
  }

make.tree.scores <- function(dna.files, tree.files = NA, dna.format = 'sequential', tree = "none",...) {
## ARGUMENTS
##   dna.files: a vector of files to analyze either on the single tree provided (not currently implemented) or on:
##   tree.files: a vector of tree files corresponding to the dna.files, same order
##   dna.format: format of the dna files, using read.dna
##   tree: a phylogeny of class phylo; not currently implemented
  require(phangorn)
  missingData <- c('n', '-', 'N')
  N = length(dna.files)
  if(is.na(tree.files[1])) {
    warning("No tree files entered; only matrix stats returned.")
	trees = F
	}
  cols.prettyNames <- c('Steps', 'Variable', 'Parsimony informative', 'CI', 'CI on informative characters only', 'Aligned length', 'Percent N or -', 'Percent of tree in terminals', 'Mean bootstraps', 'Standard deviation of bootstraps')
  cols <- c('steps', 'var', 'inf', 'ci', 'ci.inf', 'length', 'missing', 'tree.terms', 'boots.mean', 'boots.sd')
  outMat <- matrix(NA, nrow = N, ncol = length(cols), dimnames = list(names(dna.files), cols))
  for(i in 1:N) {
    message(paste('Doing file', i))
	dna.dat <- read.dna(dna.files[[i]], format = dna.format)
	dna.phyDat <- as.phyDat(dna.dat)
	dna.dat.char <- as.character(dna.phyDat)
	if(trees) {
	  tree <- read.tree(tree.files[[i]])
	  outMat[i, 'steps'] <- parsimony(tree, dna.phyDat)
	  outMat[i, 'ci'] <- sum(phangorn:::lowerBound(dna.phyDat) * attr(dna.phyDat, 'weight')) / outMat[i, 'steps']
	  outMat[i, 'length'] <- sum(attr(dna.phyDat, 'weight'))
	  outMat[i, 'boots.mean'] <- mean(as.numeric(tree$node.label), na.rm = T)
	  outMat[i, 'boots.sd'] <- sd(as.numeric(tree$node.label), na.rm = T)
	  outMat[i, 'tree.terms'] <- sum(tree$edge.length[mrca.branches(tree)]) / sum(tree$edge.length) 
	  }
	uniques <- apply(dna.dat.char, 2, 
					 function(x) {
					   uniqueSet <- unique(x)[!unique(x) %in% missingData]
					   if(length(uniqueSet) < 2) return(c(F, F))
					   else {
					     uniqueLengths <- sapply(uniqueSet, function(y) sum(x == y))
						 if(sum(uniqueLengths > 1) > 1) return(c(T, T))
                         else return(c(T, F))
                         }
                       }
                     )					   
	outMat[i, c('var', 'inf')] <- apply(uniques, 1, sum)  
    outMat[i, 'missing'] <- sum(dna.dat.char %in% missingData) / cumprod(dim(dna.dat.char))[2]
	}
  dimnames(outMat)[[2]] <- cols.prettyNames
  return(outMat)
  }

mrca.branches <- function(tree, repPattern = "_re", origPattern = "_h") {
# returns the branches subtending the mrca of the pairs in the mrcaList
  reps <- grep(repPattern, tree$tip.label, fixed = TRUE, value = TRUE)
  orig <- gsub(repPattern, origPattern, reps, fixed = TRUE)
  mrca.mat <- cbind(reps, orig)
  mrca.tr <- mrca(tree)
  out <- apply(mrca.mat, 1, function(x) which(tree$edge[, 2] == mrca.tr[x[1], x[2]]))
  return(out)
  }
  
					   
do.EST.phylo <- function(dat = d6m10p2.wRE.mat, lengths = d6m10p2.wRE$radSummary$locus.lengths, blast = blast.results.concat, t.hold = -15, minLength = 50, maxLength = 55, sampleReps = 0, makeDirs = TRUE) {
## this needs to spit out paths and loci for each analysis
  if(is.na(lengths[1])) lengths <- lengths.report(dat,0) # as written, this will fail if lengths = NA, b/c lengths.report needs a pyRAD object, whereas this function needs a pyRAD.mat object
  sizeFilteredLoci <- names(lengths[lengths >= minLength & lengths <= maxLength])
  estLoci <- loci.by.threshold(blast, t.hold)
  estLoci <- estLoci[which(estLoci %in% sizeFilteredLoci)]
  otherLoci <- sizeFilteredLoci[!sizeFilteredLoci %in% unique(blast[[1]])]
  if(makeDirs) dir.create(paste('./oak.ests.tHold.', t.hold, '.', format(Sys.time(), "%Y-%d-%m"), '/', sep = ''))
  analysisFile <- file(paste('oak.ests.tHold.', t.hold, '.', format(Sys.time(), "%Y-%d-%m"), '.bat', sep = ''), open = "wt")
  outfileName <- paste('oak.ests.tHold.', t.hold, '.m', maxLength, '.', minLength, '.', format(Sys.time(), "%Y-%d-%m"), '.phy', sep = '')
  rad2phy(dat, loci = estLoci, outfile = outfileName)
  writeLines(paste('raxmlHPC-PTHREADS -f a -T 4 -x 123555 -# 200 -s ', outfileName, ' -m GTRGAMMA -n ', outfileName, '.tre', sep = ''), con = analysisFile)
  if (sampleReps == 0) return(0) #aborts if no sample reps are requested
  for(i in 1:sampleReps) {
    message(paste("*** writing rep", i))
	sampledLoci <- sample(otherLoci, length(estLoci))
	outfileName <- paste('oak.otherLoci.tHold.', t.hold, '.m', maxLength, '.', minLength, '.rep', i, '.', format(Sys.time(), "%Y-%d-%m"), '.phy', sep = '')
	rad2phy(dat, loci = sampledLoci, outfile = outfileName)
    if(makeDirs) {
	  dirOut <- paste('./oak.others.tHold.', t.hold, '.rep', i, '.', format(Sys.time(), "%Y-%d-%m"), '/', sep = '')
	  dir.create(dirOut)
	  }
	else dirOut <- ''
	writeLines(paste('raxmlHPC-PTHREADS -f a -T 4 -x 123555 -# 200 -w ', dirOut, ' -s ', outfileName, ' -m GTRGAMMA -n ', outfileName, '.tre', sep = ''), con = analysisFile)
	}
  return(0)
  }

## Make a helpful barplot:
 barplot(sapply(-1:-30, function(x) length(loci.by.threshold(threshold = x))), names.arg=-1:-30, xlab = "log10(e-value)", ylab = "Number of loci for which best blast is this good or better")
 
  
## OLD -- Redone on 2013-02-21

## read in and concatenate blast results to EST databases
#blast.results <- lapply(dir('../../../OAK-RAD-ANNOTATIONS/QuercusFeb0513/', full = T, patt = 'blastN'), read.delim, header = F, as.is = T)
#names(blast.results)<- (dir('../../../OAK-RAD-ANNOTATIONS/QuercusFeb0513/', full = F, patt = 'blastN')
#blast.results.concat <- blast.results[[1]]
#for(i in 2:length(blast.results)) blast.results.concat <- rbind(blast.results.concat, blast.results[[i]])


## Making a data matrix for the EST-linked markers
# oaks.d6m4.mat <- rad2mat(oaks.d6m4)
# rad2phy(oaks.d6m4.mat, loci = loci.by.threshold(), outfile = 'oaks.dna.d6m4.blast.e-15.phy')
# rad2phy(oaks.d6m4.mat, loci = loci.by.threshold(threshold = -25), outfile = 'oaks.dna.d6m4.blast.e-25.phy')

## Subsampling loci to see how phylogeny based on non-EST linked markers compares
## PROBLEMATIC: mean length of the EST-linked markers is longer, so this results in a shorter data matrix
# subsampled.loci <- lapply(rep(8833,10), sample, x = dimnames(oaks.d6m4.mat)[[2]][!dimnames(oaks.d6m4.mat)[[2]] %in% blast.results.concat[[1]]])
# for(i in 1:10) rad2phy(oaks.d6m4.mat, loci = subsampled.loci[[i]], outfile = paste('oaks.dna.d6m4.blast.subsample.', i,'.phy', sep = ''))

