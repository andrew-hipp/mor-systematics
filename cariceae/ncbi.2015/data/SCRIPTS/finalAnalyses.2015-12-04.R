## all code post-Kasey Pham clustering
## 2015-12-04
## ahipp@mortonarb.org

require(ape)
require(phytools)
require(treescape)
require(parallel)
require(vegan)

readData = T
cleanData = T
analyzeAndPlot = T

if(readData) {
  message('\nREADING DATA')
  ncbi.meta <- read.delim('../METADATA/ncbi.meta.2015.09.03.tsv', row.names = 1, as.is = T)
  source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/phyloBuddy/R/make.matrices.R')
  ncbi.geneMat <- make.gene.matrix(ncbi.meta)
  tr.bipart <- lapply(dir('../TREES/', full = T), function(x) read.tree(dir(x, patt = 'bipartitions\\.', full = T)))
  tr.boots <- lapply(dir('../TREES/', full = T), function(x) read.tree(dir(x, patt = 'boot', full = T)))
  names(tr.bipart) <- names(tr.boots) <- gsub('trees.', '', dir('../TREES/'), fixed = T)
  for(i in names(tr.bipart)) tr.bipart[[i]]$tip.label <- as.character(label.elements(tr.bipart[[i]], returnNum = 1:2, returnDelim = "|", fixed = T ))
  dna.5 <- read.dna('../MATRICES/5_gene.fas', format = 'fasta')
  dna.12 <- read.dna('../MATRICES/12_gene.fas', format = 'fasta')
  dna.itsScaff <- read.dna('../MATRICES/ITS_scaffold.fas', format = 'fasta')
  alignments <- lapply(dir('../ALIGNMENTS', full = TRUE), read.dna, format = 'fasta')
  names(alignments) <- sapply(strsplit(dir('../ALIGNMENTS'), '.2015.', fixed = T), function(x) x[[1]][1])
  } # close readData

if(cleanData) {
  message('\nCLEANING TREES')
  deletesies <- c('Kobresia_royleana|AUTHOR_Zheng_1514', 'Kobresia_macrantha|AUTHOR_Zheng_1509',
                  'Kobresia_royleana|AUTHOR_Zheng_1544', 'Kobresia_royleana|AUTHOR_Zheng_1543',
                  'Kobresia_royleana|AUTHOR_Zheng_1515', 'Schoenoxiphium_caricoides|Gehrke_3250|Williams_7265',
                 'Carex_atrofusca|AUTHOR_Taberlet_1397', "Carex_siderosticta|Leveille_Bourret_4776")
  og <- 'Carex_siderosticta|Waterway_7159'
  #ogPatt <- ('sidero|tsiangii|scaposa|lingii|scaposa|densifimbriata|glossostigma|longshengensis|kwangsiensis|okamotoi|ciliatomarginata|pachygyna|subcapitata|tumidula')
  for(i in names(tr.bipart)) {
    message(paste('doing bipartition', i))
    tr.bipart[[i]] <- drop.tip(tr.bipart[[i]], deletesies)
    tr.bipart[[i]] <- drop.tip(ladderize(root(tr.bipart[[i]], og)), og)
    message(paste('.... doing bootstraps', i))
    tr.boots[[i]] <- lapply(tr.boots[[i]], function(x) drop.tip(ladderize(root(drop.tip(x, deletesies), og)), og))
    class(tr.boots[[i]]) <- 'multiPhylo'
  }
} # close cleanData

if(analyzeAndPlot) {

message('\nANALYZING AND MAKING FIGURES')

## SUPPLEMENTAL FIGURES: ALL MATRICES, INDIVIDUAL TREES
lapply(names(tr.bipart), function(i) {
  pdf(paste("SUPPLEMENT.FIG-S3", i, 'pdf', sep = '.'), 10, 45)
  plot(tr.bipart[[i]], cex = 0.3, show.node.label = T)
  dev.off()
  }
  )

## SUPPLEMENTAL FIGURE: DNA MATRIX IMAGES
jpeg('SUPPLEMENT.FIG-S2.DNA.matrices.2015-12-09.jpg', 2550, 3300, quality = 90)
layout(matrix(1:3, 3, 1))
image(dna.12, show.labels = F, legend = F, main = '12 genes matrix')
image(dna.5, show.labels = F, legend = F, main = '5 genes matrix')
image(dna.itsScaff, show.labels = F, legend = F, main = 'ITS scaffold')
dev.off()

## SUPPLEMENTAL TABLE 1: INDIVIDUAL ALIGNMENT STATS
tempAlignments <- alignments
tempAlignments$dna.12 = dna.12
tempAlignments$dna.5 = dna.5
tempAlignments$dna.itsScaff = dna.itsScaff

table.suppl.alignments <- cbind(
  t(sapply(tempAlignments, dim)),
  missingDatProportion = sapply(tempAlignments, function(x) round(sum(as.character(x) %in% c('n', 'N', '-')) / (dim(x)[1] * dim(x)[2]), 3))
  )
rm(tempAlignments)
dimnames(table.suppl.alignments)[[2]][1:2] <- c('numberOfTips','width.bp')
write.csv(table.suppl.alignments, 'SUPPLEMENT.TABLE-S1.alignmentsStats.csv')

## METHODS: How many vouchers with > 1 organism?
ncbi.meta.split <- split(ncbi.meta, ncbi.meta$cleanedVoucher)
ncbi.organisms <- sapply(ncbi.meta.split, function(x) length(unique(x$organism)))
sum(ncbi.organisms == 1) # 3967
sum(ncbi.organisms == 2) # 18
sum(ncbi.organisms > 2)  # 0

## METHODS: how any unique gene regions?
unique(ncbi.meta$cleanedGeneRegion) # 60 gene regions

## FIGURE XX
## gene matrix plot
source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/phyloBuddy/R/plot.geneMat.R')
ncbi.geneVec <- colnames(ncbi.geneMat[5:64])
plot.geneMat(ncbi.geneMat, tr.bipart$'12_gene', ncbi.geneVec,
                                        pdfTitle = paste('FIG4.12-region.geneMat.plot', i, format(Sys.time(), "%Y.%m.%d.%H%M.pdf"), sep = '.'),
                                        pdfW = 7, pdfH = 7, pch = 19, panes = c(3, 1, 2))

## FIGURE 1
## sequences by date plot
pdf('FIG1.sequence.plots.byDate.2015-11-30.pdf',7,5)
par(mar = c(5,5,1,1) + 0.1)
plot(y=1:7494, x=sort(as.Date(ncbi.meta$create_date, format = "%d-%b-%y"), decreasing=F), pch = 19, cex = 0.5, ylab = 'Cumulative number of sequences deposited\nin NCBI Genbank through 1 March 2015', xlab = 'Year of sequence record creation in NCBI Genbank')
dev.off()

## FIGURE 2
## gene sampling
pdf('FIG2.geneSampling.2015-11-30.pdf', 3.375, 7)
par(mar = c(5,6,0,1) + 0.1)
barplot(sort(table(ncbi.meta$cleanedGeneRegion)[sort(unique(ncbi.meta$cleanedGeneRegion), decreasing = T)], decreasing=T), horiz=T, xlim = c(0, 2000), las = 1, xlab = "Number of sequences submitted to\nNCBI GenBank per gene region", cex.lab = 0.6, cex.axis = 0.5, cex = 0.5)
dev.off()

## FIGURE 3
## gene regions per individual
ncbi.genes <- sapply(ncbi.meta.split, function(x) length(unique(x$cleanedGeneRegion)))
hist(ncbi.genes)
table(ncbi.genes)
pdf('FIG3.genesPerIndlBarplot.2015-12-09.pdf', 3.75, 3.75)
barplot(table(ncbi.genes), ylab = 'Carex specimens in NCBI', xlab = 'Gene regions sampled per individual')
abline(v = mean(ncbi.genes), lty = 'dashed')
dev.off()

## FIGURE XX
## ITS vs ETS vs combined bootstrap comparison
source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/phyloBuddy/R/match.nodes.R')
its.ets.tips <- intersect(intersect(tr.bipart$ITS_only$tip.label, tr.bipart$ETS_only$tip.label), tr.bipart$ITS_ETS$tip.label)
its.ets.matched <- match.nodes(drop.tip(tr.bipart$ITS_ETS, which(!tr.bipart$ITS_ETS$tip.label %in% its.ets.tips)),
                               list(its = drop.tip(tr.bipart$ITS_only, which(!tr.bipart$ITS_only$tip.label %in% its.ets.tips)),
                                    ets = drop.tip(tr.bipart$ETS_only, which(!tr.bipart$ETS_only$tip.label %in% its.ets.tips))),
                               plotBoots = FALSE)
pdf('FIGXX.its.ets.bootComparison.2015-12-04.pdf', 3.375, 4)
matplot(its.ets.matched$mat.boots[order(its.ets.matched$mat.boots[, 1], decreasing = T), ], type = 'l',
        lwd = c(2,0.5,0.5), lty = c('solid'), col = c('black', 'red', 'blue'),
        xlab = 'node (arbitrary)', ylab='Bootstrap', cex.axis = 0.5)
legend(x=0, y=20, legend=c('ITS + ETS (concatenated matrix)', 'ITS only', 'ETS only'),
      bty = 'n', lwd = c(2, 0.5, 0.5), lty = c('solid'), col = c('black', 'red', 'blue'), cex = 0.5)
dev.off()

## FIGURES Sxx and xx
## taxon disparity
source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/phyloBuddy/R/trees.R')
disparity.mats <- lapply(tr.bipart, function(x) summary.by.elements(x, delim = '[|_]', returnDelim = '_', returnNum = 1:2))
pdf('SUPPLEMENT.XX.disparityMatricesHistograms.2015-12-04.pdf', 8.5, 11)
layout(matrix(1:6, 3, 2))
for(i in names(disparity.mats)) {
    a = hist(disparity.mats[[i]]$disparity.mat[disparity.mats[[i]]$disparity.mat[, 'count'] > 1, 'disparity'], 20,
             main = paste('Dataset:', i), xlab = 'Taxonomic disparity', ylab = 'Number of unique species with > 1 tip')
    text(max(a$breaks), 0.90 * max(a$counts),
        paste("N = ", sum(disparity.mats[[i]]$disparity.mat[, 'count']), " tips",
              "\nN = ", dim(disparity.mats[[i]]$disparity.mat)[1], " species",
              '\n', round(mean(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), 1), " +/- ", round(sd(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), 2), ' (s.d.) tips per species',
              "\n", sum(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), " spp represented by > 1 tip",
              "\n", round(100 * sum((disparity.mats[[i]]$disparity.mat[, 'count'] > 1) & (disparity.mats[[i]]$disparity.mat[, 'disparity'] == 0)) / sum(disparity.mats[[i]]$disparity.mat[, 'count'] > 1)), "% monophyletic",
              sep = ''), pos = 2)
}
rm(a)
dev.off()

pdf('FIGXX.disparityMatricesHistogramsForMS.2015-12-10.pdf', 3.375, 9)
layout(matrix(1:3, 3, 1))
for(i in c('12_gene', '5_gene', 'ITS_scaffold')) {
    a = hist(disparity.mats[[i]]$disparity.mat[disparity.mats[[i]]$disparity.mat[, 'count'] > 1, 'disparity'], 20,
             main = paste('Dataset:', gsub('gene', 'region', gsub("_", " ", i, fixed = T))), xlab = 'Taxonomic disparity', ylab = 'Number of unique species with > 1 tip')
    text(max(a$breaks), 0.85 * max(a$counts),
        paste("N = ", sum(disparity.mats[[i]]$disparity.mat[, 'count']), " tips",
              "\nN = ", dim(disparity.mats[[i]]$disparity.mat)[1], " species",
              '\n', round(mean(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), 1), " +/- ", round(sd(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), 2), ' (s.d.) tips per species',
              "\n", sum(disparity.mats[[i]]$disparity.mat[, 'count'] > 1), " spp represented by > 1 tip",
              "\n", round(100 * sum((disparity.mats[[i]]$disparity.mat[, 'count'] > 1) & (disparity.mats[[i]]$disparity.mat[, 'disparity'] == 0)) / sum(disparity.mats[[i]]$disparity.mat[, 'count'] > 1)), "% of species monophyletic",
              sep = ''), pos = 2, cex = 0.6)
}
rm(a)
dev.off()

## FIGURE XX
## ordination -- this takes a bit of time, so be careful before starting
### and now using topo, from ape
a = vector('list')
for(i in 1:302) {
  for(j in (i+1): 303) {
    a[[length(a) + 1]] <- c(i, j)}}
ordination.topo <- mclapply(a, function(x) dist.topo(ordination.trees[[x[1]]], ordination.trees[[x[2]]]), mc.cores = 28)
ordination.treesUsed <- matrix(NA, 45753, 2, dimnames = list(1:45753, c('tree1', 'tree2')))
counter = 1
for(i in 1:302) {
  for(j in (i+1): 303) {
    ordination.treesUsed[counter, ] <- c(i, j)
    counter <- counter + 1
    }}
save(ordination.topo, ordination.treesUsed, file = 'ordination.topo.v2.Rdata')
ordination.topo.mat <- matrix(NA, 303,303)
for(i in 1:length(ordination.topo)) ordination.topo.mat[ordination.treesUsed[i,2], ordination.treesUsed[i,1]] <- ordination.topo[[i]]
ordination.topo.mds <- monoMDS(as.dist(ordination.topo.mat))

pdf('FIGXX.dist.topo.plot.2015-12-10.pdf', 7, 5)
plot(ordination.topo.mds$points, pch = c(rep(21,3), rep(24, 300)), col = c('red', 'blue', 'gray', rep('red', 100), rep('blue', 100), rep('gray', 100)))
points(ordination.topo.mds$points[1:3,], pch = 21, bg = c('red', 'blue', 'gray'), cex = 1)
legend(-2, -1, c('5-gene ML', '12-gene ML', 'ITS scaffold ML', '5-gene bootstraps', '12-gene bootstraps', 'ITS scaffold bootstraps'),
       pch = c(21,21,21,24,24,24), pt.cex = 1,
       pt.bg = c('red', 'blue', 'gray', NA,NA,NA), col = c(rep('black', 3), 'red', 'blue', 'gray'), bty = 'n', cex = 0.6)
dev.off()

### and now using topo, from ape
a = vector('list')
for(i in 1:302) {
  for(j in (i+1): 303) {
    a[[length(a) + 1]] <- c(i, j)}}
ordination.topo <- mclapply(a, function(x) dist.topo(ordination.trees[[x[1]]], ordination.trees[[x[2]]]), mc.cores = 28)
ordination.treesUsed <- matrix(NA, 45753, 2, dimnames = list(1:45753, c('tree1', 'tree2')))
counter = 1
for(i in 1:302) {
  for(j in (i+1): 303) {
    ordination.treesUsed[counter, ] <- c(i, j)
    counter <- counter + 1
    }}
save(ordination.topo, ordination.treesUsed, file = 'ordination.topo.v2.Rdata')
ordination.topo.mat <- matrix(NA, 303,303)
for(i in 1:length(ordination.topo)) ordination.topo.mat[ordination.treesUsed[i,2], ordination.treesUsed[i,1]] <- ordination.topo[[i]]
ordination.topo.mds <- monoMDS(as.dist(ordination.topo.mat))

pdf('FIGXX.dist.topo.plot.2015-12-10.pdf', 7, 5)
plot(ordination.topo.mds$points, pch = c(rep(21,3), rep(24, 300)), col = c('red', 'blue', 'gray', rep('red', 100), rep('blue', 100), rep('gray', 100)))
points(ordination.topo.mds$points[1:3,], pch = 21, bg = c('red', 'blue', 'gray'), cex = 1)
legend(-2, -1, c('5-gene ML', '12-gene ML', 'ITS scaffold ML', '5-gene bootstraps', '12-gene bootstraps', 'ITS scaffold bootstraps'),
       pch = c(21,21,21,24,24,24), pt.cex = 1,
       pt.bg = c('red', 'blue', 'gray', NA,NA,NA), col = c(rep('black', 3), 'red', 'blue', 'gray'), bty = 'n', cex = 0.6)
dev.off()


## first using treescape...
ordination.tips <- intersect(intersect(tr.bipart$'5_gene'$tip.label, tr.bipart$'12_gene'$tip.label), tr.bipart$ITS_scaffold$tip.label) ## should match tips in ITS_scaffold, which it does
ordination.trees <- c(c(tr.bipart$'5_gene', tr.bipart$'12_gene', tr.bipart$ITS_scaffold), tr.boots$'5_gene', tr.boots$'12_gene', tr.boots$ITS_scaffold)
for(i in 1:length(ordination.trees)) ordination.trees[[i]] <- drop.tip(ordination.trees[[i]], ordination.trees[[i]]$tip.label[!ordination.trees[[i]]$tip.label %in% ordination.tips])
ordination.treescape <- treescape(ordination.trees, nf = 3)
pdf('FIGXX.treescape.plot.2015-12-04.pdf', 7, 5)
plot(ordination.treescape$pco$li[, 1:2], pch = c(rep(21,3), rep(24, 300)), col = c('red', 'blue', 'gray', rep('red', 100), rep('blue', 100), rep('gray', 100)))
points(ordination.treescape$pco$li[1:3, 1:2], pch = 21, bg = c('red', 'blue', 'gray'), cex = 1)
legend(-8500, -3500, c('5-gene ML', '12-gene ML', 'ITS scaffold ML', '5-gene bootstraps', '12-gene bootstraps', 'ITS scaffold bootstraps'),
       pch = c(21,21,21,24,24,24), pt.cex = 1,
       pt.bg = c('red', 'blue', 'gray', NA,NA,NA), col = c(rep('black', 3), 'red', 'blue', 'gray'), bty = 'n', cex = 0.6)
dev.off()

### and now using topo, from ape
a = vector('list')
for(i in 1:302) {
  for(j in (i+1): 303) {
    a[[length(a) + 1]] <- c(i, j)}}
ordination.topo <- mclapply(a, function(x) dist.topo(ordination.trees[[x[1]]], ordination.trees[[x[2]]]), mc.cores = 28)
ordination.treesUsed <- matrix(NA, 45753, 2, dimnames = list(1:45753, c('tree1', 'tree2')))
counter = 1
for(i in 1:302) {
  for(j in (i+1): 303) {
    ordination.treesUsed[counter, ] <- c(i, j)
    counter <- counter + 1
    }}
save(ordination.topo, ordination.treesUsed, file = 'ordination.topo.v2.Rdata')
ordination.topo.mat <- matrix(NA, 303,303)
for(i in 1:length(ordination.topo)) ordination.topo.mat[ordination.treesUsed[i,2], ordination.treesUsed[i,1]] <- ordination.topo[[i]]
ordination.topo.mds <- monoMDS(as.dist(ordination.topo.mat))

pdf('FIGXX.dist.topo.plot.2015-12-10.pdf', 7, 5)
plot(ordination.topo.mds$points, pch = c(rep(21,3), rep(24, 300)), col = c('red', 'blue', 'gray', rep('red', 100), rep('blue', 100), rep('gray', 100)))
points(ordination.topo.mds$points[1:3,], pch = 21, bg = c('red', 'blue', 'gray'), cex = 1)
legend(-2, -1, c('5-gene ML', '12-gene ML', 'ITS scaffold ML', '5-gene bootstraps', '12-gene bootstraps', 'ITS scaffold bootstraps'),
       pch = c(21,21,21,24,24,24), pt.cex = 1,
       pt.bg = c('red', 'blue', 'gray', NA,NA,NA), col = c(rep('black', 3), 'red', 'blue', 'gray'), bty = 'n', cex = 0.6)
dev.off()

require(phytools)
tr.cophylo5v12 <- cophylo(ordination.trees[[1]], ordination.trees[[2]], assoc = cbind(ordination.tips, ordination.tips), rotate = T)
save(tr.cophylo5v12, file = 'tr.cophylo5v12.Rdata')


pdf('FIG4.cophyloplot.5v12.pdf', 7, 9)
plot(tr.cophylo5v12, ftype = 'off')
text(c(-0.25, 0.25),rep(1.02, 2), c("5-region tree", '12-region tree'), cex = 0.8)
rect(-0.34,-0.01,0.34, 0.5304, border = 'blue', lwd = 1.5)
text(-0.36, 0.26, "Core Carex", col = 'blue', cex = 0.8, srt = 90)
text(0.36, 0.26, "Core Carex", col = 'blue', cex = 0.8, srt = 270)
rect(-0.38,0.535,-0.01, 0.6422, border = 'red', lwd = 1.5)
text(-0.40, 0.592, "Caricoid", col = 'red', cex = 0.8, srt = 90)
rect(0.01,0.535,0.38, 0.8717, border = 'grey', lwd = 1.5)
text(0.40, mean(c(0.8717, 0.5405)), 'Vignea', col = 'grey', cex = 0.8, srt = 270)
rect(-0.37, 0.6468, -0.01, 0.9841, border = 'grey', lwd = 1.5)
text(-0.39, mean(c(0.648, 0.9841)), 'Vignea', col = 'grey', cex = 0.8, srt = 90)
rect(0.01, .8763, 0.37, 0.9841, border = 'red', lwd = 1.5)
text(0.39, mean(c(0.875, 0.9841)), 'Caricoid', col = 'red', cex = 0.8, srt = 270)
#rect(-0.51,0.985,0.51, 1.01, border = 'green', lwd = 1.5)
dev.off()

tr.cophylo5v12.taxa <- cbind(tree.5 = read.tree(text = write.tree(tr.cophylo5v12$trees[[1]]))$tip.label, 
                             tree.12 = read.tree(text = write.tree(tr.cophylo5v12$trees[[2]]))$tip.label,
							 hline = (1:length(ordination.tips))/ length(ordination.tips))
write.csv(tr.cophylo5v12.taxa, 'tr.cophylo5v12.taxa.csv')


pdf('SUPPLEMENT.FIGx.cophyloplot.5vITSscaff.pdf', 7, 9)
require(phytools)
tr.cophylo5vITSscaff <- cophylo(ordination.trees[[1]], ordination.trees[[3]], assoc = cbind(ordination.tips, ordination.tips), rotate = T)
save(tr.cophylo5vITSscaff, file = 'tr.cophylo5vITSscaff.Rdata')
plot(tr.cophylo5vITSscaff, ftype = 'off')
dev.off()
