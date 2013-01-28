## IOS analyses, Andrew Hipp
## 25 January 2013

oaks.d6m4 <- read.pyRAD('./pyRAD_2013-01-23/c_hipp_d6m4p2_wRE.readloci.txt')
overlap.report.d6m4 <- overlap.report(oaks.d6m4$radSummary$inds.mat)
layout(matrix(1:9,3,3))
apply(overlap.report.d6m4[, c(3,5)], 1, function(x) pie(x, radius = sum(x)/32000, col = c('black', 'white'), labels =""))

library(vegan)
oaks.d6m4.jaccard.MDS <- metaMDS(oaks.d6m4$radSummary$inds.mat[-which(row.names(oaks.d6m4$radSummary$inds.mat) %in% c('>AC_h', '>MI_h')), ], 'jaccard')

cols <- rep('red', dim(oaks.d6m4.jaccard.MDS$points)[1])
cols[grep("_re", row.names(oaks.d6m4.jaccard.MDS$points), fixed = TRUE)] <- 'black'
plot(oaks.d6m4.jaccard.MDS$points, pch = 19, col= cols, cex = 2)
text(oaks.d6m4.jaccard.MDS$points, row.names(oaks.d6m4.jaccard.MDS$points), cex = 0.7, adj = -0.2)