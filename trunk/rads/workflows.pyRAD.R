## Workflows and wrappers for other functions
## from the r2pyRAD library
## Andrew Hipp, The Morton Arboretum, 2014

gen.RAD.loci.datasets <- function(rads, trees = 'none', loci = 'all', taxa = 'all', minTaxa = 3, onlyVariable = TRUE, fileBase = format(Sys.time(), "rads.%Y-%m-%d"), splitInto = 8) {
  if(!fileBase %in% dir()) lapply(paste(fileBase, ".", (0:splitInto), sep = ''), dir.create) # defaults to making a directory to hold all the files
  if(loci[1] == "all") loci <- unique(rads$locus.index)[unique(rads$locus.index) != ""]
  if(taxa[1] == "all") taxa <- unique(rads$tips)[gsub('/', '', unique(rads$tips), fixed = TRUE) != ''] # gets rid of the mock tip that is left by the current version of read.pyRAD
  if(trees[1] != 'none') taxa <- intersect(taxa, trees[[1]]$tip.label)
  locus.set <- subset.pyRAD.loci(rads, loci, taxa)
  locus.list <- locus.set$DNA[names(which(locus.set$ntaxa >= minTaxa))]
  if(onlyVariable) locus.list <- locus.list[names(which(locus.set$variable))]
  if(trees[1] != 'none') tree.vector.matrix <- matrix(NA, nrow = length(locus.list), ncol = length(trees), dimnames = list(names(locus.list), names(trees)))
  batch = 0
  for(i in names(locus.list)) {
    if(batch == splitInto) batch <- 1
	  else batch <- batch + 1
	message(paste('Writing', i))
	locus.taxa <- names(locus.list[[i]])[names(locus.list[[i]]) %in% taxa]
	write.DNAStringSet(locus.list[[i]][locus.taxa], filename = paste(fileBase, '.', batch, '/', i, '.phy', sep = ''))
	if(trees[1] != 'none') {
	  trees.out <- lapply(trees, drop.tip, tip = trees[[1]]$tip.label[!trees[[1]]$tip.label %in% locus.taxa])
	  class(trees.out) <- 'multiPhylo'
	  trees.out <- unique(trees.out)
	  message(paste('... kept', length(trees.out), 'trees'))
	  write.tree(trees.out, file = paste(fileBase, '.', batch, '/', i, '.tre', sep = ''))
	  tree.vector.matrix[i, ] <- attr(trees.out, "old2new.trees")
	  }
  }
  write.csv(tree.vector.matrix, paste(fileBase, '.0/', 'tree.vector.matrix.csv', sep = ''))
  ## AND CREATE BATCH FILE...
  ## AND EXPORT!
  ## There will also need to be a read function
  }