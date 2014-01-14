## Workflows and wrappers for other functions
## from the r2pyRAD library
## Andrew Hipp, The Morton Arboretum, 2014

gen.RAD.loci.datasets <- function(rads, trees = 'none', loci = 'all', taxa = 'all', minTaxa = 3, onlyVariable = TRUE, fileBase = format(Sys.time(), "rads.%Y-%m-%d")) {
  if(!fileBase %in% dir()) dir.create(fileBase) # defaults to making a directory to hold all the files
  if(loci[1] == "all") loci <- unique(rads$locus.index)[unique(rads$locus.index) != ""]
  if(taxa[1] == "all") taxa <- unique(rads$tips)[gsub('/', '', unique(rads$tips), fixed = TRUE) != ''] # gets rid of the mock tip that is left by the current version of read.pyRAD
  if(trees[1] != 'none') taxa <- intersect(taxa, trees[[1]]$tip.label)
  locus.set <- subset.pyRAD.loci(rads, loci, taxa)
  locus.list <- locus.set$DNA[names(which(locus.set$ntaxa >= minTaxa))]
  if(onlyVariable) locus.list <- locus.list[names(which(locus.set$variable))]
  for(i in names(locus.list)) {
    locus.taxa <- names(locus.list[[i]])[names(locus.list[i]) %in% taxa]
	write.DNAStringSet(locus.list[i][locus.taxa], filename = paste(i, '.phy', sep = ''))
	if(trees[1] != 'none') write.tree(lapply(trees, prune, tip = trees[[1]]$tip.label[!trees[[1]]$tip.label %in% locus.taxa]), file = paste(i, '.phy', sep = ''))
  }
  ## AND CREATE BATCH FILE...
  ## AND EXPORT!
  ## There will also need to be a read function
  }