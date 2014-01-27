gen.RAD.loci.datasets <-
function(rads, trees = 'none', loci = 'all', taxa = 'all', minTaxa = 4, onlyVariable = TRUE, fileBase = format(Sys.time(), "rads.%Y-%m-%d"), splitInto = 8, raxPath = "~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-AVX", header = "#!/bin/sh") {
  # create directory structure
  if(!paste(fileBase, ".", (0), sep = '') %in% dir()) lapply(paste(fileBase, ".", (0:splitInto), sep = ''), dir.create) # defaults to making a directory to hold all the files
  
  # initiate log files
  analysisFileOut <- lapply(paste(fileBase, '.0/raxml.batch.', 1:splitInto, '.', fileBase, '.sh', sep = ''), file, open = "a")
  for(i in 1:splitInto) cat(header, '\n', file = analysisFileOut[[i]])
  indexFileOut <- file(paste(fileBase, '.0/tree.index.lines.txt', sep = ''), 'a')
  
  # subset loci and trees
  if(loci[1] == "all") loci <- unique(rads$locus.index)[unique(rads$locus.index) != ""]
  if(taxa[1] == "all") taxa <- unique(rads$tips)[gsub('/', '', unique(rads$tips), fixed = TRUE) != ''] # gets rid of the mock tip that is left by the current version of read.pyRAD
  if(trees[1] != 'none') taxa <- intersect(taxa, trees[[1]]$tip.label)
  locus.set <- subset.pyRAD.loci(rads, loci, taxa)
  locus.list <- locus.set$DNA[names(which(locus.set$ntaxa >= minTaxa))]
  if(onlyVariable) locus.list <- locus.list[names(which(locus.set$variable))]
  if(trees[1] != 'none') tree.vector.matrix <- matrix(NA, nrow = length(locus.list), ncol = length(trees), dimnames = list(names(locus.list), names(trees)))
  
  # subset each locus, write them out
  batch = 0
  for(i in names(locus.list)) {
    if(batch == splitInto) batch <- 1
	  else batch <- batch + 1
	message(paste('Writing', i))
	locus.taxa <- names(locus.list[[i]])[names(locus.list[[i]]) %in% taxa]
	datFileOut <- paste(fileBase, '.', batch, '/', i, '.phy', sep = '')
	write.DNAStringSet(locus.list[[i]][locus.taxa], filename = datFileOut)
	if(trees[1] != 'none') {
	  trees.out <- lapply(trees, drop.tip, tip = trees[[1]]$tip.label[!trees[[1]]$tip.label %in% locus.taxa])
	  trees.out <- lapply(trees.out, unroot)
	  class(trees.out) <- 'multiPhylo'
	  trees.out <- unique(trees.out) # this really slows things down... if there were a way to speed this up it w/b great.
	  message(paste('... kept', length(trees.out), 'trees'))
	  treeFileOut <- paste(fileBase, '.', batch, '/', i, '.tre', sep = '')
	  write.tree(trees.out, file = treeFileOut)
	  cat(i, '\t', paste(attr(trees.out, "old.index"), collapse = '\t'), '\n', sep = '', file = indexFileOut)
	  }
	analysisLine <- paste(raxPath, "-f G -s", paste('../', datFileOut, sep = ''), "-m GTRGAMMA -z", paste('../', treeFileOut, sep = ''), "-n", paste(i, '.lnL', sep = ''))
	cat(analysisLine, '\n', file = analysisFileOut[[batch]])
    }
  }
