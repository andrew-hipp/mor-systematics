gen.RAD.loci.datasets <-
function(rads, trees = "none", loci = "all", taxa = "all", minTaxa = 4, 
                      onlyVariable = TRUE, fileBase = "DEFAULT", 
					  splitInto = 1, 
					  cores = 2,
					  raxSinglePath = "~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-AVX", 
					  raxMultiPath = "~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-PTHREADS-AVX", 
					  header = "#!/bin/sh")
					  {
  if(taxa == 'all') taxa <- row.names(rads$radSummary$inds.mat)
  # create directory structure
  if(fileBase == 'DEFAULT') fileBase <- format(Sys.time(), "rads.%Y-%m-%d")
  if(!paste(fileBase, ".", (0), sep = '') %in% dir()) lapply(paste(fileBase, ".", (0:splitInto), sep = ''), dir.create) # defaults to making a directory to hold all the files
  
  # initiate log files
  analysisFileOut <- lapply(paste(fileBase, '.0/raxml.batch.', 1:splitInto, '.', fileBase, '.sh', sep = ''), file, open = "a")
  for(i in 1:splitInto) cat(header, '\n', file = analysisFileOut[[i]])
  indexFileOut <- file(paste(fileBase, '.0/tree.index.lines.txt', sep = ''), 'a')
  
  # make full analysis output
  fullMatrixAnalysisFile <- file(paste(fileBase, '.0/fullMatrix.sh', sep = ''), open = "a")
  cat(header, '\n', file = fullMatrixAnalysisFile)
  datFileOut <- paste(fileBase, '.0/fullMatrix.phy', sep = '')
  rad2phy(rad2mat(rads), taxa, filter.by(rads, taxa = taxa, threshold = minTaxa), datFileOut)
  treeFileOut <- paste(fileBase, '.0/fullTreeSet.tre', sep = '')
  write.tree(trees, treeFileOut)
  analysisLine <- paste(raxMultiPath, "-f G -s", paste('../', datFileOut, sep = ''), "-T", cores, "-m GTRGAMMA -z", paste('../', treeFileOut, sep = ''), "-n", paste(i, '.lnL', sep = ''))
  cat(analysisLine, '\n', file = fullMatrixAnalysisFile)
  close(fullMatrixAnalysisFile)
  
  # subset loci and trees
  if(loci[1] == "all") loci <- unique(rads$locus.index)[unique(rads$locus.index) != ""]
  if(trees[1] != 'none') taxa <- intersect(taxa, trees[[1]]$tip.label)
  if(length(taxa) <- 0) error('no taxa match between your RAD and tree datasets: please check names and try again')
  locus.set <- subset.pyRAD.loci(rads, loci, taxa)
  locus.list <- locus.set$DNA[names(which(locus.set$ntaxa >= minTaxa))]
  if(onlyVariable) locus.list <- locus.list[names(which(locus.set$variable))]
  if(trees[1] != 'none') tree.vector.matrix <- matrix(NA, nrow = length(locus.list), ncol = length(trees), dimnames = list(names(locus.list), names(trees)))
  
  # subset each locus, write them out
  batch = 0
  trees <- lapply(trees, function(x) x) #this is a workaround -- lapply wasn't workind correctly on multiPhylo object from nni
  for(i in names(locus.list)) {
    if(batch == splitInto) batch <- 1
	  else batch <- batch + 1
	message(paste('Doing', i))
	locus.taxa <- names(locus.list[[i]])[names(locus.list[[i]]) %in% taxa]
	datFileOut <- paste(fileBase, '.', batch, '/', i, '.phy', sep = '')
	write.DNAStringSet(locus.list[[i]][locus.taxa], filename = datFileOut)
	if(trees[1] != 'none') {
	  trees.out <- try(lapply(trees, drop.tip, tip = trees[[1]]$tip.label[!trees[[1]]$tip.label %in% locus.taxa]))
	  if(class(trees.out) == "try-error") {
	    message('...error with drop.tip -- bailing out...')
		next
		}
	  trees.out <- try(lapply(trees.out, unroot))
	  if(class(trees.out) == "try-error") {
	    message('...error with unroot -- bailing out...')
		next
		}
	  class(trees.out) <- 'multiPhylo'
	  trees.out <- unique(trees.out) # this really slows things down... if there were a way to speed this up it w/b great.
	  message(paste('... kept', length(trees.out), 'trees'))
	  treeFileOut <- paste(fileBase, '.', batch, '/', i, '.tre', sep = '')
	  write.tree(trees.out, file = treeFileOut)
	  cat(i, '\t', paste(attr(trees.out, "old.index"), collapse = '\t'), '\n', sep = '', file = indexFileOut)
	  }
	analysisLine <- paste(raxSinglePath, "-f G -s", paste('../', datFileOut, sep = ''), "-m GTRGAMMA -z", paste('../', treeFileOut, sep = ''), "-n", paste(i, '.lnL', sep = ''))
	cat(analysisLine, '\n', file = analysisFileOut[[batch]])
    }
  # Close all log and batch files
  for (i in analysisFileOut) close(i)
  close(indexFileOut)
  }
