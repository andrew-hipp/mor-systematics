genTrees <-
function(x, N = 200, filebase = 'trial', method = c('nni', 'random'), maxmoves = 2, perms = c(length(nni(x)), as.integer(N-(length(nni(x))))), software = c('raxml', 'paup'), ...) {
  ## Arguments:
  ## x = phylo tree
  ## N = total number of trees to generate
  ## filebase = file name base; a tree file (.tre) and paup command file (.nex) will be created for both
  ## method = method for generating trees
  ## maxmoves = maximum number of rearrangements per tree for nni or spr
  ## perms = number of permutations per maxmoves
  ## ... = additional arguments to pass along to rtree.phylo or rNNI
  ## works with nni, 12 nov 10
  ## January 2014: as written, this doesn't unroot the tree. It ought to, unless you are evaluating trees in a rooted framework (e.g., not using GTR)
  ## 20 January 2014: updated to make sure all trees are unique
  if(class(x) != 'phylo') stop('This function requires a phylo object as its first argument')
  if(method[1] == 'nni') {
	for(i in seq(maxmoves)) {
	  message(paste('doing maxmoves', i))
	  if(i == 1) {
	    originalTree <- list(x)
		treeset <- c(originalTree, lapply(nni(x), function(x) x))
		}
	  # else treeset <- c(treeset, rNNI(x, i, perms[i]))
	  else treeset <- c(treeset, lapply(unique(rNNI(x, i, perms[i] * 1.5)), function(x) x))
	  treeset <- lapply(treeset, unroot)
	  class(treeset) <- 'multiPhylo'
	  treeset <- unique(treeset)[1:sum(perms[1:i], 1)]
	  # just takes the first set of uniques... chops off non-uniques presented so far
      }	# end i
	} # end if	  
  else if(method[1] == 'random') treeset = c(x, rtree.phylo(x, N, ...))
  class(treeset) <- 'multiPhylo'
  if(software[1] == 'raxml') {
    browser()
	message('writing raxml')
	write.tree(treeset, file = paste(filebase, '.trees.tre', sep = ''))
    # write.tree(x, file = paste(filebase, '.optimal.tre', sep = '')) ## no longer separating optimal from full trees
    }
  if(software[1] == 'raxml') {
    message('RAxML chosen as analysis software. Currently, you just need to run this on your own to get the site likelihoods.Try something like this:\n
	/home/andrew/code/raxml/standard-RAxML-7.7.2/raxmlHPC-PTHREADS-SSE3 -f g -T 10 -s d6m10.phy -m GTRGAMMA -z analysis.d6m10/RAxML_bestTree.d6m10.out -n d6m10.phy.reduced.siteLnL')
	}
  return(treeset)
  }
