unique.multiPhylo <- function (x, incomparables = FALSE, use.edge.length = FALSE, use.tip.label = TRUE, verbose = FALSE, ...) {
    ## edited 2014-01-13, AH
    n <- length(x)
    keep <- 1L
	tree.vector <- integer(n) # this is a vector telling which new tree each old tree belongs to
    tree.vector[1] <- 1
	for (i in 2:n) {
        if(verbose) message(paste('--- Doing tree', i))
        already.seen <- FALSE
        tree.vector[i] <- length(keep) + 1 # defaults to its own tree unless this is changed
		for (s in 1:length(keep)) { # this line changed...
			if (all.equal(x[[keep[s]]], x[[i]], use.edge.length = use.edge.length, use.tip.label = use.tip.label, ...)) ### and this line
                {
                already.seen <- TRUE
                if(verbose) message(paste('breaking out of tree', i, 'which is not unique'))
                tree.vector[i] <- s
				break
            }
        }
        if (!already.seen)
            keep <- c(keep, i)
    }
    out <- x[keep]
    attr(out, 'old2new.trees') <- tree.vector
	out
}