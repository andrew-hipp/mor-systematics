getal <-
function (data) 
{
    x <- dim(data)
    if (max(data[, 2], na.rm = TRUE) < 1e+06) 
        modulo = 1000
    if (max(data[, 2], na.rm = TRUE) < 10000) 
        modulo <- 100
    if (max(data[, 2], na.rm = TRUE) < 100) 
        modulo <- 10
    firstal <- data[, -1]%/%modulo
    secal <- data[, -1]%%modulo
    ind <- vector(length = 0)
    nbpop <- length(unique(data[, 1]))
	# two changes below were necessary... indexed wrong if 2 precedes 1 in genotype tables
    # for (i in unique(data[, 1])) {
	for (i in 1:length(unique(data[, 1]))) {
        # dum <- 1:sum(data[, 1] == i)
		dum <- 1:sum(data[, 1] == unique(data[, 1])[i])
        if (i == 1) 
            ind <- dum
        else ind <- c(ind, dum)
    }
    ind <- rep(ind, 2)
    if (x[2] == 2) 
        data.al <- data.frame(pop = rep(data[, 1], 2), ind = ind, 
            al = c(firstal, secal))
    else data.al <- data.frame(pop = rep(data[, 1], 2), ind = ind, 
        al = rbind(firstal, secal))
    names(data.al)[-c(1:2)] <- names(data)[-1]
    return(data.al)
}
