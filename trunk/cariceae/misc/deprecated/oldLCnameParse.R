  ac.punctMat <- cbind(comma = regexpr(",", ac), other =regexpr("[. ]", ac))
  ac.punctMat[ac.punctMat < 0] <- 999
  ac.whichOne <- apply(ac.punctMat, 1, function(x) which(x == min(x))) # 1 = comma first, 2 = space or period first
  ac.whichOne[sapply(ac.whichOne, length) > 1] <- 3 # 3 = one-word only
  ac.splitted <- strsplit(ac, "[. ,]")
  for(i in which(!is.na(ac))) {
    ac.temp <- ac.splitted[[i]][nchar(ac.splitted[[i]]) > 1]
	if(ac.whichOne[i] == 1) pcln[i] <- ac.temp[1]
	if(ac.whichOne[i] == 2) pcln[i] <- ac.temp[2]
	if(ac.whichOne[i] == 3) pcln[i] <- ac[i]
	}
