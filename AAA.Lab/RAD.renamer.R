## rename RAD-seq trees after you run them

dat.metaFile = 'https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/AAA.DATA/RADSpecimensFloragenexLiveCopy_2020-11-16.csv'
if(!exists(treeToDo)) {
  stop('Assign the tree you want to relabel to treeToDo and rerun')
  }

dat.meta <- read.csv(dat.metaFile)
