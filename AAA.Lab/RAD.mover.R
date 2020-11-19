## move RAD-seq data after demultiplexing
## assumes you start in the folder where the fastq are

library(magrittr)

dat.metaFile = 'https://www.dropbox.com/s/81ppkw8xqycmxz7/RADSpecimensFloragenexLiveCopy_2020-11-16.csv?raw=1'
dat.dir = '.'
seqPatt = 'fastq.gz'
stripPost = '_R1'
fileDrops = 'FGXCONTROL'

dat.meta <- read.csv(dat.metaFile)
dat.files <- dir(dat.dir, patt = seqPatt)
if(!is.null(stripPost)) {
  dat.files <- sapply(strsplit(dat.files, stripPost, fixed = T), '[', 1)
}
dat.files <- grep('FGXCONTROL', dat.files, invert = T, value = T)
