## move RAD-seq data after demultiplexing
## assumes you start in the folder where the fastq are

library(magrittr)
library(openxlsx)
dat.metaFile = 'https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/AAA.DATA/RADSpecimensFloragenexLiveCopy.xlsx'
dat.dir = './'
seqPatt = 'fastq.gz'
stripPost = '_R1'
fileDrops = 'FGXCONTROL'
matchCol = 'Extraction_Tube_no'
outPath = '/mnt/LAB_DATA/RAD-data/'
outPrefix = 'AAA.'

if(substr(dat.metaFile, nchar(dat.metaFile)-2, nchar(dat.metaFile)) == 'csv') {
  dat.meta <- read.csv(dat.metaFile)
} else dat.meta <- read.xlsx(dat.metaFile)
dat.filesOrig <- dir(dat.dir, patt = seqPatt)
dat.filesOrig <- grep('FGXCONTROL', dat.filesOrig, invert = T, value = T)
dat.files <- sapply(strsplit(dat.filesOrig, stripPost, fixed = T), '[', 1)

dat.files.row <-
  # sapply(dat.files, grep, x = dat.meta[[matchCol]])
   match(dat.files, dat.meta[[matchCol]])
if(class(dat.files.row) != 'integer') stop('Check errors: dat.files.row')
if(any(is.na(dat.files.row))) stop('NA in dat.files.row')

dat.files.species <- dat.meta$Species[dat.files.row] %>% as.character
dat.files.genus <-
  dat.files.species %>%
  as.character %>%
  strsplit(split = ' ', fixed = T) %>%
  sapply(FUN = '[', 1) %>%
  toupper

out <- paste(
  'mv ',
  dat.dir, dat.filesOrig, ' ',
  outPath, outPrefix, dat.files.genus, '/', dat.filesOrig, ' ',
  '# ', dat.files.species,
  sep = '')
