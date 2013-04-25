## Checking Blast2GO results
## A Hipp, 22 April 2013

## read in the datasets that Elisabeth sent
oak.blasts = lapply(dir('./2013-04-09.data/', full = T), read.delim)
names(oak.blasts) = c('est.others', 'quercus.all', 'ref.seq.RNA')
for(i in names(oak.blasts)) oak.blasts[[i]]$data.source <- i
concat.fields <- c("query", "subject", "subject.accessions", "evalue", "query.start", 
                  "query.end", "subject.start", "subject.end", "bit.score", "Seq..Description", "Seq..Length", 
                  "X.Hits", "min..eValue", "mean.Similarity", "X.GOs", "GOs", "Enzyme.Codes", 'data.source')
oak.blasts.new <- rbind(oak.blasts[[1]][concat.fields], oak.blasts[[2]][concat.fields])
oak.blasts.new <- rbind(oak.blasts.new, oak.blasts[[3]][concat.fields])
oak.blasts <- oak.blasts.new
rm(oak.blasts.new)

unique.loci <- sort(unique(oak.blasts$query))
stats.of.interest <- c('number of hits', 'number of unique descriptions', 'descriptions concatenated', 
                       'number of unique GOs', 'GO terms concatenated', 'QUERCUS', 'EST.OTHERS', 'REF.SEQ.RNA')
locus.summary.mat <- matrix(NA, length(unique.loci), length(stats.of.interest), dimnames = list(unique.loci, stats.of.interest))
for(i in unique.loci) {
  # browser()
  temp <- oak.blasts[oak.blasts$query == i, ]
  locus.summary.mat[i, 'QUERCUS'] <- sum(temp$data.source == 'quercus.all') 
  locus.summary.mat[i, 'EST.OTHERS'] <- sum(temp$data.source == 'est.others') 
  locus.summary.mat[i, 'REF.SEQ.RNA'] <- sum(temp$data.source == 'ref.seq.RNA') 
  locus.summary.mat[i, 1] = dim(temp)[1]
  temp.desc <- unique(temp$Seq..Description)
  missing <- unique(c(grep("NA", temp.desc, fixed = TRUE), grep("N/A", temp.desc, fixed = TRUE)))
  if(length(missing) > 0) temp.desc <- temp.desc[-missing]
  temp.desc <- temp.desc[temp.desc != "0"]
  locus.summary.mat[i, 2] <- length(temp.desc)
  locus.summary.mat[i, 3] <- paste(temp.desc, collapse = "||")
  GOs <- sort(unlist(strsplit(as.character(temp$GOs), "; ")))
  GOs <- GOs[GOs != ""]
  locus.summary.mat[i, 4] <- length(GOs)
  locus.summary.mat[i, 5] <- paste(GOs, collapse = "||")
  }