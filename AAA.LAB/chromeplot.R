# look at chrome plot
library(magrittr)
 
ch <- read.delim('chromosomelist.txt', header = F)

x <- read.table('Ficus-globosa_S21.depth.txt.gz', sep = '\t', header = F, strip.white=T)

chVector <- ch$V1[grep('Chromosome', ch$V2)] %>%
  gsub(pattern = '>', replace = '', fixed = T)
png('Ficus-globosa-depths.png', width = 1000, height = 2000)
layout(matrix(1:length(chVector), length(chVector), 1))
for(i in chVector) {
  message(paste('plotting chromosome', i))
  xtemp <- x[x$V1 == i, ]
  plot(xtemp[,2], xtemp[,3], col = ifelse(xtemp[,3] < 20,'red','black'), pch=19, xlab='position', ylab='coverage')
}
dev.off()
