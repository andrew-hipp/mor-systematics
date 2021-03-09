plot.strObj <- function(x, ...){
  require(ggplot2)
  require(data.table)
  xm <- melt(x)
  names(xm) <- c('Sample', 'Cluster', 'Probability')
  out <- ggplot(xm, aes(x=Sample,y=Probability,fill=Cluster))
  out <- out + geom_bar(stat="identity", position="stack")
  out <- out + theme_classic()
  out <- out + scale_fill_brewer(palette="Spectral")
  out <- out + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
  return(out)
}
