plot.strObj <- function(x, extras = NA, addSpBar = FALSE, ...){
# extras is a data frame, one item per row of plot.strObj
  require(ggplot2)
  require(data.table)
  if(!is.na(extras[1])) {
    x <- cbind(x, extras)
  }
  xm <- melt(x)
  xm$sample <- factor(xm$sample, levels = unique(xm$sample))
  names(xm) <- c('Sample', 'Cluster', 'Probability', names(extras))
  out <- ggplot(xm, aes(x=Sample,y=Probability,fill=Cluster))
  out <- out + geom_bar(stat="identity", position="stack")
  out <- out + theme_classic()
  out <- out + scale_fill_brewer(palette="Spectral")
  out <- out + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
  if(addSpBar) {
    out <- out +
      geom_bar(
        mapping = aes(x = Species, y = 1.1, fill = Species),
        stat = "identity", width = 1
      ) # close geom_bar
  } # close addPopBar
  return(out)
}
