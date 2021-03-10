plot.strObj <- function(x, species = NULL,
                        clusterNames = NULL, addSpBar = FALSE, ...){
# species is a vector, one item per row of plot.strObj
# clusterNames is an optional named vector of form c(cluster1 = 'xx', ...)
  require(ggplot2)
  require(data.table)
  if(!identical(NULL, clusterNames)) {
    names(x)[match(names(clusterNames), names(x))] <- clusterNames
  }
  if(!identical(NULL, species)) {
    x <- cbind(x, Species = species)
  }
  xm <- melt(x)
  xm$sample <- factor(xm$sample, levels = unique(xm$sample))
  names(xm) <- c('Sample', 'Species', 'Cluster', 'Probability')
  out <- ggplot(xm, aes(x=Sample,y=Probability,fill=Cluster))
  out <- out + geom_bar(stat="identity", position="stack")
  out <- out + theme_classic()
  out <- out + scale_fill_brewer(palette="Spectral")
  out <- out + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
  # if(spFacet) out <- out + facet_grid(~Species)
  if(addSpBar) {
    out <- out +
      geom_segment(
        data = x,
        mapping = aes(x = min(Species), xend = max(Species),
                      y = 1.1, yend = 1.1, fill = Species), #deleted colour = Species
        width = 1
      ) # close geom_bar
  } # close addPopBar
  return(out)
}
