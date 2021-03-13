plot.strObj <- function(x, species = NULL,
                        clusterNames = NULL,
                        addSpBar = FALSE,
                        colorVect = NULL,
                        add.x.lab = FALSE,
                        ...){
# species is a vector, one item per row of plot.strObj
# clusterNames is an optional named vector of form c(cluster1 = 'xx', ...)
  if(!identical(species, NULL)) {
    warning('species vector not implemented; dropped for now')
    species <- NULL
  }
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
  if(!identical(NULL, species)) {
    names(xm) <- c('Sample', 'Species', 'Cluster', 'Probability')
  } else names(xm) <- c('Sample', 'Cluster', 'Probability')
  out <- ggplot(xm, aes(x=Sample,y=Probability,fill=Cluster))
  out <- out + geom_bar(stat="identity", position="stack")
  # out <- out + theme_classic()
  if(identical(colorVect, NULL)) {
    out <- out + scale_fill_brewer(palette="Spectral")
  } else out <- out + scale_fill_manual(values = colorVect)
  if(add.x.lab) {
    out <- out + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
  } else out <- out + theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()
                          )
  if(addSpBar) {
    warning('addSpBar currently not working properly')
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
