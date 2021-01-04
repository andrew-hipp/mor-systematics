ic.barplot <- function(x, ...) {
## x is output from informationCriterion
  windows(10,3.5)
  par(mar  = c(5, 14, 4, 2) + 0.1)
  ticks = barplot(x$AICcwi, horiz=T, xlab = 'AICc weight', axes = 1, ...)
  axis(2, labels  = (x$names), at = ticks, las  = 1, tick= F, cex.axis = 0.6)
  }
