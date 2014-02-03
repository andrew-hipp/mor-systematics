\name{plot.rankedPartitionedRAD}
\alias{plot.rankedPartitionedRAD}
\title{Plot a \code{rankedPartitionedRAD} object}
\description{This function creates a plot with 1-3 panels
			 describing the number of loci supporting, disfavoring,
			 or ambiguous for each tree, plotted against log-likelihood}
\usage{plot.rankedPartitionedRAD(x, fileprefix = NULL, widthScalar = 0.85, panels = c("bestMat", "worstMat", "doubleCountMat"), squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3), primeTreeColor = "red", primeTreeCharacter = 19, filebase = paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."), paste(c("minT", "rangeL", "diffL", "noDoubles"), x$params, collapse = "_", sep = ""), ".pdf", sep = ""), ...)}
\arguments{
  \item{x}{an object of class \code{rankedPartitionedRAD}}
  \item{fileprefix}{prefix of the files to save; as written, you have to provide some prefix
                    to identify the files, but the remainder of the file name defaults to a
					combination of parameters and date}
  \item{widthScalar}{a multiplier on the width of the entire plot to get the plotting panels square
					 or not as you prefer}
  \item{panels}{a character vector describing the order of panels; 
				can be any of "bestMat", "worstMat", "doubleCountMat" in any order}
  \item{squareSize}{the size of the panels in inches; defaults to 5" for a 2-panel plot, 
					3" for a 3-panel plot, only because that works well in my trials}
  \item{primeTreeColor}{the color of the dot representing the first tree in the matrix, which 
						should be your original tree (unpermuted) if you've been following the 
						workflow from \code{genTrees} through this final and exciting plot}
  \item{primeTreeCharacter}{the character of the dot representing the first tree in the matrix}
  \item{filebase}{the base of the file name for exported plots}
  \item{\dots}{additional arguments passed to \code{plot}}
}
\details{This function as written only prints to pdf}
\value{No value. This function is called for its side-effect of printing plots to a file.}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
	\code{\link{rank.partitionedRAD}},
	\code{\link{plot.locus.dist}}
}
\keyword{hplot}