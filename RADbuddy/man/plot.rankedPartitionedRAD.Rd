\name{plot.rankedPartitionedRAD}
\alias{plot.rankedPartitionedRAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.rankedPartitionedRAD(x, fileprefix = NULL, widthScalar = 0.85, panels = c("bestMat", "worstMat", "doubleCountMat"), squareSize = switch(as.character(length(panels)), `2` = 5, `3` = 3), primeTreeColor = "red", primeTreeCharacter = 19, filebase = paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."), paste(c("minT", "rangeL", "diffL", "noDoubles"), x$params, collapse = "_", sep = ""), ".pdf", sep = ""), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{fileprefix}{
%%     ~~Describe \code{fileprefix} here~~
}
  \item{widthScalar}{
%%     ~~Describe \code{widthScalar} here~~
}
  \item{panels}{
%%     ~~Describe \code{panels} here~~
}
  \item{squareSize}{
%%     ~~Describe \code{squareSize} here~~
}
  \item{primeTreeColor}{
%%     ~~Describe \code{primeTreeColor} here~~
}
  \item{primeTreeCharacter}{
%%     ~~Describe \code{primeTreeCharacter} here~~
}
  \item{filebase}{
%%     ~~Describe \code{filebase} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, fileprefix = NULL, widthScalar = 0.85, panels = c("bestMat", 
    "worstMat", "doubleCountMat"), squareSize = switch(as.character(length(panels)), 
    `2` = 5, `3` = 3), primeTreeColor = "red", primeTreeCharacter = 19, 
    filebase = paste(format(Sys.time(), "rad.partitioned.\%Y-\%m-\%d."), 
        paste(c("minT", "rangeL", "diffL", "noDoubles"), x$params, 
            collapse = "_", sep = ""), ".pdf", sep = ""), ...) 
{
    if (class(x) != "rankedPartitionedRAD") 
        warning("Not the expected object class; this function may misbehave")
    if (!is.null(fileprefix)) 
        pdf(paste(fileprefix, filebase, sep = "."), width = squareSize * 
            length(panels) * widthScalar, height = squareSize)
    layout(matrix(seq(length(panels)), 1, length(panels)))
    for (i in panels) {
        temp.ylab = switch(i, bestMat = "Number of loci supporting tree", 
            worstMat = "Number of loci disfavoring tree", doubleCountMat = "Loci overlapping, excluded")
        plot(colSums(x$radMat), colSums(x[[i]]), xlab = "Tree log-likelihood, summed over loci", 
            ylab = temp.ylab, ...)
        points(colSums(x$radMat)[1], colSums(x[[i]])[1], pch = primeTreeCharacter, 
            col = primeTreeColor)
    }
    dev.off()
    return("done")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
