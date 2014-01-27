\name{rank.partitionedRAD}
\alias{rank.partitionedRAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rank.partitionedRAD(radMat, criterion = c("lnL.threshold", "percentile"), minTrees = 10, min.overall.diff.lnL = 1.5, threshold.lnL = 0, discardDoubleCounts = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{radMat}{
%%     ~~Describe \code{radMat} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{minTrees}{
%%     ~~Describe \code{minTrees} here~~
}
  \item{min.overall.diff.lnL}{
%%     ~~Describe \code{min.overall.diff.lnL} here~~
}
  \item{threshold.lnL}{
%%     ~~Describe \code{threshold.lnL} here~~
}
  \item{discardDoubleCounts}{
%%     ~~Describe \code{discardDoubleCounts} here~~
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
function (radMat, criterion = c("lnL.threshold", "percentile"), 
    minTrees = 10, min.overall.diff.lnL = 1.5, threshold.lnL = 0, 
    discardDoubleCounts = TRUE) 
{
    radMat.preLnLDiff <- radMat.preMinTrees <- NA
    if (minTrees > 1) {
        nTrees <- apply(radMat, 1, function(x) length(unique(x)))
        radMat.preMinTrees <- radMat
        radMat <- radMat[nTrees >= minTrees, ]
    }
    if (min.overall.diff.lnL > 0) {
        lnL.diff <- apply(radMat, 1, function(x) abs(diff(range(x))))
        radMat.preLnLDiff <- radMat
        radMat <- radMat[(lnL.diff >= min.overall.diff.lnL), 
            ]
    }
    if (criterion[1] == "lnL.threshold") {
        bestMat <- t(apply(radMat, 1, function(x) abs(x - max(x)) <= 
            threshold.lnL))
        worstMat <- t(apply(radMat, 1, function(x) abs(x - min(x)) <= 
            threshold.lnL))
        doubleCountMat <- bestMat & worstMat
    }
    if (discardDoubleCounts) {
        doubleCounts <- apply(doubleCountMat, 1, sum) > 0
        bestMat <- bestMat[!doubleCounts, ]
        worstMat <- worstMat[!doubleCounts, ]
        radMat.preDoubleCounts <- radMat
        radMat <- radMat[!doubleCounts, ]
    }
    out <- list(bestMat = bestMat, worstMat = worstMat, doubleCountMat = doubleCountMat, 
        radMat = radMat, radMat.preLnLDiff = radMat.preLnLDiff, 
        radMat.preMinTrees = radMat.preMinTrees, radMat.preDoubleCounts = radMat.preDoubleCounts, 
        params = c(minTrees, min.overall.diff.lnL, threshold.lnL, 
            discardDoubleCounts))
    class(out) <- "rankedPartitionedRAD"
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
