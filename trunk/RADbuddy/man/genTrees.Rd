\name{genTrees}
\alias{genTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genTrees(x, N = 200, filebase = "trial", method = c("nni", "random"), maxmoves = 2, perms = c(length(nni(x)), as.integer(N - (length(nni(x))))), software = c("raxml", "paup"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{filebase}{
%%     ~~Describe \code{filebase} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{maxmoves}{
%%     ~~Describe \code{maxmoves} here~~
}
  \item{perms}{
%%     ~~Describe \code{perms} here~~
}
  \item{software}{
%%     ~~Describe \code{software} here~~
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
function (x, N = 200, filebase = "trial", method = c("nni", "random"), 
    maxmoves = 2, perms = c(length(nni(x)), as.integer(N - (length(nni(x))))), 
    software = c("raxml", "paup"), ...) 
{
    if (class(x) != "phylo") 
        stop("This function requires a phylo object as its first argument")
    if (method[1] == "nni") {
        for (i in seq(maxmoves)) {
            message(paste("doing maxmoves", i))
            if (i == 1) {
                originalTree <- list(x)
                treeset <- c(originalTree, lapply(nni(x), function(x) x))
            }
            else treeset <- c(treeset, lapply(unique(rNNI(x, 
                i, perms[i] * 1.5)), function(x) x))
            treeset <- lapply(treeset, unroot)
            class(treeset) <- "multiPhylo"
            treeset <- unique(treeset)[1:sum(perms[1:i], 1)]
        }
    }
    else if (method[1] == "random") 
        treeset = c(x, rtree.phylo(x, N, ...))
    class(treeset) <- "multiPhylo"
    if (software[1] == "raxml") {
        browser()
        message("writing raxml")
        write.tree(treeset, file = paste(filebase, ".trees.tre", 
            sep = ""))
    }
    if (software[1] == "raxml") {
        message("RAxML chosen as analysis software. Currently, you just need to run this on your own to get the site likelihoods.Try something like this:\n\n\t/home/andrew/code/raxml/standard-RAxML-7.7.2/raxmlHPC-PTHREADS-SSE3 -f g -T 10 -s d6m10.phy -m GTRGAMMA -z analysis.d6m10/RAxML_bestTree.d6m10.out -n d6m10.phy.reduced.siteLnL")
    }
    return(treeset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
