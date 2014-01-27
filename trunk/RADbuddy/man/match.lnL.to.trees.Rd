\name{match.lnL.to.trees}
\alias{match.lnL.to.trees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
match.lnL.to.trees(directory = "getwd()", lnLprefix = "RAxML_info.", lnLsuffix = ".lnL", treeIndexFile = "tree.index.lines.txt", locus.names = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{
%%     ~~Describe \code{directory} here~~
}
  \item{lnLprefix}{
%%     ~~Describe \code{lnLprefix} here~~
}
  \item{lnLsuffix}{
%%     ~~Describe \code{lnLsuffix} here~~
}
  \item{treeIndexFile}{
%%     ~~Describe \code{treeIndexFile} here~~
}
  \item{locus.names}{
%%     ~~Describe \code{locus.names} here~~
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
function (directory = "getwd()", lnLprefix = "RAxML_info.", lnLsuffix = ".lnL", 
    treeIndexFile = "tree.index.lines.txt", locus.names = NULL, 
    ...) 
{
    treeIndex <- read.delim(paste(directory, "/", treeIndexFile, 
        sep = ""), as.is = T, header = F, row.names = 1)
    if (is.null(locus.names)) 
        locus.names <- row.names(treeIndex)
    logfile = file(format(Sys.time(), "match.lnL.\%Y-\%m-\%d.log.txt"), 
        open = "a")
    lnL.list <- lapply(paste(directory, "/", lnLprefix, locus.names, 
        lnLsuffix, sep = ""), get.raxml.treeLikelihoods)
    close(logfile)
    names(lnL.list) <- locus.names
    raxml.worked <- names(which(lnL.list != "FAIL"))
    lnL.list <- lnL.list[raxml.worked]
    treeIndex <- treeIndex[raxml.worked, ]
    locus.names <- raxml.worked
    out.mat <- matrix(NA, nrow = length(locus.names), ncol = dim(treeIndex)[2], 
        dimnames = list(locus.names, NULL))
    for (i in locus.names) {
        error.out = try(names(lnL.list[[i]]) <- try(unique(as.character(treeIndex[i, 
            ]))))
        if (class(error.out) == "try-error") 
            message(paste("assigning names to locus", i, "failed"))
        else out.mat[i, ] <- lnL.list[[i]][as.character(treeIndex[i, 
            ])]
    }
    out.mat <- out.mat[-c(which(apply(out.mat, 1, function(x) any(is.na(x))))), 
        ]
    class(out.mat) <- "partitionedRAD"
    return(out.mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
