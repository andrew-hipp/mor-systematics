\name{rad2phy}
\alias{rad2phy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rad2phy(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], outfile = "pyMat.out.phy", padding = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pyDat}{
%%     ~~Describe \code{pyDat} here~~
}
  \item{inds}{
%%     ~~Describe \code{inds} here~~
}
  \item{loci}{
%%     ~~Describe \code{loci} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{padding}{
%%     ~~Describe \code{padding} here~~
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
function (pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], 
    outfile = "pyMat.out.phy", padding = 50) 
{
    if (class(pyDat) != "rad.mat") 
        warning("I'm expecting output from rad.mat")
    outfile = file(outfile, "wt")
    open(outfile)
    cat(paste(length(inds), sum(sapply(pyDat[inds[1], loci], 
        nchar)), "\n"), file = outfile)
    for (i in inds) {
        message(paste("Writing DNA line for individual", i))
        cat(i, file = outfile)
        cat(paste(rep(" ", padding - nchar(i)), collapse = ""), 
            file = outfile)
        cat(paste(pyDat[i, loci], collapse = ""), file = outfile)
        cat("\n", file = outfile)
    }
    close(outfile)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
