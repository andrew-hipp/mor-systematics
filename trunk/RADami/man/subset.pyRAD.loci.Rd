\name{subset.pyRAD.loci}
\alias{subset.pyRAD.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subset.pyRAD.loci(x, loci, taxa, format = "DNAStringSet", reportInterval = 20, nucVarType = c("verystrict", "strict", "relaxed"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{loci}{
%%     ~~Describe \code{loci} here~~
}
  \item{taxa}{
%%     ~~Describe \code{taxa} here~~
}
  \item{format}{
%%     ~~Describe \code{format} here~~
}
  \item{reportInterval}{
%%     ~~Describe \code{reportInterval} here~~
}
  \item{nucVarType}{
%%     ~~Describe \code{nucVarType} here~~
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
function (x, loci, taxa, format = "DNAStringSet", reportInterval = 20, 
    nucVarType = c("verystrict", "strict", "relaxed"), ...) 
{
    excludedNucs <- switch(nucVarType[1], verystrict = 5:17, 
        strict = 5:17, relaxed = 15:17)
    nucThresh <- switch(nucVarType[1], verystrict = 2, strict = 1, 
        relaxed = 1)
    out <- list(DNA = structure(vector("list", length(loci)), 
        names = loci), variable = structure(logical(length(loci)), 
        names = loci), ntaxa = structure(integer(length(loci)), 
        names = loci))
    inds.vector <- x$tips \%in\% taxa
    counter = 0
    start.time <- Sys.time()
    for (i in loci) {
        counter <- counter + 1
        if (counter/reportInterval - counter\%/\%reportInterval == 
            0) {
            message(paste("... subsetting", counter, "of", length(loci), 
                "-- Estimated time remaining =", round(((Sys.time() - 
                  start.time)/counter) * (length(loci) - counter), 
                  1), attr(Sys.time() - start.time, "units")))
        }
        seq.index <- x$locus.index == i & inds.vector
        out$DNA[[i]] <- DNAStringSet(x$seqs[seq.index])
        names(out$DNA[[i]]) <- x$tips[seq.index]
        out$variable[i] <- any(apply(consensusMatrix(out$DNA[[i]])[-c(excludedNucs), 
            ], 2, function(x) sum(x >= nucThresh) > 1))
        out$ntaxa[i] <- sum(seq.index)
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
