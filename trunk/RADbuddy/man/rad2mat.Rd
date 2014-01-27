\name{rad2mat}
\alias{rad2mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rad2mat(pyDat, fill.N = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pyDat}{
%%     ~~Describe \code{pyDat} here~~
}
  \item{fill.N}{
%%     ~~Describe \code{fill.N} here~~
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
function (pyDat, fill.N = TRUE) 
{
    if (!"radSummary" \%in\% names(pyDat)) 
        pyDat$radSummary <- summary(pyDat)
    loci <- dimnames(pyDat$radSummary$inds.mat)[[2]]
    inds <- dimnames(pyDat$radSummary$inds.mat)[[1]]
    out <- matrix(NA, length(inds), length(loci), dimnames = list(inds, 
        loci))
    for (i in loci) {
        message(paste("Doing locus", i))
        out[, i] <- pyDat$radSummary$seqs.per.locus[[i]][inds]
        if (fill.N) 
            out[inds[!inds \%in\% pyDat$radSummary$tips.per.locus[[i]]], 
                i] <- paste(rep("N", pyDat$radSummary$locus.lengths[[i]]), 
                collapse = "")
    }
    class(out) <- "rad.mat"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
