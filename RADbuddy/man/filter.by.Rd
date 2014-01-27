\name{filter.by}
\alias{filter.by}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filter.by(dat, taxa, threshold = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{taxa}{
%%     ~~Describe \code{taxa} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function (dat, taxa, threshold = "all") 
{
    if (!class(dat) \%in\% c("summary.pyRAD.loci", "pyRAD.loci")) 
        stop("This function only works with summary.pyRAD.loci datatypes")
    if (class(dat) == "pyRAD.loci") 
        dat.mat <- dat$radSummary$inds.mat[taxa, ]
    else dat.mat <- dat$inds.mat[taxa, ]
    if (threshold == "all") 
        threshold <- length(taxa)
    return(names(which(apply(dat.mat, 2, sum) == threshold)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
