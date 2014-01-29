\name{unique.multiPhylo}
\alias{unique.multiPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unique.multiPhylo(x, incomparables = FALSE, use.edge.length = FALSE, use.tip.label = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{incomparables}{
%%     ~~Describe \code{incomparables} here~~
}
  \item{use.edge.length}{
%%     ~~Describe \code{use.edge.length} here~~
}
  \item{use.tip.label}{
%%     ~~Describe \code{use.tip.label} here~~
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
function (x, incomparables = FALSE, use.edge.length = FALSE, 
    use.tip.label = TRUE, ...) 
{
    n <- length(x)
    keep <- 1L
    old.index <- seq_len(n)
    for (i in 2:n) {
        already.seen <- FALSE
        for (j in keep) {
            if (all.equal(x[[j]], x[[i]], use.edge.length = use.edge.length, 
                use.tip.label = use.tip.label)) {
                already.seen <- TRUE
                old.index[i] <- j
                break
            }
        }
        if (!already.seen) 
            keep <- c(keep, i)
    }
    res <- x[keep]
    attr(res, "old.index") <- old.index
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
