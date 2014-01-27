\name{overlap.report}
\alias{overlap.report}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overlap.report(dat, repPattern = "_re", origPattern = "_h")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{repPattern}{
%%     ~~Describe \code{repPattern} here~~
}
  \item{origPattern}{
%%     ~~Describe \code{origPattern} here~~
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
function (dat, repPattern = "_re", origPattern = "_h") 
{
    if (class(dat) == "summary.pyRAD.loci") 
        dat <- dat$inds.mat
    reps <- grep(repPattern, row.names(dat), fixed = TRUE, value = TRUE)
    orig <- gsub(repPattern, origPattern, reps, fixed = TRUE)
    outCols <- c("Original", "Replicate", "Intersection", "Union", 
        "Not intersection", "Overlap proportion", "Original loci replicated", 
        "Increase, original to replicate")
    out <- matrix(NA, nrow = length(orig), ncol = length(outCols), 
        dimnames = list(orig, outCols))
    for (i in 1:length(orig)) {
        message(paste("Doing names", orig[i]))
        out[orig[i], "Original"] <- sum(dat[orig[i], ])
        out[orig[i], "Replicate"] <- sum(dat[reps[i], ])
        out[orig[i], "Intersection"] <- sum(colSums(dat[c(orig[i], 
            reps[i]), ]) == 2, na.rm = T)
        out[orig[i], "Union"] <- sum(colSums(dat[c(orig[i], reps[i]), 
            ]) \%in\% 1:2, na.rm = T)
    }
    out[, "Not intersection"] <- out[, "Union"] - out[, "Intersection"]
    out[, "Overlap proportion"] <- round(out[, "Intersection"]/out[, 
        "Union"], 3)
    out[, "Original loci replicated"] <- round(out[, "Intersection"]/out[, 
        "Original"], 3)
    out[, "Increase, original to replicate"] <- round(out[, "Replicate"]/out[, 
        "Original"], 3)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
