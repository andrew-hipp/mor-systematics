\name{lengths.report}
\alias{lengths.report}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lengths.report(dat, numtodo = 10, reportInterval = 2000, high.mem = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{numtodo}{
%%     ~~Describe \code{numtodo} here~~
}
  \item{reportInterval}{
%%     ~~Describe \code{reportInterval} here~~
}
  \item{high.mem}{
%%     ~~Describe \code{high.mem} here~~
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
function (dat, numtodo = 10, reportInterval = 2000, high.mem = TRUE) 
{
    if (class(dat) != "pyRAD.loci") 
        stop("This function runs on a pyRAD data object")
    last.lines <- dat$cons - 1
    num.loci <- length(last.lines)
    datSeqs <- as.character(dat$seqs)
    if (high.mem) 
        block.lengths <- sapply(datSeqs[last.lines][1:ifelse(numtodo < 
            1, num.loci, numtodo)], function(x) nchar(as.character(x)))
    else {
        block.lengths = integer(0)
        for (i in 1:num.loci) block.lengths = c(block.lengths, 
            nchar(as.character(datSeqs[last.lines[i]])))
        if (i/reportInterval - i\%/\%reportInterval == 0) {
            message(paste("...", i, "of", num.loci, "-- Estimated time remaining =", 
                ((Sys.time() - start.time)/i) * (num.loci - i), 
                attr(Sys.time() - start.time, "units")))
        }
    }
    names(block.lengths) <- dat$locus.index[last.lines]
    return(block.lengths)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
