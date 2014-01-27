\name{read.pyRAD}
\alias{read.pyRAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.pyRAD(filename, reportInterval = 20000, breakLinesSeparate = FALSE, doSummary = TRUE, makeSeqMat = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{reportInterval}{
%%     ~~Describe \code{reportInterval} here~~
}
  \item{breakLinesSeparate}{
%%     ~~Describe \code{breakLinesSeparate} here~~
}
  \item{doSummary}{
%%     ~~Describe \code{doSummary} here~~
}
  \item{makeSeqMat}{
%%     ~~Describe \code{makeSeqMat} here~~
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
function (filename, reportInterval = 20000, breakLinesSeparate = FALSE, 
    doSummary = TRUE, makeSeqMat = TRUE, ...) 
{
    message("Reading data...")
    dat <- readLines(filename, ...)
    dat.breakLines <- dat.consensusLines <- grep("//", dat, fixed = TRUE)
    dat.breakLines.vector <- dat[dat.breakLines]
    message("Splitting data...")
    dat.split <- strsplit(dat, " {1,100}")
    dat.names <- sapply(dat.split, function(x) x[1])
    dat.seqs <- sapply(dat.split, function(x) x[2])
    dat.firstLocusLines <- c(1, (dat.breakLines[1:(length(dat.breakLines) - 
        1)] + 1))
    if (breakLinesSeparate) {
        dat.lastLocusLines <- dat.breakLines - 2
        dat.consensusLines <- dat.breakLines - 1
    }
    else dat.lastLocusLines <- dat.breakLines - 1
    dat.locus.index <- character(length(dat))
    locusCounter <- 1
    message("Assigning locus number...")
    start.time <- Sys.time()
    for (i in 1:length(dat.firstLocusLines)) {
        dat.locus.index[dat.firstLocusLines[i]:dat.lastLocusLines[i]] <- names(dat.breakLines.vector)[i] <- paste("locus.", 
            i, sep = "")
        if (i/reportInterval - i\%/\%reportInterval == 0) {
            message(paste("...", i, "of", length(dat.firstLocusLines), 
                "-- Estimated time remaining =", round(((Sys.time() - 
                  start.time)/i) * (length(dat.firstLocusLines) - 
                  i), 1), attr(Sys.time() - start.time, "units")))
        }
    }
    out = list(tips = dat.names, seqs = dat.seqs, breaks = dat.breakLines, 
        break.vectors = dat.breakLines.vector, cons = dat.consensusLines, 
        locus.index = dat.locus.index, file.read = filename, 
        timestamp = date())
    class(out) <- "pyRAD.loci"
    if (doSummary) 
        out$radSummary <- summary(out)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
