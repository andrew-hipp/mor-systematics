\name{consensus.pyRAD}
\alias{consensus.pyRAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
consensus.pyRAD(pyIn, from = NA, to = NA, fastaNames = T, writeFile = "rads.con.1_100.txt", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pyIn}{
%%     ~~Describe \code{pyIn} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{fastaNames}{
%%     ~~Describe \code{fastaNames} here~~
}
  \item{writeFile}{
%%     ~~Describe \code{writeFile} here~~
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
function (pyIn, from = NA, to = NA, fastaNames = T, writeFile = "rads.con.1_100.txt", 
    ...) 
{
    if (class(pyIn) != "pyRAD.loci") 
        stop("pyRAD input required, from read.pyRAD")
    require(Biostrings)
    allLoci <- unique(as.character(pyIn$locus.index))
    allLoci <- allLoci[!allLoci == ""]
    if (!is.na(from)) 
        allLoci <- allLoci[from:to]
    seqs <- as.character(pyIn$seqs)
    loc.index <- as.character(pyIn$locus.index)
    out <- character(0)
    for (i in allLoci) {
        message(i)
        out <- c(out, consensusString(DNAStringSet(gsub("-", 
            "N", seqs[loc.index == i])), ...))
    }
    if (fastaNames) 
        allLoci <- paste(">", allLoci, sep = "")
    names(out) <- allLoci
    if (!is.na(writeFile)) 
        write.table(out, writeFile, sep = "\n", quote = F, col.names = F)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
