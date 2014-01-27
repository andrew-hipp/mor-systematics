\name{write.DNAStringSet}
\alias{write.DNAStringSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.DNAStringSet(x, format = "phylip", padding = 30, filename = "DNAStringSetOut.phy")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{format}{
%%     ~~Describe \code{format} here~~
}
  \item{padding}{
%%     ~~Describe \code{padding} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
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
function (x, format = "phylip", padding = 30, filename = "DNAStringSetOut.phy") 
{
    x.width <- width(x)[1]
    x <- as.character(x)
    for (i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", 
        (padding - nchar(names(x)[i]))), collapse = ""), x[i], 
        sep = "")
    writeLines(c(paste(length(x), x.width), x), filename)
    return(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
