\name{summary.pyRAD.loci}
\alias{summary.pyRAD.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.pyRAD.loci(object, var.only = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{var.only}{
%%     ~~Describe \code{var.only} here~~
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
function (object, var.only = FALSE, ...) 
{
    message("\nDoing a pyRAD summary")
    reportInterval <- 2000
    locus.names <- as.character(unique(object$locus.index))
    locus.names <- locus.names[locus.names != ""]
    variable.loci <- locus.names[!is.na(object$seqs[object$breaks])]
    if (var.only) 
        locus.names <- variable.loci
    num.loci <- length(locus.names)
    tip.names <- as.character(unique(object$tips[-c(object$breaks, 
        object$cons)]))
    message("Splitting tips by locus name...")
    tips.per.locus <- split(object$tips, object$locus.index)[locus.names]
    seqs.per.locus <- split(object$seqs, object$locus.index)[locus.names]
    num.inds.per.locus <- sapply(tips.per.locus, length)
    inds.mat <- matrix(NA, nrow = length(tip.names), ncol = num.loci, 
        dimnames = list(tip.names, locus.names))
    message("Making tips matrix...")
    start.time <- Sys.time()
    for (i in seq(num.loci)) {
        temp <- try(tip.names \%in\% tips.per.locus[[locus.names[i]]])
        if (class(temp) != "try-error") {
            inds.mat[, locus.names[i]] <- temp
            names(seqs.per.locus[[i]]) <- tips.per.locus[[i]]
        }
        else (message(paste("Error occurred with locus", locus.names[i])))
        if (i/reportInterval - i\%/\%reportInterval == 0) {
            message(paste("...", i, "of", num.loci, "-- Estimated time remaining =", 
                ((Sys.time() - start.time)/i) * (num.loci - i), 
                attr(Sys.time() - start.time, "units")))
        }
    }
    out <- list(num.loci = num.loci, tips.per.locus = tips.per.locus, 
        break.vectors = object$break.vectors, seqs.per.locus = seqs.per.locus, 
        num.inds.per.locus = num.inds.per.locus, variable.loci = variable.loci, 
        inds.mat = inds.mat, locus.lengths = lengths.report(object, 
            0))
    class(out) <- "summary.pyRAD.loci"
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
