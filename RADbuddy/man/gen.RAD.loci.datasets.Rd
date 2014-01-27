\name{gen.RAD.loci.datasets}
\alias{gen.RAD.loci.datasets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gen.RAD.loci.datasets(rads, trees = "none", loci = "all", taxa = "all", minTaxa = 4, onlyVariable = TRUE, fileBase = format(Sys.time(), "rads.%Y-%m-%d"), splitInto = 8, raxPath = "~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-AVX", header = "#!/bin/sh")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rads}{
%%     ~~Describe \code{rads} here~~
}
  \item{trees}{
%%     ~~Describe \code{trees} here~~
}
  \item{loci}{
%%     ~~Describe \code{loci} here~~
}
  \item{taxa}{
%%     ~~Describe \code{taxa} here~~
}
  \item{minTaxa}{
%%     ~~Describe \code{minTaxa} here~~
}
  \item{onlyVariable}{
%%     ~~Describe \code{onlyVariable} here~~
}
  \item{fileBase}{
%%     ~~Describe \code{fileBase} here~~
}
  \item{splitInto}{
%%     ~~Describe \code{splitInto} here~~
}
  \item{raxPath}{
%%     ~~Describe \code{raxPath} here~~
}
  \item{header}{
%%     ~~Describe \code{header} here~~
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
function (rads, trees = "none", loci = "all", taxa = "all", minTaxa = 4, 
    onlyVariable = TRUE, fileBase = format(Sys.time(), "rads.\%Y-\%m-\%d"), 
    splitInto = 8, raxPath = "~/code/raxml/standard-RAxML-8.0.2/raxmlHPC-AVX", 
    header = "#!/bin/sh") 
{
    if (!paste(fileBase, ".", (0), sep = "") \%in\% dir()) 
        lapply(paste(fileBase, ".", (0:splitInto), sep = ""), 
            dir.create)
    analysisFileOut <- lapply(paste(fileBase, ".0/raxml.batch.", 
        1:splitInto, ".", fileBase, ".sh", sep = ""), file, open = "a")
    for (i in 1:splitInto) cat(header, "\n", file = analysisFileOut[[i]])
    indexFileOut <- file(paste(fileBase, ".0/tree.index.lines.txt", 
        sep = ""), "a")
    if (loci[1] == "all") 
        loci <- unique(rads$locus.index)[unique(rads$locus.index) != 
            ""]
    if (taxa[1] == "all") 
        taxa <- unique(rads$tips)[gsub("/", "", unique(rads$tips), 
            fixed = TRUE) != ""]
    if (trees[1] != "none") 
        taxa <- intersect(taxa, trees[[1]]$tip.label)
    locus.set <- subset.pyRAD.loci(rads, loci, taxa)
    locus.list <- locus.set$DNA[names(which(locus.set$ntaxa >= 
        minTaxa))]
    if (onlyVariable) 
        locus.list <- locus.list[names(which(locus.set$variable))]
    if (trees[1] != "none") 
        tree.vector.matrix <- matrix(NA, nrow = length(locus.list), 
            ncol = length(trees), dimnames = list(names(locus.list), 
                names(trees)))
    batch = 0
    for (i in names(locus.list)) {
        if (batch == splitInto) 
            batch <- 1
        else batch <- batch + 1
        message(paste("Writing", i))
        locus.taxa <- names(locus.list[[i]])[names(locus.list[[i]]) \%in\% 
            taxa]
        datFileOut <- paste(fileBase, ".", batch, "/", i, ".phy", 
            sep = "")
        write.DNAStringSet(locus.list[[i]][locus.taxa], filename = datFileOut)
        if (trees[1] != "none") {
            trees.out <- lapply(trees, drop.tip, tip = trees[[1]]$tip.label[!trees[[1]]$tip.label \%in\% 
                locus.taxa])
            trees.out <- lapply(trees.out, unroot)
            class(trees.out) <- "multiPhylo"
            trees.out <- unique(trees.out)
            message(paste("... kept", length(trees.out), "trees"))
            treeFileOut <- paste(fileBase, ".", batch, "/", i, 
                ".tre", sep = "")
            write.tree(trees.out, file = treeFileOut)
            cat(i, "\t", paste(attr(trees.out, "old.index"), 
                collapse = "\t"), "\n", sep = "", file = indexFileOut)
        }
        analysisLine <- paste(raxPath, "-f G -s", paste("../", 
            datFileOut, sep = ""), "-m GTRGAMMA -z", paste("../", 
            treeFileOut, sep = ""), "-n", paste(i, ".lnL", sep = ""))
        cat(analysisLine, "\n", file = analysisFileOut[[batch]])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
