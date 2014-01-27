\name{compare.all.trees}
\alias{compare.all.trees}
\title{Compare all trees in a treeset, and return a matrix (lower diagonal only) indicating which pairwise comparisons are identical}
\usage{
compare.all.trees(treeset, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treeset}{an object of class \code{\multiPhylo} or a list of \code{\phylo} objects}
  \item{\dots}{additional arguments to be passed along to \code{\all.equal.phylo}}
}
\details{Slow, tedious one-on-one comparison of trees, depending on the implementation in \code{\ape}. 
         As written, this function is sensitive to rooting.}
\value{A matrix of \code{\TRUE} and \code{\FALSE} values depending on whether the compared trees are identical in topology}
\author{Andrew Hipp}
\seealso{
  \code{\link{unique.multiPhylo}}
}
