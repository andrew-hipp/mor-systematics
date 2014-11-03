## functions for manipulating and querying a SP2 checklist

subset.checklist <- function(x, rank = "Species", usage = "accepted", terms = "Carex|Cymophyllus|Kobresia|Uncinia|Schoenoxiphium") {
  x$Rank == rank & x$Usage == usage & 1:dim(x)[1] %in% grep(terms, x$Term.name)
  }