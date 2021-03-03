cbind2 <- function (x, fill = NA)
{
    rows <-
    unique(
      unlist(
        lapply(x, row.names)))
    columns <- sum(sapply(x, function(y) dim(y)[2]))
    out <- matrix(fill, length(rows), columns, dimnames = list(rows,
        NULL))
    index = 1
    for (i in 1:length(x)) {
        out[row.names(x[[i]]), index:(index + dim(x[[i]])[2] - 1)] <- x[[i]]
        index = index + dim(x[[i]])[2]
    }
    out
}
