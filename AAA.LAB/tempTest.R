## to test Q. regarding polarity of contrasts

makeRat = function (x = div, y = rates)
{
    do = sample(c(0, 1), dim(x)[1], replace = T)
    for (i in 1:dim(x)[1]) {
        if (!exists('out'))
            out <- numeric(0)
        if (do[i] == 0)
            out <- c(out, log(rates[i, 1]/rates[i, 2]), log(div[i,
                1]/div[i, 2]))
        if (do[i] == 1)
            out <- c(out, log(rates[i, 2]/rates[i, 1]), log(div[i,
                2]/div[i, 1]))
    }
    out <- matrix(out, dim(x)[1], byrow = T)
    out <- data.frame(out)
    names(out) <- c('x','y')
    out
    }
