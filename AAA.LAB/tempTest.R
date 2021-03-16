## to test Q. regarding polarity of contrasts
library(magrittr)

nNodes = 50
nPerms = 1000 ## not actually used yet... JP? shall we simulate the Type I?
multiplier = runif(50, 0, 200)

div <- matrix(abs(rnorm(n = nNodes*2)), nNodes, byrow = T, dimnames = list(NULL, c('x', 'y'))) * multiplier
div <- apply(div, 1, sort) %>% t %>% data.frame(names = c('x', 'y'))
names(div) <- c('x', 'y')
rates <- matrix(abs(rnorm(n = nNodes*2)), nNodes, byrow = T, dimnames = list(NULL, c('x', 'y'))) * multiplier
rates <- apply(rates, 1, sort) %>% t %>% data.frame
names(rates) <- c('x', 'y')

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

empRat.lm <- function(plotIt = TRUE) {
  temp <- makeRat()
  out <- lm(y ~ x + 0, temp)
  if(plotIt) plot(y ~ x, temp, pch = 19)
  abline(out)
  out
}

message('
  to simulate drawing the same nodes over and over,
  randomly selecting direction, run this:

  summary(empRat.lm())')
