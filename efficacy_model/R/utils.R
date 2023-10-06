logit <- function(x){
  x / (1 + x)
}

revlogit <- function(x){
  x / (1 - x)
}

draw_lognormal <- function(n, mean, sd){
  exp(rnorm(n, mean, sd))
}

draw_logit <- function(n, mean, sd){
  draw <- exp(rnorm(n, mean, sd))
  logit(draw)
}