stan_logit <- function(x){
  log(x / (1 - x))
}

stan_inv_logit <- function(x){
  1 / (1 + exp(-x))
}