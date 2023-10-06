draw_prior <- function(priors, N){
  list(
    duration1 = rlnorm(N, log(priors$prior_duration1[1]), priors$prior_duration1[2]),
    duration2 = rlnorm(N, log(priors$prior_duration2[1]), priors$prior_duration2[2]),
    logit_rho1 = rnorm(N, priors$prior_logit_rho1[1], priors$prior_logit_rho1[2]),
    logit_rho2 = rnorm(N, priors$prior_logit_rho2[1], priors$prior_logit_rho2[2]),
    sigma = rlnorm(1, log(priors$prior_sigma[1]), priors$prior_sigma[2]),
    global_duration1 = rlnorm(1, log(priors$prior_duration1[1]), priors$prior_duration1[2]),
    global_duration1_sd = rgamma(1, priors$prior_re_sd[1], priors$prior_re_sd[2]),
    global_duration2 = rlnorm(1, log(priors$prior_duration2[1]), priors$prior_duration2[2]),
    global_duration2_sd = rgamma(1, priors$prior_re_sd[1], priors$prior_re_sd[2]),
    global_rho1 = rnorm(1, priors$prior_logit_rho1[1], priors$prior_logit_rho1[2]),
    global_rho1_sd = rgamma(1, priors$prior_re_sd[1], priors$prior_re_sd[2]),
    global_rho2 = rnorm(1, priors$prior_logit_rho2[1], priors$prior_logit_rho2[2]),
    global_rho2_sd = rgamma(1, priors$prior_re_sd[1], priors$prior_re_sd[2])
  )
}