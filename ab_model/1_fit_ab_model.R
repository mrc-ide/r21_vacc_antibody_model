# STAN modelling fitting of antibody dynamics

# Most code can only be run with access to original trial data

seed <- 258180
set.seed(seed)

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Function
source("ab_model/R/stan_utils.R")
source("ab_model/R/draw_priors.R")

### Load in processed trial data ###############################################

# Generates the following object:
# data - complete dataset

# Number of individuals
N <- length(unique(data$individual))
# Number of data
R <- nrow(data)

# Isolate peaks. Stan doesn't handle NA, so missing values set to Inf and shouldn't
# be used in the ab titre estimation
peaks <- filter(data, 
                t == 0 | (t == 365 & phase == 2) | (t == t_boost2 & phase == 3)) |>
  complete(individual,  nesting(t, phase), fill = list(ab = Inf)) |>
  arrange(individual, phase)
peak1 <- filter(peaks, phase  == 1)$ab
peak2 <- filter(peaks, phase  == 2)$ab
peak3 <- filter(peaks, phase  == 3)$ab
################################################################################

### Priors #####################################################################
# These get passed to distributions within the STAN programme
priors <- list(
  prior_duration1 = c(100, 0.66),
  prior_duration2 = c(1825, 1),
  prior_logit_rho1 = c(stan_logit(0.5), 1),
  prior_logit_rho2 = c(stan_logit(0.5), 1),
  prior_sigma = c(1, 1),
  prior_re_sd = c(5, 1)
)
################################################################################

### Stan data ##################################################################
stan_data <- c(
  list(
    N = N,
    R = R,
    individual = data$individual,
    ab = data$ab,
    t = data$t,
    phase = data$phase,
    peak1 = peak1,
    peak2 = peak2,
    peak3 = peak3
  ),
  priors
)

################################################################################

### Initial values #############################################################
# Define initial values by drawing from the prior
init = list(
  draw_prior(priors, N),
  draw_prior(priors, N),
  draw_prior(priors, N),
  draw_prior(priors, N)
)
################################################################################

### Model fit ##################################################################
fit <- stan(file = "ab_model/ab_model.stan",
            data = stan_data,
            chains = 4,
            iter = 10000,
            init = init,
            seed = seed,
            control = list(adapt_delta = 0.99))

################################################################################