# Visualise stan output

# Most code can only be run with access to original trial data

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(rstan)
library(patchwork)

source("ab_model/R/antibody_model.R")
source("ab_model/R/stan_utils.R")

# Timing of first peak
first_peak_t <- 84
# Time of first boost (after first peak)
t_boost1 <- 365
t_boost2 <- 365 * 2

# Load processed trial data
# Loads the following objects:
# - data - complete dataset
# - stan_data - dataset to fit in STAN
# - fit - output of model fitting in STAN

# Fit processing
fit_summary <- summary(fit)$summary |>
  data.frame() |>
  rownames_to_column(var = "par")

### Individual fits ############################################################

# Individual model predictions
indiv_pred_pd <- list()
for(i in 1:length(unique(data$individual))){
  indiv_fit <- fit_summary |>
    filter(grepl(paste0("[[]", i, "[]]"), par))
  
  td <- filter(data, individual == i)
  ts <- 0:max(td$t)
  phase <- data.frame(t = ts) |>
    left_join(td, by = "t") |>
    fill(phase, .direction = "down") |>
    pull(phase)
  
  indiv_pred_pd[[i]] <-
    data.frame(t = ts,
               individual = i,
               ab =
                 antibody_titre(ts,
                                phase,
                                peak1 = stan_data$peak1[i],
                                peak2 = stan_data$peak2[i],
                                peak3 = stan_data$peak3[i],
                                duration1 = indiv_fit[1, "X50."],
                                duration2 = indiv_fit[2, "X50."],
                                rho1 = indiv_fit[5, "X50."],
                                rho2 = indiv_fit[6, "X50."]),
               abl1 =
                 antibody_titre(ts,
                                phase,
                                peak1 = stan_data$peak1[i],
                                peak2 = stan_data$peak2[i],
                                peak3 = stan_data$peak3[i],
                                duration1 = indiv_fit[1, "X2.5."],
                                duration2 = indiv_fit[2, "X2.5."],
                                rho1 = indiv_fit[5, "X97.5."],
                                rho2 = indiv_fit[6, "X97.5."]),
               abl2 =
                 antibody_titre(ts, 
                                phase,
                                peak1 = stan_data$peak1[i],
                                peak2 = stan_data$peak2[i],
                                peak3 = stan_data$peak3[i],
                                duration1 = indiv_fit[1, "X25."],
                                duration2 = indiv_fit[2, "X25."],
                                rho1 = indiv_fit[5, "X75."],
                                rho2 = indiv_fit[6, "X75."]),
               abu1 =
                 antibody_titre(ts, 
                                phase,
                                peak1 = stan_data$peak1[i],
                                peak2 = stan_data$peak2[i],
                                peak3 = stan_data$peak3[i],
                                duration1 = indiv_fit[1, "X97.5."],
                                duration2 = indiv_fit[2, "X97.5."],
                                rho1 = indiv_fit[5, "X2.5."],
                                rho2 = indiv_fit[6, "X2.5."]),
               abu2 =
                 antibody_titre(ts, 
                                phase,
                                peak1 = stan_data$peak1[i],
                                peak2 = stan_data$peak2[i],
                                peak3 = stan_data$peak3[i],
                                duration1 = indiv_fit[1, "X75."],
                                duration2 = indiv_fit[2, "X75."],
                                rho1 = indiv_fit[5, "X25."],
                                rho2 = indiv_fit[6, "X25."]))
  
}
indiv_pred_pd <- dplyr::bind_rows(indiv_pred_pd)
# Individual observations
indiv_dat_pd <- data.frame(
  individual = stan_data$individual,
  t = stan_data$t,
  ab = stan_data$ab
)
# Individual plot
indiv_pred_plot <- ggplot() +
  geom_ribbon(data = indiv_pred_pd, aes(x = t, ymin = abl2, ymax = abu2), fill = "red", alpha = 0.4) + 
  geom_ribbon(data = indiv_pred_pd, aes(x = t, ymin = abl1, ymax = abu1), fill = "red", alpha = 0.2) +
  geom_line(data = indiv_pred_pd, aes(x = t, y = ab), col = "darkred") +
  geom_point(data = indiv_dat_pd, aes(x = t, y = ab), col = "darkslategray4", size = 0.8) +
  facet_wrap(~ individual, scales = "free_y") + 
  ylab("Anti-NANP IgG (EU/mL)") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 365 * 5, 365), labels = 0:5, name = "Time (years)") +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 6))

################################################################################

### Priors for dose-response model #############################################
# Sample global parameters
global_par_names <- c("global_duration1", 
                      "global_duration2", 
                      "global_rho1",
                      "global_rho2")
samp <- extract(fit, pars = global_par_names)

# Fit normal distributions to posteriors
prior_inputs1 <- lapply(samp[1:2], function(x){
  MASS::fitdistr(x, densfun = "lognormal")$estimate
})
prior_inputs2 <- lapply(samp[3:4], function(x){
  MASS::fitdistr(x, densfun = "normal")$estimate
})
prior_inputs <- c(prior_inputs1, prior_inputs2)
names(prior_inputs) <- global_par_names
prior_inputs_df <- bind_rows(prior_inputs) |>
  mutate(par = global_par_names) |>
  select(par, mean, sd, meanlog, sdlog) |>
  rename(prior_mean = mean,
         prior_sd = sd,
         prior_meanlog = meanlog,
         prior_sdlog = sdlog)
#write.csv(prior_inputs_df, "ab_model/output/ab_prior_parameters.csv")

# Plot posterior and prior to check:
samp_pd <- list()
dist_pd <- list()
for(par in global_par_names){
  samp_pd[[par]] <- data.frame(par = par, y = samp[[par]])
  
  m <- prior_inputs[[par]][1]
  sd <- prior_inputs[[par]][2]
  
  f <- dnorm
  if(grepl("duration", par)){
    f <- dlnorm
  }
  
  dist_pd[[par]] <- data.frame(
    par = par, 
    x = seq(
      min(samp_pd[[par]]$y),
      max(samp_pd[[par]]$y),
      length.out = 1000
    )
  ) |>
    mutate(d = f(x,
                 mean = m,
                 sd = sd
    )
    )
}
samp_pd <- bind_rows(samp_pd)
dist_pd <- bind_rows(dist_pd)

ggplot() + 
  geom_histogram(data = samp_pd, aes(x = y, y = ..density..), bins = 100, fill = "grey40") +
  geom_line(data = dist_pd, aes(x = x, y = d), col = "deeppink") +
  facet_wrap(~ par, scales = "free", nrow = 2) + 
  theme_bw()
################################################################################

### Global Fit #################################################################
global_fit <- fit_summary |>
  filter(grepl("global", par))

# Following page 5 of: Frequentist performances of Bayesian prediction intervals for
# random-effects meta-analysis
qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)

# Summarise the geometric mean and quantiles
summary_geom <- function(x){
  x <- x[is.finite(x)]
  c(
    exp(quantile(log(x), c(0.025, 0.25))),
    mean = exp(mean(log(x))),
    exp(quantile(log(x), c(0.75, 0.975)))
  )
}

ab_global_par <- data.frame(
  peak1 = summary_geom(stan_data$peak1),
  peak2 =  summary_geom(stan_data$peak2),
  peak3 =  summary_geom(stan_data$peak3),
  duration1 = qlnorm(qs,
                     log(global_fit[global_fit$par == "global_duration1", "mean"]),
                     global_fit[global_fit$par == "global_duration1_sd", "mean"]),
  duration2 = qlnorm(qs,
                     log(global_fit[global_fit$par == "global_duration2", "mean"]),
                     global_fit[global_fit$par == "global_duration2_sd", "mean"]),
  rho1 = stan_inv_logit(qnorm(qs,
                              global_fit[global_fit$par == "global_rho1", "mean"],
                              global_fit[global_fit$par == "global_rho1_sd", "mean"])),
  rho2 = stan_inv_logit(qnorm(qs,
                              global_fit[global_fit$par == "global_rho2", "mean"],
                              global_fit[global_fit$par == "global_rho2_sd", "mean"]))
)

# Define t and phases for pre boost 2 and post boost 2 groups
t1 <- 1:(365*5)
phase1 <- ifelse(t1 >= t_boost1, 2, 1)
group1 <- rep("Pre boost 2", length(t1))

t2 <- t1[t1 >= t_boost2]
phase2 <- rep(3, length(t2))
group2 <- rep("Post boost 2", length(t2))

t <- c(t1, t2)
phase <- c(phase1, phase2)
group <- c(group1, group2)

# Predict across quantiles
global_pred <- apply(ab_global_par, 1, function(x, t, phase){
  antibody_titre(t, 
                 phase,
                 peak1 = x[1],
                 peak2 = x[2],
                 peak3 = x[3],
                 duration1 = x[4],
                 duration2 = x[5],
                 rho1 = x[6],
                 rho2 = x[7])
}, t = t, phase = phase) |>
  data.frame()

colnames(global_pred) <- c("abll", "abl", "ab", "abu", "abuu")
global_pred$t <- t + first_peak_t
global_pred$group <- group

# Observations
global_dat_pd <- data |>
  mutate(scenario = ifelse(phase == 3, "Post boost 2", "Pre boost 2")) |>
  group_by(t, scenario) |>
  summarise(y = exp(mean(log(ab))),
            ymin = exp(quantile(log(ab), 0.025)),
            ymax = exp(quantile(log(ab), 0.975))) |>
  mutate(t = t + first_peak_t) |>
  ungroup()

# Plotting data for pre boost 2
global_pred_1B <- global_pred |>
  filter(group == "Pre boost 2")
global_dat_pd_1B <- global_dat_pd |>
  filter(scenario == "Pre boost 2")

global_pred_plot_1B <-
  ggplot() + 
  geom_ribbon(data = global_pred_1B, aes(x = t, ymin = abll, ymax = abuu), alpha = 0.2, fill = "darkorchid3") + 
  geom_ribbon(data = global_pred_1B, aes(x = t, ymin = abl, ymax = abu), alpha = 0.4, fill = "darkorchid3") +
  geom_line(data = global_pred_1B, aes(x = t, y = ab), lwd = 1, colour = "darkorchid4") +
  geom_linerange(data = global_dat_pd_1B, aes(x = t, ymin = ymin, ymax = ymax)) +
  geom_point(data = global_dat_pd_1B, aes(x = t, y = y)) +
  scale_fill_discrete(guide = "none") +
  scale_colour_discrete(guide = "none") +
  scale_y_log10(limits = c(0.1, NA)) + 
  scale_x_continuous(breaks = seq(0, 365 * 5, 365), labels = 0:5, name = "Time (years)") +
  ylab("Anti-NANP IgG (EU/mL)") +
  theme_bw()

## Plotting data for post boost 2
global_pred_2B <- global_pred |>
  filter(group == "Post boost 2" | (group == "Pre boost 2" & t < t_boost2 + first_peak_t))
global_dat_pd_2B <- global_dat_pd |>
  filter(scenario == "Post boost 2" | (scenario == "Pre boost 2" & t < t_boost2 + first_peak_t))

global_pred_plot_2B <-
  ggplot() + 
  geom_ribbon(data = global_pred_2B, aes(x = t, ymin = abll, ymax = abuu), alpha = 0.2, fill = "springgreen3") + 
  geom_ribbon(data = global_pred_2B, aes(x = t, ymin = abl, ymax = abu), alpha = 0.4, fill = "springgreen3") +
  geom_line(data = global_pred_2B, aes(x = t, y = ab), lwd = 1, col = "springgreen4") +
  geom_linerange(data = global_dat_pd_2B, aes(x = t, ymin = ymin, ymax = ymax)) +
  geom_point(data = global_dat_pd_2B, aes(x = t, y = y)) +
  scale_fill_discrete(guide = "none") +
  scale_colour_discrete(guide = "none") +
  scale_y_log10(limits = c(0.1, NA)) + 
  scale_x_continuous(breaks = seq(0, 365 * 5, 365), labels = 0:5, name = "Time (years)") +
  ylab("Anti-NANP IgG (EU/mL)") +
  theme_bw()

global_pred_plot <- global_pred_plot_1B | global_pred_plot_2B

#write.csv(ab_global_par, "ab_model/output/antibody_global_fit.csv")
################################################################################

### Extract parameters for malariasimulation ###################################
logpeak1 <- log(stan_data$peak1[stan_data$peak1 < Inf])
logpeak2 <- log(stan_data$peak2[stan_data$peak2 < Inf])
logpeak3 <- log(stan_data$peak3[stan_data$peak3 < Inf])

r21_malarisimulation_parameters <- data.frame(
  par = c("r21_cs",
          "r21_cs_boost",
          "r21_cs_boost2",
          "r21_rho",
          "r21_rho_boost",
          "r21_ds",
          "r21_dl"),
  mu = c(mean(logpeak1),
         mean(logpeak2),
         mean(logpeak3),
         global_fit[global_fit$par == "global_rho1", "mean"],
         global_fit[global_fit$par == "global_rho2", "mean"],
         log(global_fit[global_fit$par == "global_duration1", "mean"]),
         log(global_fit[global_fit$par == "global_duration2", "mean"])),
  sd = c(sd(logpeak1),
         sd(logpeak2),
         sd(logpeak3),
         global_fit[global_fit$par == "global_rho1_sd", "mean"],
         global_fit[global_fit$par == "global_rho2_sd", "mean"],
         global_fit[global_fit$par == "global_duration1_sd", "mean"],
         global_fit[global_fit$par == "global_duration2_sd", "mean"])
)

# To check: The following transformations take place within malariasimulation:
rtss_p <- malariasimulation::get_parameters()
# Peaks
mean(exp(rnorm(100000, rtss_p$rtss_cs[1], rtss_p$rtss_cs[2])))
mean(exp(rnorm(100000, rtss_p$rtss_cs_boost[1], rtss_p$rtss_cs_boost[2])))
mean(exp(rnorm(100000, r21_malarisimulation_parameters[1,2], r21_malarisimulation_parameters[1,3])))
mean(exp(rnorm(100000, r21_malarisimulation_parameters[2,2], r21_malarisimulation_parameters[2,3])))
# Rhos
mean(malariasimulation:::invlogit(rnorm(100000, rtss_p$rtss_rho[1], rtss_p$rtss_rho[2])))
mean(malariasimulation:::invlogit(rnorm(100000, rtss_p$rtss_rho_boost[1], rtss_p$rtss_rho_boost[2])))
mean(malariasimulation:::invlogit(rnorm(100000, r21_malarisimulation_parameters[3,2], r21_malarisimulation_parameters[3,3])))
mean(malariasimulation:::invlogit(rnorm(100000, r21_malarisimulation_parameters[4,2], r21_malarisimulation_parameters[4,3])))
# Durations
mean(exp(rnorm(100000, rtss_p$rtss_ds[1], rtss_p$rtss_ds[2])))
mean(exp(rnorm(100000, rtss_p$rtss_dl[1], rtss_p$rtss_dl[2])))
mean(exp(rnorm(100000, r21_malarisimulation_parameters[5,2], r21_malarisimulation_parameters[5,3])))
mean(exp(rnorm(100000, r21_malarisimulation_parameters[6,2], r21_malarisimulation_parameters[6,3])))

#write.csv(r21_malarisimulation_parameters, "ab_model/output/r21_malarisimulation_parameters.csv", row.names = FALSE)

### Get draws for modelling uncertainty in malariasimulation
set.seed(13022023)
n_draws <- 50
draws <- as.matrix(fit) |>
  as.data.frame() |>
  select(contains("global"))
draws <- draws[sample(1:nrow(draws), n_draws),] |>
  mutate(draw = 1:n_draws) |>
  mutate(global_duration1 = log(global_duration1),
         global_duration2 = log(global_duration2))

r21_malarisimulation_parameter_draws <- data.frame(
  par = rep(
    c("r21_cs",
      "r21_cs_boost",
      "r21_cs_boost2",
      "r21_rho",
      "r21_rho_boost",
      "r21_ds",
      "r21_dl"),
    each = n_draws),
  draw = rep(1:n_draws, 7),
  mu = c(
    rep(mean(logpeak1), n_draws),
    rep(mean(logpeak2), n_draws),
    rep(mean(logpeak3), n_draws),
    draws$global_rho1,
    draws$global_rho2,
    draws$global_duration1,
    draws$global_duration2
  ),
  sd = c(
    rep(sd(logpeak1), n_draws),
    rep(sd(logpeak2), n_draws),
    rep(sd(logpeak3), n_draws),
    draws$global_rho1_sd,
    draws$global_rho2_sd,
    draws$global_duration1_sd,
    draws$global_duration2_sd
  )
)

# Check draws are in line with best fit parameters (means and sds)
ggplot() + 
  geom_histogram(data = r21_malarisimulation_parameter_draws, aes(x = mu)) + 
  geom_vline(data = r21_malarisimulation_parameters, aes(xintercept = mu), col = "deeppink") +
  facet_wrap(~par, scale = "free")

ggplot() + 
  geom_histogram(data = r21_malarisimulation_parameter_draws, aes(x = sd)) + 
  geom_vline(data = r21_malarisimulation_parameters, aes(xintercept = sd), col = "deeppink") +
  facet_wrap(~par, scale = "free")

#write.csv(r21_malarisimulation_parameter_draws, "ab_model/output/r21_malarisimulation_parameter_draws.csv", row.names = FALSE)
