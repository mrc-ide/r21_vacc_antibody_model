# Fit model relating vaccine efficacy to antibody titres

# Most code can only be run with access to original trial data

# Parameters to be estimated: 
# EIR, r_clin, gamma_llin, v_max, alpha, beta
# Seasonality parameters are estimated in a first step by fitting to incidence in 
# the unvaccinated group

library(malariaEquilibrium)
library(cpp11)
library(tidyverse)
library(mets)
library(drjacoby)
library(umbrella)
library(ggplot2)

source("efficacy_model/R/antibody_model.R")
source("efficacy_model/R/vaccine_efficacy_model.R")
source("efficacy_model/R/eq.R")
source("efficacy_model/R/season.R")
source("efficacy_model/R/hazard.R")
source("efficacy_model/R/efficacy_fitting.R")

### Read in fixed parameters and data #####################################################

### Fitted seasonality parameters

# 3-year seasonality profile over whole trial
seas_profile <- readRDS("efficacy_model/output/seas_profile.RDS")
# Average parameters for each 1 year cycle
seas_profile2 <- readRDS("efficacy_model/output/seas_profile2.RDS")
# For forward projection, combine the 2

### Antibody parameters 
ab_parms <- read.csv("efficacy_model/data/antibody_global_fit.csv")   # geometric mean titer
ab_parms_mean <- filter(ab_parms, X == "mean")
peak <- c(ab_parms_mean$peak1, ab_parms_mean$peak2, ab_parms_mean$peak3)
t_boost <- c(364, 729)

# Processed trial data for fitting: not publicly available

# Generates the following objects:
# data - complete dataset
# events - subset of dataset of all participants experiencing a clinical case
# time_to_event1_itn - subset of dataset for participants in vaccine group with bednet use
# time_to_event1_no_itn - subset of dataset for participants in vaccine group without bednet use
# time_to_event2_itn - subset of dataset for participants in control group with bednet use
# time_to_event2_no_itn - subset of dataset for participants in control group without bednet use

### Extract age groups and bednet coverage from dataset #####################################

# Find average age at enrollment for each age group
age_to_model <- group_by(data, age_group) %>%
  summarise(median = median(age_at_enrollment),
            n_id = n_distinct(id)) %>%
  arrange(median)

itn_use <- data %>%
  select(id, subgroup, age_group, bednet) %>%
  distinct() %>%
  group_by(subgroup, age_group, bednet) %>%
  tally() %>%
  mutate(bednet = str_replace(bednet, '0', "no_bednet"),
         bednet = str_replace(bednet, '1', "bednet")) %>%
  pivot_wider(id_cols = c(subgroup, age_group), names_from = bednet, values_from = n) %>%
  mutate(no_bednet = ifelse(is.na(no_bednet), 0, no_bednet),
         itn_use = bednet/(bednet+no_bednet)) %>%
  select(subgroup, age_group, itn_use) %>%
  arrange(subgroup, age_group)

# Averaged over all age groups
itn_use_summary <- data %>%
  select(id, subgroup, bednet) %>%
  distinct() %>%
  group_by(subgroup, bednet) %>%
  tally() %>%
  mutate(bednet = str_replace(bednet, '0', "no_bednet"),
         bednet = str_replace(bednet, '1', "bednet")) %>%
  pivot_wider(id_cols = c(subgroup), names_from = bednet, values_from = n) %>%
  mutate(no_bednet = ifelse(is.na(no_bednet), 0, no_bednet),
         itn_use = bednet/(bednet+no_bednet)) %>%
  select(subgroup, itn_use) %>%
  arrange(subgroup)

### Estimate seasonality parameters from incidence in unvaccinated cohort ------

# Cases are grouped into monthly periods, starting from the day when the
# first participant enrolled (28 days after 3rd dose)
# All cases occurring within a month are summarised at the mid-point of that
# month since first enrollment

control_cases <- filter(events, group == 3) %>%
  # Recalculate calendar date for event 
  mutate(calendar_date = as.Date(eventtime, origin = start_date)) %>%
  group_by(calendar_date, start_date) %>%
  summarise(cases=n()) %>%
  ungroup() %>%
  mutate(time_since_first_enrolment = as.numeric(calendar_date - min(start_date)),
         month = ceiling(time_since_first_enrolment/30)-1,
         month_midpoint = month+0.5,
         day = round(month_midpoint * 365/12,0)) %>%
  group_by(month, month_midpoint, day) %>%
  summarise(cases=sum(cases)) %>%
  ungroup()

# Fit seasonality in each year separately using umbrella package
control_cases3 <- control_cases %>%
  mutate(year = floor(day/365),
         day = day %% 365) %>%
  select(day, year, cases)

# Fit year 1
# Fit Fourier series
seas_parms1 <- fit_fourier(rainfall = control_cases3$cases[control_cases3$year == 0], 
                           t = control_cases3$day[control_cases3$year == 0],
                           floor = 0)
seas_parms2 <- fit_fourier(rainfall = control_cases3$cases[control_cases3$year == 1], 
                           t = control_cases3$day[control_cases3$year == 1],
                           floor = 0)
seas_parms3 <- fit_fourier(rainfall = control_cases3$cases[control_cases3$year == 2], 
                           t = control_cases3$day[control_cases3$year == 2],
                           floor = 0)

# Check fit
predict_seas1 <- fourier_predict(coef = seas_parms1$coefficients, t = 1:365, floor = 0.01)
predict_seas2 <- fourier_predict(coef = seas_parms2$coefficients, t = 1:365, floor = 0.01)
predict_seas3 <- fourier_predict(coef = seas_parms3$coefficients, t = 1:365, floor = 0.01)

seas_profile_not_scaled <- c(predict_seas1$profile,predict_seas2$profile,predict_seas3$profile)
seas_profile_scaled <- length(seas_profile_not_scaled)/sum(seas_profile_not_scaled) * seas_profile_not_scaled
#saveRDS(seas_profile_scaled, "./efficacy_model/output/seas_profile.RDS")

# Plot scaled seasonality profile
ggplot() +
  geom_line(aes(x = 1:(3*365),
                y = seas_profile)) +
  #geom_point(aes(x=control_cases$day, 
  #               y =control_cases$cases),
  #           colour = "red") +
  labs(x = "Days", y = "Seasonality profile") +
  theme_classic()

### Run MCMC ###################################################################
# Can only be run with original trial data

df_params <- data.frame(
  name = c("EIR", "gamma_llin", "r_clin", "v_max", "alpha", "beta",
           "rho1", "rho2", "d_s", "d_l"),
  min = c(0, 0, 0, 0, 0, 0,
          -Inf, -Inf, 0, 0),
  max = c(Inf, 10, Inf, 1, Inf, 11159,
          Inf, Inf, Inf, Inf))

mcmc <- drjacoby::run_mcmc(
  data = list(time_to_event1_itn = time_to_event1_itn,
              time_to_event1_no_itn = time_to_event1_no_itn,
              time_to_event2_itn = time_to_event2_itn,
              time_to_event2_no_itn = time_to_event2_no_itn),
  df_params = df_params,
  loglike = efficacy_loglike,
  logprior = efficacy_logprior,
  burnin = 2000,
  samples = 25000, 
  chains = 4,
  misc = list(seasonality = seas_profile,
              age_to_model = age_to_model,
              peak = peak, t_boost = t_boost))

# Check traces
plot_par(mcmc, show = "EIR", phase = "burnin")
plot_par(mcmc, show = "EIR", phase = "sampling")
plot_par(mcmc, show = "gamma_llin", phase = "burnin")
plot_par(mcmc, show = "gamma_llin", phase = "sampling")
plot_par(mcmc, show = "r_clin", phase = "burnin")
plot_par(mcmc, show = "r_clin", phase = "sampling")
plot_par(mcmc, show = "v_max", phase = "burnin")
plot_par(mcmc, show = "v_max", phase = "sampling")
plot_par(mcmc, show = "alpha", phase = "burnin")
plot_par(mcmc, show = "alpha", phase = "sampling")
plot_par(mcmc, show = "beta", phase = "burnin")
plot_par(mcmc, show = "beta", phase = "sampling")

plot_cor(mcmc, "EIR", "r_clin")

# Convergence diagnostic
mcmc$diagnostics$rhat
# Effective sample size
mcmc$diagnostics$ess

### Posterior estimates ##########################################################

posterior_estimates <- mcmc$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
apply(posterior_estimates,2,round, 2)

# Save medians
eir_med <- posterior_estimates[1,1]
gamma_llin_med <- posterior_estimates[2,1]
r_clin_med <- posterior_estimates[3,1]
v_max_med <- posterior_estimates[4,1]
alpha_med <- posterior_estimates[5,1]
beta_med <- posterior_estimates[6,1]
rho_med <- c(1 / (1 + exp(-posterior_estimates[7,1])), 
             1 / (1 + exp(-posterior_estimates[8,1])),
             1 / (1 + exp(-posterior_estimates[8,1])))
d_s_med <- posterior_estimates[9,1]
d_l_med <- posterior_estimates[10,1]

param_summary <- data.frame(parameter = rownames(posterior_estimates),
                            median = posterior_estimates[,1],
                            cri_lo = posterior_estimates[,2],
                            cri_hi = posterior_estimates[,3])
param_summary[param_summary$parameter=="rho1",-1] <- 1 / 
  (1 + exp(-param_summary[param_summary$parameter=="rho1",-1]))
param_summary[param_summary$parameter=="rho2",-1] <- 1 / 
  (1 + exp(-param_summary[param_summary$parameter=="rho2",-1]))
# write.csv(param_summary, "./efficacy_model/output/efficacy_parameters_summary.csv",
#           row.names = FALSE)

# Draw 50 samples
set.seed(1000)
param_draws <- sample_chains(mcmc, 50)
# Transform rhos
param_draws$rho1 <- 1 / (1 + exp(-param_draws$rho1))
param_draws$rho2 <- 1 / (1 + exp(-param_draws$rho2))

# Save for malariasimulation runs
#write.csv(param_draws[,c("v_max", "alpha", "beta")], "./efficacy_model/output/efficacy_parameters.csv",
#          row.names = FALSE)
#write.csv(param_draws, "./efficacy_model/output/efficacy_parameters_full.csv",
#          row.names = FALSE)

### Plot fit (Figure 1B) #######################################################

# Load relevant parameters (can be run without full data access)
param_draws <- read.csv("./efficacy_model/output/efficacy_parameters_full.csv")
# Create new seasonality profile assuming average following the 3 year period
seas_profile_combi <- c(seas_profile, rep(seas_profile2, 9))
age_at_enrollment <- median(data$age_at_enrollment)   
param_summary <- read.csv("./efficacy_model/output/efficacy_parameters_summary.csv")
rownames(param_summary) <- param_summary$parameter
eir_med <- param_summary$median[param_summary$parameter == "EIR"]
gamma_llin_med <- param_summary$median[param_summary$parameter == "gamma_llin"]
r_clin_med <- param_summary$median[param_summary$parameter == "r_clin"]
v_max_med <- param_summary$median[param_summary$parameter == "v_max"]
alpha_med <- param_summary$median[param_summary$parameter == "alpha"]
beta_med <- param_summary$median[param_summary$parameter == "beta"]
rho1_med <- param_summary$median[param_summary$parameter == "rho1"]
rho2_med <- param_summary$median[param_summary$parameter == "rho2"]
d_s_med <- param_summary$median[param_summary$parameter == "d_s"]
d_l_med <- param_summary$median[param_summary$parameter == "d_l"]
rho_med <- c(rho1_med, rho2_med, rho2_med) 

# Fit with median posterior parameter values
fitted_efficacy_med <- simulate_trial_hazards(eir = eir_med,
                                              age_at_enrollment = age_at_enrollment,
                                              gamma_llin = gamma_llin_med,
                                              r_clin = r_clin_med,
                                              peak = peak,
                                              rho = rho_med,
                                              t_boost = t_boost,
                                              d_s = d_s_med, d_l = d_l_med,
                                              v_max = v_max_med,
                                              alpha = alpha_med, beta = beta_med,
                                              age =  1:(365 * 12), s2 = 1.16,
                                              season = seas_profile_combi)

# Take weighted average of clinical hazard between ITN users and non-users
c2b_med <- fitted_efficacy_med$c2b_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]) +
  fitted_efficacy_med$c2b_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]
c2a_med <- fitted_efficacy_med$c2a_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]) +
  fitted_efficacy_med$c2a_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]
c3_med <- fitted_efficacy_med$c3_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) +
  fitted_efficacy_med$c3_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]

# Fit with lower percentile
rho_lower <- c(param_summary["rho1",4], param_summary["rho2",4], param_summary["rho2",4])
fitted_efficacy_lower <- simulate_trial_hazards(eir = param_summary["EIR",3],
                                                age_at_enrollment = age_at_enrollment,
                                                gamma_llin = param_summary["gamma_llin",3],
                                                r_clin = param_summary["r_clin",3],
                                                peak = peak,
                                                rho = rho_lower,
                                                t_boost = t_boost,
                                                d_s = param_summary["d_s",3],
                                                d_l = param_summary["d_l",3],
                                                v_max = param_summary["v_max",3],
                                                alpha = param_summary["alpha",4],
                                                beta = param_summary["beta",4],
                                                age =  1:(365 * 12), s2 = 1.16,
                                                season = seas_profile_combi)

# Take weighted average of clinical hazard between ITN users and non-users
c2b_lower <- fitted_efficacy_lower$c2b_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]) +
  fitted_efficacy_lower$c2b_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]
c2a_lower <- fitted_efficacy_lower$c2a_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]) +
  fitted_efficacy_lower$c2a_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]
c3_lower <- fitted_efficacy_lower$c3_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) +
  fitted_efficacy_lower$c3_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]

# Fit with upper percentile
rho_upper <- c(param_summary["rho1",3], param_summary["rho2",3],param_summary["rho2",3])
fitted_efficacy_upper <- simulate_trial_hazards(eir = param_summary["EIR",4],
                                                age_at_enrollment = age_at_enrollment,
                                                gamma_llin = param_summary["gamma_llin",4],
                                                r_clin = param_summary["r_clin",4],
                                                peak = peak,
                                                rho = rho_upper,
                                                t_boost = t_boost,
                                                d_s = param_summary["d_s",4],
                                                d_l = param_summary["d_l",4],
                                                v_max = param_summary["v_max",4],
                                                alpha = param_summary["alpha",3],
                                                beta = param_summary["beta",3],
                                                age =  1:(365 * 12), s2 = 1.16,
                                                season = seas_profile_combi)

# Take weighted average of clinical hazard between ITN users and non-users
c2b_upper <- fitted_efficacy_upper$c2b_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]) +
  fitted_efficacy_upper$c2b_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2b"]
c2a_upper <- fitted_efficacy_upper$c2a_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]) +
  fitted_efficacy_upper$c2a_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="2a"]
c3_upper <- fitted_efficacy_upper$c3_no_itn*(1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) +
  fitted_efficacy_upper$c3_itn*itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]

efficacy_est1 <- data.frame(time = 1:length(c2b_med),
                            ve_med = (1-c2b_med/c3_med),
                            ve_lower = (1-c2b_lower/c3_lower),
                            ve_upper = (1-c2b_upper/c3_upper)) %>%
  filter(time<(5*365))
efficacy_est2 <- data.frame(time = 1:length(c2a_med),
                            ve_med = (1-c2a_med/c3_med),
                            ve_lower = (1-c2a_lower/c3_lower),
                            ve_upper = (1-c2a_upper/c3_upper)) %>%
  filter(time<(5*365))

fitp1 <- ggplot() +
  geom_ribbon(data = efficacy_est1,
              aes(x= time/365,ymin = ve_lower * 100, ymax = ve_upper * 100),
              fill = "red", alpha = 0.2) +
  geom_ribbon(data = efficacy_est2,
              aes(x= time/365,ymin = ve_lower * 100, ymax = ve_upper * 100),
              fill = "steelblue3", alpha = 0.2) +
  geom_line(data = efficacy_est2,
            aes(x = time/365,
                y= ve_med * 100, linetype = "Model"),
            col = "steelblue3", linewidth=1) +
  geom_line(data = efficacy_est1,
            aes(x = time/365,
                y= ve_med * 100, linetype = "Model"),
            col = "red", linewidth=1) +
  geom_pointrange(data= filter(vaccine_efficacy_df, comparison != "group2a (high-dose booster2) vs group3"),
                  aes(x=midpoint/365, y = efficacy*100,
                      ymin = efficacy_lo*100,
                      ymax = efficacy_hi*100, col = "No second booster", shape = "Data"),
                  size = 1, linewidth=0.8) +
  geom_pointrange(data = filter(vaccine_efficacy_df, comparison == "group2a (high-dose booster2) vs group3"),
                  aes(x=(midpoint+25)/365, y = efficacy*100,
                      ymin = efficacy_lo*100,
                      ymax = efficacy_hi*100, col = "With second booster", shape = "Data"),
                  size = 1, linewidth=0.8) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks = c(-25, 0, 25, 50, 75, 100)) +
  scale_colour_manual(values = c("No second booster" = "red3",
                                 "With second booster" = "steelblue4")) +
  geom_vline(xintercept = 1064/365, linetype = "dashed", col = "grey") +
  coord_cartesian(xlim = c(0,5), ylim = c(-25,100), expand = FALSE) +  # maybe remove
  ylab("Clinical efficacy (%)") + xlab("Time after third dose (years)") +
  labs(col = "", linetype = "", shape = "", tag = "B") +
  theme_classic() +
  guides(fill = guide_legend(byrow = TRUE), linetype = "none", shape = "none") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        legend.text = element_text(size = 10),
        legend.position = "none",
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, 'cm'))

### Plot other outputs for Figure 1 ###########################################

# Antibody titre with 2 boosters
ab_med <- antibody_titre(1:(365*15),
                         peak = peak,
                         rho = rho_med,
                         t_boost = t_boost,
                         d_s = d_s_med,
                         d_l = d_l_med)
# Antibody titre with a single booster
ab_med2 <- antibody_titre(1:(365*15),
                          peak = peak[1:2],
                          rho = rho_med[1:2],
                          t_boost = t_boost[1],
                          d_s = d_s_med,
                          d_l = d_l_med)

# Efficacy against infection over time
dose_response_med <- data.frame(count = 1:length(ab_med),
                                ab_2boost = ab_med,
                                efficacy_inf_2boost = 
                                  vaccine_efficacy(ab_med, v_max = v_max_med, 
                                                   alpha = alpha_med, beta = beta_med),
                                ab_1boost = ab_med2,
                                efficacy_inf_1boost = 
                                  vaccine_efficacy(ab_med2, v_max = v_max_med, 
                                                   alpha = alpha_med, beta = beta_med))

# Dose-response relationship between vaccine efficacy against infection and antibody titre (Figure 1C)
fitp3 <- ggplot() +
  geom_line(data=dose_response_med, 
            aes(x = ab_2boost, y = efficacy_inf_2boost*100), size = 1.3) +
  ylab("Efficacy against infection (%)\n") +
  xlab("Anti-CSP antibody titre (EU/mL)") +
  scale_x_continuous(trans="log10", breaks = c(10,100,10000)) +
  coord_trans(x="log10", ylim = c(0, 100),
              xlim=c(min(dose_response_med$ab_2boost),55656), 
              expand = FALSE,
              clip = "off") +
  labs(tag = "C") +
  theme_classic()  +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14)) 

# Efficacy against infection over time (Figure 1D)
fitp2 <- ggplot(dose_response_med) +
  geom_line(aes(x = ((1:length(ab_2boost))+28)/365, 
                y= efficacy_inf_2boost*100),  size = 1.3, col = "steelblue3") +
  geom_line(aes(x = ((1:length(ab_1boost))+28)/365, 
                y= efficacy_inf_1boost*100), size = 1.3, col = "red") +
  coord_cartesian(xlim = c(0,5), ylim = c(0, 100), expand = FALSE) +
  ylab("Efficacy against infection (%)") + xlab("Time after third dose (years)") +
  labs(tag = "D") +
  theme_classic() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14))

# Antibody fit (Figure 1A)
# Outputs from ab_model scripts

# Adjust to start at day 28 on the plot (after 3rd dose)
global_dat_pd_1B <- readRDS("efficacy_model/data/global_dat_pd_1B.RDS")  %>%
  mutate(t=t-56)
global_pred_1B <- readRDS("efficacy_model/data/global_pred_1B.RDS") %>%
  mutate(t=t-56) %>%
  filter(t <= (5*365)) 
global_dat_pd_2B <- readRDS("efficacy_model/data/global_dat_pd_2B.RDS") %>%
  mutate(t=t-56)
global_pred_2B <- readRDS("efficacy_model/data/global_pred_2B.RDS") %>%
  mutate(t=t-56) %>%
  filter(t <= (5*365))

global_pred_2B2 <- filter(global_pred_2B, group == "Post boost 2")
global_dat_pd_2B2 <- filter(global_dat_pd_2B, scenario == "Post boost 2")

fitp4 <- ggplot() + 
  geom_ribbon(data = global_pred_1B, aes(x = t, ymin = abll, ymax = abuu), 
              alpha = 0.1, fill = "red") +  # darkorchid3
  geom_line(data = global_pred_1B, aes(x = t, y = ab), lwd = 1, colour = "red") +   # darkorchid4
  geom_ribbon(data = global_pred_2B2, aes(x = t, ymin = abll, ymax = abuu), 
              alpha = 0.1, fill = "steelblue3") +   # springgreen3
  geom_line(data = global_pred_2B2, aes(x = t, y = ab), lwd = 1, col = "steelblue3") +
  geom_linerange(data = global_dat_pd_1B, aes(x = t, ymin = ymin, ymax = ymax,
                                              col = "No second booster"),
                 linewidth=0.8) +
  geom_point(data = global_dat_pd_1B, aes(x = t, y=y, 
                                          col = "No second booster"), size = 4) +
  geom_linerange(data = global_dat_pd_2B2, aes(x = t, ymin = ymin, ymax = ymax,
                                               col = "With second booster"),
                 linewidth=0.8) +
  geom_point(data = global_dat_pd_2B2, aes(x = t, y=y,
                                           col = "With second booster"),
             size = 4, alpha = 0.8) +
  scale_fill_discrete(guide = "none") +
  scale_colour_manual(values = c("No second booster" = "red3",
                                 "With second booster" = "steelblue4")) +
  scale_y_log10(limits = c(0.1, NA), breaks = c(1,100,10000)) + 
  scale_x_continuous(breaks = seq(0, 365 * 5, 365), labels = 0:5, limits=c(0,5*365),
                     name = "Time after third dose (years)") +
  ylab("Anti-CSP antibody titre (EU/mL)") +
  labs(tag = "A") +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        legend.text = element_text(size = 10),
        legend.position=c(.25,.15),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.001, 'cm'))

# Combined plot

fig1 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(fitp4, fitp3, nrow = 2), 
                                   gridExtra::arrangeGrob(fitp1, fitp2, nrow = 2), 
                                   ncol=2)
#ggsave("efficacy_model/output/fig1.png", fig1, width = 27, height = 18, units="cm")
