# Model validation to Phase 3 data

# Most code can only be run with access to original trial data

# Steps: 
# 1. Fit seasonality to control group data
# 2. Fit EIR to control group data
# 3. Project vaccine impact using fitted model
# 4. Plot comparison

library(malariaEquilibrium)
library(cpp11)
library(tidyverse)
library(mets)
library(drjacoby)
library(umbrella)

source("efficacy_model/R/antibody_model.R")
source("efficacy_model/R/vaccine_efficacy_model.R")
source("efficacy_model/R/eq.R")
source("efficacy_model/R/season.R")
source("efficacy_model/R/hazard.R")
source("efficacy_model/R/efficacy_fitting.R")

### Read in fixed parameters and data ##########################################

param_draws <- read.csv("./efficacy_model/output/efficacy_parameters_full.csv")
fix_gamma_llin <- median(param_draws$gamma_llin)
fix_r_clin <- median(param_draws$r_clin)
ab_parms <- read.csv("efficacy_model/data/antibody_global_fit.csv")
ab_parms_mean <- filter(ab_parms, X == "mean")
peak <- c(ab_parms_mean$peak1, ab_parms_mean$peak2)

# Data from control arm in Phase 3 trial sites: not publicly available

# Generates the following objects:
# data - complete dataset
# events - subset of dataset of all participants experiencing a clinical case
# time_to_event2_itn - subset of dataset for participants in control group with bednet use
# time_to_event2_no_itn - subset of dataset for participants in control group without bednet use
# time_at_booster - median day of receiving booster following enrollment in each site

# Extract age groups and bednet coverage 

# Find average age at enrollment for each age group (in days)
age_to_model <- group_by(data, site_id, age_group) %>%
  summarise(median = median(age_at_enrollment*30),
            n_id = n_distinct(id)) %>%
  arrange(median)

itn_use <- data %>%
  select(id, site_id, subgroup, age_group, bednet) %>%
  distinct() %>%
  group_by(subgroup, site_id, age_group, bednet) %>%
  tally() %>%
  mutate(bednet = str_replace(bednet, '0', "no_bednet"),
         bednet = str_replace(bednet, '1', "bednet")) %>%
  pivot_wider(id_cols = c(subgroup,site_id, age_group), names_from = bednet, values_from = n) %>%
  mutate(no_bednet = ifelse(is.na(no_bednet), 0, no_bednet),
         itn_use = bednet/(bednet+no_bednet)) %>%
  select(subgroup, site_id, age_group, itn_use) %>%
  arrange(subgroup, age_group)

# Averaged over all age groups
itn_use_summary <- data %>%
  select(id, site_id, subgroup, bednet) %>%
  distinct() %>%
  group_by(subgroup, site_id, bednet) %>%
  tally() %>%
  mutate(bednet = str_replace(bednet, '0', "no_bednet"),
         bednet = str_replace(bednet, '1', "bednet")) %>%
  pivot_wider(id_cols = c(subgroup, site_id), names_from = bednet, values_from = n) %>%
  mutate(no_bednet = ifelse(is.na(no_bednet), 0, no_bednet),
         itn_use = bednet/(bednet+no_bednet)) %>%
  select(subgroup, site_id, itn_use) %>%
  arrange(subgroup)

### 1. Estimate seasonality parameters from incidence in unvaccinated cohort ------

# Plot cases in control group 
control_cases <- filter(data, status == 1) %>%
  # Recalculate calendar date for event 
  mutate(calendar_date = as.Date(eventtime, origin = start_date)) %>%
  group_by(site_id, eventtime, calendar_date, start_date) %>%
  summarise(cases=n())

sites <- unique(control_cases$site_id)  
control_cases_by_site <- list()
floors <- c(0.0001, 0, 0.0001, 0.0001, 0)

# Seasonality is based on calendar time so the starting point is the first date of enrollment
# in each site
for(i in 1:length(unique(sites))) {
  
  control_cases_by_site[[i]] <- filter(control_cases, site_id == sites[i]) %>%
    group_by(site_id) %>%
    mutate(time_since_first_enrolment = as.numeric(calendar_date - min(start_date)),
           month = ceiling(time_since_first_enrolment/30)-1,
           month_midpoint = month+0.5,
           day = round(month_midpoint * 365/12,0)) %>%
    group_by(site_id, month, month_midpoint, day) %>%
    summarise(cases = sum(cases)) %>%
    ungroup()
  
}
control_cases_by_site <- do.call("rbind", control_cases_by_site)

control_cases_by_site <- left_join(control_cases_by_site, time_at_booster,
                                   by = "site_id")

control_cases_by_site$site_id[control_cases_by_site$site_id == "Nanoro & Sigle"] <- "Nanoro"
control_cases_by_site$site_id <- factor(control_cases_by_site$site_id, 
                                        levels = c("Bagamoyo", "Kilifi", "Dande",
                                                   "Bougouni", "Nanoro"))
sites[sites=="Nanoro & Sigle"]  <- "Nanoro"

control_cases_by_site <- control_cases_by_site %>%
  group_by(site_id) %>%
  mutate(scaled = (cases-min(cases))/(max(cases)-min(cases)))

# Use umbrella package for fitting

# Fit seasonality in each year separately
control_cases3 <- control_cases_by_site %>%
  group_by(site_id) %>%
  mutate(year2 = floor(day/365),
         day = (day %% 365)+1) %>%
  select(day, year2, cases)

# For second year include last point from year 1 in fitting and drop this after
control_cases3_year2 <- filter(control_cases3, year2 == 1 | day == 351) %>%
  mutate(day2 = ifelse(day == 351, 15, day+30))

# Fit Fourier series with floor = 0.01
seas_profile <- data.frame(site_id = NA, day = NA, seas_profile = NA,
                           seas_profile_scaled = NA)

seas_parms_list <- list()

for (i in 1:length(sites)) {
  # year 1
  seas_parms1 <- fit_fourier(rainfall = control_cases3$cases[
    control_cases3$year2 == 0 & 
      control_cases3$site_id == sites[i]], 
    t = control_cases3$day[
      control_cases3$year2 == 0 &
        control_cases3$site_id == sites[i]],
    floor = floors[i])

  seas_parms2 <- fit_fourier(rainfall = control_cases3_year2$cases[
    control_cases3_year2$site_id == sites[i]],
    t = control_cases3_year2$day2[
      control_cases3_year2$site_id == sites[i]],
    floor = floors[i])
  
  seas_parms_list[[i]] <- data.frame(rbind(seas_parms1$coefficients, seas_parms2$coefficients)) %>%
    mutate(site = sites[i],
           year=c(1,2))
  
  # Check fit
  predict_seas1 <- fourier_predict(coef = seas_parms1$coefficients, t = 1:365,
                                   floor = 0.0001)
  predict_seas2 <- fourier_predict(coef = seas_parms2$coefficients, t = 1:365,
                                   floor = 0.0001)
  
  predict_seas2_profile <-  predict_seas2$profile[-c(1:30)]
  
  seas_profile1 <- c(predict_seas1$profile,predict_seas2_profile)[
    1:max(365*control_cases3$year2[control_cases3$site_id == sites[i]]+
            control_cases3$day[control_cases3$site_id == sites[i]])]
  seas_profile_scaled <- length(seas_profile1)/sum(seas_profile1) * seas_profile1
  seas_parms_list[[i]]$scaler <- length(seas_profile1)/sum(seas_profile1)
  
  print(sum(seas_profile_scaled))
  print(mean(seas_profile_scaled))
  
  seas_profile <- rbind(seas_profile,
                        data.frame(site_id = sites[i],
                                   day = 1:max(365*control_cases3$year2[control_cases3$site_id == sites[i]]+
                                                 control_cases3$day[control_cases3$site_id == sites[i]]),
                                   seas_profile = seas_profile1,
                                   seas_profile_scaled = seas_profile_scaled))
  
}
seas_profile <- drop_na(seas_profile)
seas_parms_list <- do.call("rbind", seas_parms_list)
#saveRDS(seas_profile, "./efficacy_model/output/validation_seas_profile.RDS")

### Split data into site files ################################################

# For t_boost, assuming no effect until 28 days after median booster timing 
# (same assumption as in fitting to Phase 2 data)

seas_profile <- readRDS("./efficacy_model/output/validation_seas_profile.RDS")

# Bagamoyo
time_to_event2_itn_bagamoyo <- filter(time_to_event2_itn, site_id == "Bagamoyo") %>%
  select(-site_id, -site_type)
time_to_event2_no_itn_bagamoyo <- filter(time_to_event2_no_itn, site_id == "Bagamoyo") %>%
  select(-site_id, -site_type)
seas_profile_bagamoyo <- seas_profile$seas_profile_scaled[seas_profile$site_id == "Bagamoyo"]
itn_use_summary_bagamoyo <- filter(itn_use_summary, site_id == "Bagamoyo")
age_to_model_bagamoyo <- filter(age_to_model, site_id == "Bagamoyo")
t_boost_bagamoyo <- time_at_booster$median_timing[time_at_booster$site_id == "Bagamoyo"]+28   
peak_bagamoyo <- 1*peak

# Kilifi
time_to_event2_itn_kilifi <- filter(time_to_event2_itn, site_id == "Kilifi") %>%
  select(-site_id, -site_type)
time_to_event2_no_itn_kilifi <- filter(time_to_event2_no_itn, site_id == "Kilifi") %>%
  select(-site_id, -site_type)
seas_profile_kilifi <- seas_profile$seas_profile_scaled[seas_profile$site_id == "Kilifi"]
itn_use_summary_kilifi <- filter(itn_use_summary, site_id == "Kilifi")
age_to_model_kilifi <- filter(age_to_model, site_id == "Kilifi")
t_boost_kilifi <- time_at_booster$median_timing[time_at_booster$site_id == "Kilifi"]+28   
peak_kilifi <- 1*peak

# Dande
time_to_event2_itn_dande <- filter(time_to_event2_itn, site_id == "Dande") %>%
  select(-site_id, -site_type)
time_to_event2_no_itn_dande <- filter(time_to_event2_no_itn, site_id == "Dande") %>%
  select(-site_id, -site_type)
seas_profile_dande <- seas_profile$seas_profile_scaled[seas_profile$site_id == "Dande"]
itn_use_summary_dande <- filter(itn_use_summary, site_id == "Dande")
age_to_model_dande <- filter(age_to_model, site_id == "Dande")
t_boost_dande <- time_at_booster$median_timing[time_at_booster$site_id == "Dande"]+28   
peak_dande <- 1*peak

# Bougouni
time_to_event2_itn_bougouni <- filter(time_to_event2_itn, site_id == "Bougouni") %>%
  select(-site_id, -site_type)
time_to_event2_no_itn_bougouni <- filter(time_to_event2_no_itn, site_id == "Bougouni") %>%
  select(-site_id, -site_type)
seas_profile_bougouni <- seas_profile$seas_profile_scaled[seas_profile$site_id == "Bougouni"]
itn_use_summary_bougouni <- filter(itn_use_summary, site_id == "Bougouni")
age_to_model_bougouni <- filter(age_to_model, site_id == "Bougouni")
t_boost_bougouni <- time_at_booster$median_timing[time_at_booster$site_id == "Bougouni"]+28   
peak_bougouni <- 1*peak

# Nanoro
time_to_event2_itn_nanoro <- filter(time_to_event2_itn, site_id == "Nanoro & Sigle") %>%
  select(-site_id, -site_type)
time_to_event2_no_itn_nanoro <- filter(time_to_event2_no_itn, site_id == "Nanoro & Sigle") %>%
  select(-site_id, -site_type)
seas_profile_nanoro <- seas_profile$seas_profile_scaled[seas_profile$site_id == "Nanoro"]
itn_use_summary_nanoro <- filter(itn_use_summary, site_id == "Nanoro & Sigle")
age_to_model_nanoro <- filter(age_to_model, site_id == "Nanoro & Sigle")
t_boost_nanoro <- time_at_booster$median_timing[time_at_booster$site_id == "Nanoro & Sigle"]+28   
peak_nanoro <- 1*peak

### 2. MCMC to estimate EIR in each site #######################################

# Bagamoyo
df_params <- data.frame(
  name = c("EIR"),
  min = c(0),
  max = c(100/365))

mcmc_bagamoyo <- drjacoby::run_mcmc(
  data = list(time_to_event2_itn = time_to_event2_itn_bagamoyo,
              time_to_event2_no_itn = time_to_event2_no_itn_bagamoyo),
  df_params = df_params,
  loglike = validation_loglike,
  logprior = validation_logprior,
  burnin = 1000,
  samples = 1000, 
  chains = 4,
  misc = list(seasonality = seas_profile_bagamoyo,
              age_to_model = age_to_model_bagamoyo,
              gamma_llin = fix_gamma_llin, r_clin = fix_r_clin))

# Kilifi
df_params <- data.frame(
  name = c("EIR"),
  min = c(0),
  max = c(100/365))

mcmc_kilifi <- drjacoby::run_mcmc(
  data = list(time_to_event2_itn = time_to_event2_itn_kilifi,
              time_to_event2_no_itn = time_to_event2_no_itn_kilifi),
  df_params = df_params,
  loglike = validation_loglike,
  logprior = validation_logprior,
  burnin = 1000,
  samples = 1000, 
  chains = 4,
  misc = list(seasonality = seas_profile_kilifi,
              age_to_model = age_to_model_kilifi,
              gamma_llin = fix_gamma_llin, r_clin = fix_r_clin))

# Dande
df_params <- data.frame(
  name = c("EIR"),
  min = c(0),
  max = c(100/365))

mcmc_dande <- drjacoby::run_mcmc(
  data = list(time_to_event2_itn = time_to_event2_itn_dande,
              time_to_event2_no_itn = time_to_event2_no_itn_dande),
  df_params = df_params,
  loglike = validation_loglike,
  logprior = validation_logprior,
  burnin = 1000,
  samples = 1000, 
  chains = 4,
  misc = list(seasonality = seas_profile_dande,
              age_to_model = age_to_model_dande,
              gamma_llin = fix_gamma_llin, r_clin = fix_r_clin))

# Bougouni
df_params <- data.frame(
  name = c("EIR"),
  min = c(0),
  max = c(100/365))

mcmc_bougouni <- drjacoby::run_mcmc(
  data = list(time_to_event2_itn = time_to_event2_itn_bougouni,
              time_to_event2_no_itn = time_to_event2_no_itn_bougouni),
  df_params = df_params,
  loglike = validation_loglike,
  logprior = validation_logprior,
  burnin = 1000,
  samples = 1000, 
  chains = 4,
  misc = list(seasonality = seas_profile_bougouni,
              age_to_model = age_to_model_bougouni,
              gamma_llin = fix_gamma_llin, r_clin = fix_r_clin))

# Nanoro
df_params <- data.frame(
  name = c("EIR"),
  min = c(0),
  max = c(100/365))

mcmc_nanoro <- drjacoby::run_mcmc(
  data = list(time_to_event2_itn = time_to_event2_itn_nanoro,
              time_to_event2_no_itn = time_to_event2_no_itn_nanoro),
  df_params = df_params,
  loglike = validation_loglike,
  logprior = validation_logprior,
  burnin = 1000,
  samples = 1000, 
  chains = 4,
  misc = list(seasonality = seas_profile_nanoro,
              age_to_model = age_to_model_nanoro,
              gamma_llin = fix_gamma_llin, r_clin = fix_r_clin))

### MCMC diagnostics and extract EIR #############################################

# Bagamoyo
plot_par(mcmc_bagamoyo, show = "EIR", phase = "burnin")
plot_par(mcmc_bagamoyo, show = "EIR", phase = "sampling")
mcmc_bagamoyo$diagnostics$rhat
mcmc_bagamoyo$diagnostics$ess

posterior_estimates_bagamoyo <- mcmc_bagamoyo$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
eir_med_bagamoyo <- posterior_estimates_bagamoyo[1,1] 
set.seed(123)
eir_draws_bagamoyo <- sample_chains(mcmc_bagamoyo, 50)

# Kilifi
plot_par(mcmc_kilifi, show = "EIR", phase = "burnin")
plot_par(mcmc_kilifi, show = "EIR", phase = "sampling")
mcmc_kilifi$diagnostics$rhat
mcmc_kilifi$diagnostics$ess

posterior_estimates_kilifi <- mcmc_kilifi$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
eir_med_kilifi <- posterior_estimates_kilifi[1,1] 
set.seed(123)
eir_draws_kilifi <- sample_chains(mcmc_kilifi, 50)

# Dande

plot_par(mcmc_dande, show = "EIR", phase = "burnin")
plot_par(mcmc_dande, show = "EIR", phase = "sampling")
mcmc_dande$diagnostics$rhat
mcmc_dande$diagnostics$ess

posterior_estimates_dande <- mcmc_dande$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
eir_med_dande <- posterior_estimates_dande[1,1] 
set.seed(123)
eir_draws_dande <- sample_chains(mcmc_dande, 50)

# Bougouni
plot_par(mcmc_bougouni, show = "EIR", phase = "burnin")
plot_par(mcmc_bougouni, show = "EIR", phase = "sampling")
mcmc_bougouni$diagnostics$rhat
mcmc_bougouni$diagnostics$ess

posterior_estimates_bougouni <- mcmc_bougouni$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
eir_med_bougouni <- posterior_estimates_bougouni[1,1] 
set.seed(123)
eir_draws_bougouni <- sample_chains(mcmc_bougouni, 50)

# Nanoro
plot_par(mcmc_nanoro, show = "EIR", phase = "burnin")
plot_par(mcmc_nanoro, show = "EIR", phase = "sampling")
mcmc_nanoro$diagnostics$rhat
mcmc_nanoro$diagnostics$ess

posterior_estimates_nanoro <- mcmc_nanoro$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  summarise(across(everything(), quantile, c(0.5, 0.025, 0.975))) |>
  signif(3) |>
  t()
eir_med_nanoro <- posterior_estimates_nanoro[1,1] 
set.seed(123)
eir_draws_nanoro <- sample_chains(mcmc_nanoro, 50)

#saveRDS(eir_draws_bagamoyo, "efficacy_model/output/eir_draws_bagamoyo.RDS")
#saveRDS(eir_draws_kilifi, "efficacy_model/output/eir_draws_kilifi.RDS")
#saveRDS(eir_draws_dande, "efficacy_model/output/eir_draws_dande.RDS")
#saveRDS(eir_draws_bougouni, "efficacy_model/output/eir_draws_bougouni.RDS")
#saveRDS(eir_draws_nanoro, "efficacy_model/output/eir_draws_nanoro.RDS")

### Functions for running validation #########################################

run_draws <- function(eir_draws, age_at_enrollment,
                      site_peak, site_t_boost, seasonality_profile, 
                      itn_use_summary, vary_malariasimulation_parms = FALSE,
                      malariasimulation_parmset) {
  
  # Run draws
  rho <- list()
  fitted_efficacy_draws <- list()
  c2b <- list()
  c2a <- list()
  c3 <- list()
  
  for(i in 1:nrow(param_draws)) {
    rho[[i]] <- c(param_draws$rho1[i], param_draws$rho2[i])
    
    fitted_efficacy_draws[[i]] <- simulate_trial_hazards(eir =  eir_draws[i],
                                                         age_at_enrollment =  age_at_enrollment,
                                                         gamma_llin = fix_gamma_llin,
                                                         r_clin = fix_r_clin,
                                                         peak = site_peak, 
                                                         rho = rho[[i]], 
                                                         t_boost = site_t_boost, 
                                                         d_s =  param_draws$d_s[i], 
                                                         d_l = param_draws$d_l[i],
                                                         v_max = param_draws$v_max[i],
                                                         alpha = param_draws$alpha[i], 
                                                         beta = param_draws$beta[i],
                                                         age =  1:(12*365), s2 = 1.16,
                                                         season = seasonality_profile, 
                                                         vary_malariasimulation_parms = vary_malariasimulation_parms,
                                                         malariasimulation_parmset = malariasimulation_parmset)
    
    # Take weighted average of clinical hazard between ITN users and non-users
    
    c2b[[i]] <- fitted_efficacy_draws[[i]]$c2b_no_itn*
      (1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) + 
      fitted_efficacy_draws[[i]]$c2b_itn*
      itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]
    
    c2a[[i]] <- fitted_efficacy_draws[[i]]$c2a_no_itn*
      (1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) + 
      fitted_efficacy_draws[[i]]$c2a_itn*
      itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]
    
    c3[[i]] <- fitted_efficacy_draws[[i]]$c3_no_itn*
      (1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) + 
      fitted_efficacy_draws[[i]]$c3_itn*
      itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]
    
  }
  c2b_draws <- do.call("cbind", c2b) %>% 
    as_tibble() %>%
    mutate(time = 1:3650) %>%
    pivot_longer(cols = -time, names_to = "iteration", values_to = "c2b")
  c2a_draws <- do.call("cbind", c2a) %>% 
    as_tibble() %>%
    mutate(time = 1:3650) %>%
    pivot_longer(cols =-time, names_to = "iteration", values_to = "c2a")
  c3_draws <- do.call("cbind", c3) %>% 
    as_tibble() %>%
    mutate(time = 1:3650) %>%
    pivot_longer(cols = -time, names_to = "iteration", values_to = "c3")
  
  fitted_efficacy_draws <- left_join(left_join(c2b_draws, c2a_draws), c3_draws) %>%
    mutate(ve_c2b = (1-c2b/c3),
           ve_c2a = (1-c2a/c3))
  
  return(fitted_efficacy_draws)
  
  
}

run_median <- function(eir_med, age_at_enrollment,
                       site_peak, site_t_boost, seasonality_profile, 
                       itn_use_summary, vary_malariasimulation_parms = FALSE,
                       malariasimulation_parmset) {
  
  # Run median
  fitted_efficacy_med <- simulate_trial_hazards(eir = eir_med,
                                                age_at_enrollment = age_at_enrollment,
                                                gamma_llin = fix_gamma_llin,
                                                r_clin = fix_r_clin,
                                                peak = site_peak,
                                                rho = c(median(param_draws$rho1),
                                                        median(param_draws$rho2)),
                                                t_boost = site_t_boost,
                                                d_s = median(param_draws$d_s),
                                                d_l = median(param_draws$d_l),
                                                v_max = median(param_draws$v_max),
                                                alpha = median(param_draws$alpha),
                                                beta = median(param_draws$beta),
                                                age =  1:(365 * 12), s2 = 1.16,
                                                season = seasonality_profile, 
                                                vary_malariasimulation_parms = vary_malariasimulation_parms,
                                                malariasimulation_parmset = malariasimulation_parmset)
  
  # Take weighted average of clinical hazard between ITN users and non-users
  c2b_med <- fitted_efficacy_med$c2b_no_itn*
    (1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) +
    fitted_efficacy_med$c2b_itn*
    itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]
  
  c3_med <- fitted_efficacy_med$c3_no_itn*
    (1-itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]) +
    fitted_efficacy_med$c3_itn*
    itn_use_summary$itn_use[itn_use_summary$subgroup=="3"]
  
  output <- data.frame(c2b_med = c2b_med, c3_med = c3_med) %>%
    mutate(time = 1:length(c3_med),
           ve_c2b_med = (1-c2b_med/c3_med))
  return(output)
  
}

### 3. Project vaccine impact using model ##########################################################

# Vaccine efficacy is calculated over 12 months for Bagamoyo, Kilifi and Dande
# and 18 months for Bougouni and Nanoro

# Calculate vaccine efficacy over same period as in trial with median fitted EIR
# and including uncertainty in malaria model parameters

eir_draws_bagamoyo <- readRDS("efficacy_model/output/eir_draws_bagamoyo.RDS")
eir_draws_kilifi <- readRDS("efficacy_model/output/eir_draws_kilifi.RDS")
eir_draws_dande <- readRDS("efficacy_model/output/eir_draws_dande.RDS")
eir_draws_bougouni <- readRDS("efficacy_model/output/eir_draws_bougouni.RDS")
eir_draws_nanoro <- readRDS("efficacy_model/output/eir_draws_nanoro.RDS")

# Load malariasimulation parameter draws
urlfile <-'https://raw.githubusercontent.com/mrc-ide/malariasimulation/master/data-raw/parameter_draws.csv'
ms_draws<-read.csv(urlfile)
drawID <- readRDS("efficacy_model/data/malariasimulation_drawID.rds")   # select same 50 draws used in malariasimulation runs
ms_draws <- filter(ms_draws, draw %in% drawID)  
ms_draws <- ms_draws %>% pivot_wider(id_cols = draw, names_from = parameter, values_from = value)

## Bagamoyo

fitted_efficacy_median_bagamoyo <- list()
for(i in 1:nrow(ms_draws)) {
  fitted_efficacy_median_bagamoyo[[i]] <- run_median(eir_med=0.00706,
                                                     age_at_enrollment= age_to_model_bagamoyo$median,
                                                     site_peak=peak_bagamoyo,
                                                     site_t_boost=t_boost_bagamoyo,
                                                     seasonality_profile=seas_profile_bagamoyo,
                                                     itn_use_summary=itn_use_summary_bagamoyo,
                                                     vary_malariasimulation_parms = TRUE,
                                                     malariasimulation_parmset = ms_draws[i,])
}
fitted_efficacy_median_bagamoyo <- do.call("rbind", fitted_efficacy_median_bagamoyo)
fitted_efficacy_median_bagamoyo$drawID <- rep(drawID, each = nrow(fitted_efficacy_median_bagamoyo)/nrow(ms_draws))

## Kilifi

fitted_efficacy_median_kilifi <- list()
for(i in 1:nrow(ms_draws)) {
  fitted_efficacy_median_kilifi[[i]] <- run_median(eir_med=0.0068,  
                                                   age_at_enrollment=age_to_model_kilifi$median,
                                                   site_peak=peak_kilifi,
                                                   site_t_boost=t_boost_kilifi,
                                                   seasonality_profile=seas_profile_kilifi,
                                                   itn_use_summary=itn_use_summary_kilifi,
                                                   vary_malariasimulation_parms = TRUE,
                                                   malariasimulation_parmset = ms_draws[i,])
}
fitted_efficacy_median_kilifi <- do.call("rbind", fitted_efficacy_median_kilifi)
fitted_efficacy_median_kilifi$drawID <- rep(drawID, each = nrow(fitted_efficacy_median_kilifi)/nrow(ms_draws))

## Dande

fitted_efficacy_median_dande <- list()
for(i in 1:nrow(ms_draws)) {
  fitted_efficacy_median_dande[[i]] <- run_median(eir_med=0.0107,
                                                  age_at_enrollment=age_to_model_dande$median,
                                                  site_peak=peak_dande,
                                                  site_t_boost=t_boost_dande,
                                                  seasonality_profile=seas_profile_dande,
                                                  itn_use_summary=itn_use_summary_dande,
                                                  vary_malariasimulation_parms = TRUE,
                                                  malariasimulation_parmset = ms_draws[i,])
}
fitted_efficacy_median_dande <- do.call("rbind", fitted_efficacy_median_dande)
fitted_efficacy_median_dande$drawID <- rep(drawID, each = nrow(fitted_efficacy_median_dande)/nrow(ms_draws))

## Bougouni

fitted_efficacy_median_bougouni <- list()
for(i in 1:nrow(ms_draws)) {
  fitted_efficacy_median_bougouni[[i]] <- run_median(eir_med=0.0134,  
                                                     age_at_enrollment=age_to_model_bougouni$median,
                                                     site_peak=peak_bougouni,
                                                     site_t_boost=t_boost_bougouni,
                                                     seasonality_profile=seas_profile_bougouni,
                                                     itn_use_summary=itn_use_summary_bougouni,
                                                     vary_malariasimulation_parms = TRUE,
                                                     malariasimulation_parmset = ms_draws[i,])
}
fitted_efficacy_median_bougouni <- do.call("rbind", fitted_efficacy_median_bougouni)
fitted_efficacy_median_bougouni$drawID <- rep(drawID, each = nrow(fitted_efficacy_median_bougouni)/nrow(ms_draws))

# Nanoro
fitted_efficacy_median_nanoro <- list()
for(i in 1:nrow(ms_draws)) {
  fitted_efficacy_median_nanoro[[i]] <- run_median(eir_med=0.0563,
                                                   age_at_enrollment=age_to_model_nanoro$median,
                                                   site_peak=peak_nanoro,
                                                   site_t_boost=t_boost_nanoro,
                                                   seasonality_profile=seas_profile_nanoro,
                                                   itn_use_summary=itn_use_summary_nanoro,
                                                   vary_malariasimulation_parms = TRUE,
                                                   malariasimulation_parmset = ms_draws[i,])
}
fitted_efficacy_median_nanoro <- do.call("rbind", fitted_efficacy_median_nanoro)
fitted_efficacy_median_nanoro$drawID <- rep(drawID, each = nrow(fitted_efficacy_median_nanoro)/nrow(ms_draws))

model_efficacy_malariasimulation_draws <- rbind(
  cbind(site = "Bagamoyo", 
        fitted_efficacy_median_bagamoyo[fitted_efficacy_median_bagamoyo$time <= 365,]),
  cbind(site = "Kilifi", 
        fitted_efficacy_median_kilifi[fitted_efficacy_median_kilifi$time <= 365,]),
  cbind(site = "Dande", 
        fitted_efficacy_median_dande[fitted_efficacy_median_dande$time <= 365,]),
  cbind(site = "Bougouni", 
        fitted_efficacy_median_bougouni[fitted_efficacy_median_bougouni$time <= round(18*365/12),]),
  cbind(site = "Nanoro", 
        fitted_efficacy_median_nanoro[fitted_efficacy_median_nanoro$time <= round(18*365/12),])
)
#saveRDS(model_efficacy_malariasimulation_draws, "./efficacy_model/output/validation_efficacy_draws.RDS")

### 4. Validation plot (Figure 2) ##############################################

model_draws <- readRDS("./efficacy_model/output/validation_efficacy_draws.RDS")
model_draws$site <- factor(model_draws$site, levels = 
                                        c("Bagamoyo", "Kilifi", "Dande",
                                          "Bougouni", "Nanoro")) 

model_draws_summary <- model_draws %>%
  group_by(site) %>%
  summarise(ve_median = quantile(ve_c2b_med, 0.5),
            ve_lo = quantile(ve_c2b_med, 0.025),
            ve_hi = quantile(ve_c2b_med, 0.975)) 


ggplot() +
  geom_point(data = model_draws_summary, 
             aes(x=as.numeric(site)-0.2, y = ve_median*100, col = "Model"), 
             size = 4) +
  geom_errorbar(data = model_draws_summary, 
                aes(x=as.numeric(site)-0.2, y = ve_median*100, ymin = ve_lo * 100, ymax = ve_hi*100), 
                width = 0.2, col = "#E69F00") +
  scale_colour_manual(values = c("Model" = "#E69F00")) +
  ylab("Efficacy against clinical malaria (%)") +
  xlab("") +
  scale_x_continuous(breaks = c(1,2,3,4,5), labels = model_draws_summary$site) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(col = "", shape = "") +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(order = 1)) +
  ylim(0,100) +
  theme_classic() +
  coord_flip() + 
  theme(axis.text.y = element_text(hjust = 0), legend.position = "bottom",
        axis.text = element_text(size = 11)) 
#ggsave("fig2.png", width = 14, height = 10, units="cm")
