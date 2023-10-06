### Simulation function ########################################################

simulate_trial_hazards <- function(eir, age_at_enrollment, 
                                   gamma_llin, r_clin, 
                                   peak, rho, t_boost, d_s, d_l,
                                   v_max, alpha, beta, age, s2 = 1.16, 
                                   season = rep(1, 365), vary_malariasimulation_parms = FALSE,
                                   malariasimulation_parmset) {
  
  # Calculate vaccine efficacy by antibody titre
  # for those without (group 2b)/with a second booster (group 2a) 
  ab_2b <-  antibody_titre(1:(365*12), peak = peak[1:2], rho = rho[1:2],
                        t_boost = t_boost[1], d_s = d_s, d_l = d_l)
  vx_2b <- vaccine_efficacy(ab = ab_2b, v_max = v_max, alpha = alpha, beta = beta)
  
  ab_2a <-  antibody_titre(1:(365*12), peak = peak, rho = rho,
                           t_boost = t_boost, d_s = d_s, d_l = d_l)
  vx_2a <- vaccine_efficacy(ab = ab_2a, v_max = v_max, alpha = alpha, beta = beta)
  
  # Fitted model parameters
  if(vary_malariasimulation_parms == FALSE) {
    p <- load_parameter_set("Jamie_parameters.rds")
  } else if(vary_malariasimulation_parms == TRUE) {
    p <- load_parameter_set("Jamie_parameters.rds")
    p2 <- malariasimulation_parmset
    
    # Overwrite parameters in p with those from p2:
    p$ub <- p2$ub
    p$b0 <- p2$b0
    p$IB0 <- p2$ib0 
    p$kb <- p2$kb
    p$uc <- p2$uc 
    p$PM <- p2$pcm
    p$dm <- p2$rm  
    p$phi0 <- p2$phi0 
    p$phi1 <- p2$phi1 
    p$IC0 <- p2$ic0 
    p$kc <- p2$kc
  }

  p$s2 <- s2 
  
  # Heterogeneity
  n <- 10
  gh <- statmod::gauss.quad.prob(n = n, dist = "normal")
  zeta <- exp(-p$s2 * 0.5 + sqrt(p$s2) * gh$nodes)
  weight <- gh$weights
  
  # Maternally acquired immunity 
  ica20_no_itn <- get_maternal_start(eir = eir, gamma_llin = 1, season = season,
                                     rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                     b0 = p$b0, b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb,
                                     uc = p$uc, dc = p$dc)
  
  ica20_itn <- get_maternal_start(eir = eir, gamma_llin = gamma_llin, season = season,
                                  rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                  b0 = p$b0, b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb,
                                  uc = p$uc, dc = p$dc)
  
  icm_no_itn <- get_icm(age = age, ica20 = ica20_no_itn, pm = p$PM, dm = p$dm)
  icm_itn <- get_icm(age = age, ica20 = ica20_itn, pm = p$PM, dm = p$dm)
  
  # Simulate clinical hazard
  
  # For vaccine cohort 2b, without bednets
  c2b_no_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = 1, vx = vx_2b, season = season,
                                            icm = icm_no_itn, age = age, n = n, zeta = zeta, weight = weight,
                                            rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                            b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                            dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                            kc = p$kc)
  
  # For vaccine cohort 2b, with bednets
  c2b_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                         gamma_llin = gamma_llin, vx = vx_2b, season = season,
                                         icm = icm_itn, age = age, n = n, zeta = zeta, weight = weight,
                                         rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                         b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                         dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                         kc = p$kc)
  
  # For vaccine cohort 2a, without bednets
  c2a_no_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = 1, vx = vx_2a, season = season,
                                            icm = icm_no_itn, age = age, n = n, zeta = zeta, weight = weight,
                                            rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                            b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                            dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                            kc = p$kc)
  
  # For vaccine cohort 2a, with bednets
  c2a_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                         gamma_llin = gamma_llin, vx = vx_2a, season = season,
                                         icm = icm_itn, age = age, n = n, zeta = zeta, weight = weight,
                                         rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                         b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                         dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                         kc = p$kc)
  
  
  # For unvaccinated cohort, without bednets
  c3_no_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = 1, vx = rep(0, length(age)), 
                                            season = season, icm = icm_no_itn, age = age, 
                                            n = n, zeta = zeta, weight = weight, 
                                            rho = p$rho, a0 = p$a0, ub = p$ub,
                                            db = p$db, b0 = p$b0, b1 = p$b1, v1 = p$v1,
                                            ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                            phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc)
  
  # For unvaccinated cohort, with bednets
  c3_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                         gamma_llin = gamma_llin, vx = rep(0, length(age)), 
                                         season = season, icm = icm_itn, age = age, 
                                         n = n, zeta = zeta, weight = weight, 
                                         rho = p$rho, a0 = p$a0, ub = p$ub,
                                         db = p$db, b0 = p$b0, b1 = p$b1, v1 = p$v1,
                                         ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                         phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc)
  
  
  return(list(c2b_itn = c2b_itn, 
              c2b_no_itn = c2b_no_itn, 
              c2a_itn = c2a_itn, 
              c2a_no_itn = c2a_no_itn, 
              c3_itn = c3_itn,
              c3_no_itn = c3_no_itn))
  
}

simulate_control_hazards <- function(eir, age_at_enrollment, 
                                   gamma_llin, r_clin, age, s2 = 1.16, 
                                   season = rep(1, 365)) {
  
  # Fitted model parameters
  p <- load_parameter_set("Jamie_parameters.rds")
  p$s2 <- s2 
  
  # Heterogeneity
  n <- 10
  gh <- statmod::gauss.quad.prob(n = n, dist = "normal")
  zeta <- exp(-p$s2 * 0.5 + sqrt(p$s2) * gh$nodes)
  weight <- gh$weights
  
  # Maternally acquired immunity 
  ica20_no_itn <- get_maternal_start(eir = eir, gamma_llin = 1, season = season,
                                     rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                     b0 = p$b0, b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb,
                                     uc = p$uc, dc = p$dc)
  
  ica20_itn <- get_maternal_start(eir = eir, gamma_llin = gamma_llin, season = season,
                                  rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                  b0 = p$b0, b1 = p$b1, v1 = p$v1, ib0 = p$IB0, kb = p$kb,
                                  uc = p$uc, dc = p$dc)
  
  icm_no_itn <- get_icm(age = age, ica20 = ica20_no_itn, pm = p$PM, dm = p$dm)
  icm_itn <- get_icm(age = age, ica20 = ica20_itn, pm = p$PM, dm = p$dm)
  
  # Simulate clinical hazard
  
  # For unvaccinated cohort, without bednets
  c3_no_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = 1, vx = rep(0, length(age)), 
                                            season = season, icm = icm_no_itn, age = age, 
                                            n = n, zeta = zeta, weight = weight, 
                                            rho = p$rho, a0 = p$a0, ub = p$ub,
                                            db = p$db, b0 = p$b0, b1 = p$b1, v1 = p$v1,
                                            ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                            phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc)
  
  # For unvaccinated cohort, with bednets
  c3_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                         gamma_llin = gamma_llin, vx = rep(0, length(age)), 
                                         season = season, icm = icm_itn, age = age, 
                                         n = n, zeta = zeta, weight = weight, 
                                         rho = p$rho, a0 = p$a0, ub = p$ub,
                                         db = p$db, b0 = p$b0, b1 = p$b1, v1 = p$v1,
                                         ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                         phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc)

  
  return(list(c3_itn = c3_itn,
              c3_no_itn = c3_no_itn))
  
}
### Likelihood function ########################################################

efficacy_loglike <- function(params, data, misc){
  
  ### Simulate clinical hazard ###
  
  # Set up parameters
  EIR <- as.numeric(params["EIR"]) 
  gamma_llin <- as.numeric(params["gamma_llin"]) 
  r_clin <- as.numeric(params["r_clin"])
  v_max <- as.numeric(params["v_max"])
  alpha <- as.numeric(params["alpha"])
  beta <- as.numeric(params["beta"])
  rho1 <-  1 / (1 + exp(-as.numeric(params["rho1"])))
  rho2 <- 1 / (1 + exp(-as.numeric(params["rho2"])))
  rho <- c(rho1, rho2, rho2)
  d_s <- as.numeric(params["d_s"]) 
  d_l <- as.numeric(params["d_l"]) 
  
  seasonality <- misc$seasonality
  age_to_model <- misc$age_to_model
  peak <- misc$peak
  t_boost <- misc$t_boost
  
  # Simulate hazard and calculate likelihood for each arm, age group 
  # and bednet group individually
  
  likelihood_hazard_sim <- list()
  predicted_clinical_hazard_vacc2b_itn <- list()
  predicted_clinical_hazard_vacc2b_no_itn <- list()
  predicted_clinical_hazard_vacc2a_itn <- list()
  predicted_clinical_hazard_vacc2a_no_itn <- list()
  predicted_clinical_hazard_novacc_itn <- list()
  predicted_clinical_hazard_novacc_no_itn <- list()
  ll_individual_vacc2b_itn <- list()
  ll_individual_vacc2b_no_itn <- list()
  ll_individual_vacc2a_itn <- list()
  ll_individual_vacc2a_no_itn <- list()
  ll_individual_novacc_itn <- list()
  ll_individual_novacc_no_itn <- list()
  ll_vacc2b_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_vacc2b_no_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_vacc2a_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_vacc2a_no_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_novacc_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_novacc_no_itn <- rep(0, length(unique(age_to_model$age_group)))
  
  for (i in 1:length(unique(age_to_model$age_group))) {
    
    likelihood_hazard_sim[[i]] <- simulate_trial_hazards(eir = EIR, 
                                                         age_at_enrollment = age_to_model$median[i],
                                                         gamma_llin = gamma_llin, 
                                                         r_clin = r_clin,
                                                         peak = peak, 
                                                         rho = rho, 
                                                         t_boost = t_boost, 
                                                         d_s = d_s, d_l = d_l,
                                                         v_max = v_max,
                                                         alpha = alpha, 
                                                         beta = beta,
                                                         age =  1:(365 * 12), s2 = 1.16,
                                                         season = seasonality)
    
    
    predicted_clinical_hazard_vacc2b_itn[[i]] <- likelihood_hazard_sim[[i]]$c2b_itn
    predicted_clinical_hazard_vacc2b_no_itn[[i]] <- likelihood_hazard_sim[[i]]$c2b_no_itn
    predicted_clinical_hazard_vacc2a_itn[[i]] <- likelihood_hazard_sim[[i]]$c2a_itn
    predicted_clinical_hazard_vacc2a_no_itn[[i]] <- likelihood_hazard_sim[[i]]$c2a_no_itn
    predicted_clinical_hazard_novacc_itn[[i]] <- likelihood_hazard_sim[[i]]$c3_itn
    predicted_clinical_hazard_novacc_no_itn[[i]] <- likelihood_hazard_sim[[i]]$c3_no_itn
    
    ### Calculate log-likelihood ###
    
    # Log-likelihood for vaccinated individuals in group 2b (no second booster)
    
    # Participants with bednets
    ll_individual_vacc2b_itn[[i]] <- data$time_to_event1_itn %>%
      filter(age_group == age_to_model$age_group[i] &
               subgroup == "2b") %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_vacc2b_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_vacc2b_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_vacc2b_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_vacc2b_itn[[i]])[eventtime]))))
    
    ll_vacc2b_itn[i] <- sum(ll_individual_vacc2b_itn[[i]]$ll)
    
    # Participants without bednets
    ll_individual_vacc2b_no_itn[[i]] <- data$time_to_event1_no_itn %>%
      filter(age_group == age_to_model$age_group[i] &
               subgroup == "2b") %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_vacc2b_no_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_vacc2b_no_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_vacc2b_no_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_vacc2b_no_itn[[i]])[eventtime]))))
    
    
    ll_vacc2b_no_itn[i] <- sum(ll_individual_vacc2b_no_itn[[i]]$ll)
    
    # Log-likelihood for vaccinated individuals in group 2a (with second booster)
    
    # Participants with bednets
    ll_individual_vacc2a_itn[[i]] <- data$time_to_event1_itn %>%
      filter(age_group == age_to_model$age_group[i] &
               subgroup == "2a") %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_vacc2a_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_vacc2a_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_vacc2a_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_vacc2a_itn[[i]])[eventtime]))))
    
    
    ll_vacc2a_itn[i] <- sum(ll_individual_vacc2a_itn[[i]]$ll)
    
    # Participants without bednets
    ll_individual_vacc2a_no_itn[[i]] <- data$time_to_event1_no_itn %>%
      filter(age_group == age_to_model$age_group[i] &
               subgroup == "2a") %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_vacc2a_no_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_vacc2a_no_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_vacc2a_no_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_vacc2a_no_itn[[i]])[eventtime]))))
    
    
    ll_vacc2a_no_itn[i] <- sum(ll_individual_vacc2a_no_itn[[i]]$ll)
    
    # Log-likelihood for unvaccinated individuals 
    
    # Participants with bednets
    ll_individual_novacc_itn[[i]] <- data$time_to_event2_itn %>%
      filter(age_group == age_to_model$age_group[i]) %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_novacc_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_novacc_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_novacc_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_novacc_itn[[i]])[eventtime]))))
    
    
    ll_novacc_itn[i] <- sum(ll_individual_novacc_itn[[i]]$ll)
    
    # Participants without bednets
    ll_individual_novacc_no_itn[[i]] <- data$time_to_event2_no_itn %>%
      filter(age_group == age_to_model$age_group[i]) %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_novacc_no_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[eventtime]))))
    
    
    ll_novacc_no_itn[i] <- sum(ll_individual_novacc_no_itn[[i]]$ll)
    
  }
  
  # Total log likelihood for all participants across all age groups
  ll <- sum(ll_vacc2b_itn + ll_vacc2b_no_itn + ll_vacc2a_itn + ll_vacc2a_no_itn +
              ll_novacc_itn + ll_novacc_no_itn)
  
  return(ll)
}

validation_loglike <- function(params, data, misc){
  
  ### Simulate clinical hazard ###
  
  # Set up parameters
  EIR <- as.numeric(params["EIR"]) 
  gamma_llin <- misc$gamma_llin
  r_clin <- misc$r_clin
  # gamma_llin and r_clin are fixed
  
  seasonality <- misc$seasonality
  age_to_model <- misc$age_to_model
  
  # Simulate hazard and calculate likelihood for each arm, age group 
  # and bednet group individually
  
  likelihood_hazard_sim <- list()
  predicted_clinical_hazard_novacc_itn <- list()
  predicted_clinical_hazard_novacc_no_itn <- list()
  ll_individual_novacc_itn <- list()
  ll_individual_novacc_no_itn <- list()
  ll_novacc_itn <- rep(0, length(unique(age_to_model$age_group)))
  ll_novacc_no_itn <- rep(0, length(unique(age_to_model$age_group)))
  
  for (i in 1:length(unique(age_to_model$age_group))) {
    
    likelihood_hazard_sim[[i]] <- simulate_control_hazards(eir = EIR, 
                                                         age_at_enrollment = age_to_model$median[i],
                                                         gamma_llin = gamma_llin, 
                                                         r_clin = r_clin,
                                                         age =  1:(365 * 12), s2 = 1.16,
                                                         season = seasonality)
    
    predicted_clinical_hazard_novacc_itn[[i]] <- likelihood_hazard_sim[[i]]$c3_itn
    predicted_clinical_hazard_novacc_no_itn[[i]] <- likelihood_hazard_sim[[i]]$c3_no_itn
    
    ### Calculate log-likelihood ###
    
    # Log-likelihood for unvaccinated individuals 
    
    # Participants with bednets
    ll_individual_novacc_itn[[i]] <- data$time_to_event2_itn %>%
      filter(age_group == age_to_model$age_group[i]) %>%
      mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                   cumsum(predicted_clinical_hazard_novacc_itn[[i]])[(eventtime-1)]),
             ll = case_when(status == 1 ~ 
                              log(predicted_clinical_hazard_novacc_itn[[i]][eventtime]) +
                              log(exp(cumsum(predicted_clinical_hazard_novacc_itn[[i]])[(eventtime+7)]-
                                        subtr_cumsum)),
                            status == 0 ~ 
                              log(exp(-cumsum(predicted_clinical_hazard_novacc_itn[[i]])[eventtime]))))
    
    
    ll_novacc_itn[i] <- sum(ll_individual_novacc_itn[[i]]$ll)
    
    # Participants without bednets
    
    if(nrow(data$time_to_event2_no_itn) == 0) {
      ll_novacc_no_itn[i] <- 0
      
    } else if (nrow(data$time_to_event2_no_itn) != 0) {
      
      
      ll_individual_novacc_no_itn[[i]] <- data$time_to_event2_no_itn %>%
        filter(age_group == age_to_model$age_group[i]) %>%
        mutate(subtr_cumsum = ifelse(eventtime == 1, 0, 
                                     cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[(eventtime-1)]),
               ll = case_when(status == 1 ~ 
                                log(predicted_clinical_hazard_novacc_no_itn[[i]][eventtime]) +
                                log(exp(cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[(eventtime+7)]-
                                          subtr_cumsum)),
                              status == 0 ~ 
                                log(exp(-cumsum(predicted_clinical_hazard_novacc_no_itn[[i]])[eventtime]))))
      
      
      ll_novacc_no_itn[i] <- sum(ll_individual_novacc_no_itn[[i]]$ll)
      
    }
   
  }
  
  # Total log likelihood for all participants across all age groups
  ll <- sum(ll_novacc_itn + ll_novacc_no_itn)
  
  return(ll)
}

### Prior function #############################################################

efficacy_logprior <- function(params, misc){
  
  # extract parameter values
  EIR <- params["EIR"]                  # log-normal
  gamma_llin <- params["gamma_llin"]    # uniform
  r_clin <- params["r_clin"]            # gamma
  v_max <- params["v_max"]              # beta
  alpha <- params["alpha"]              # gamma
  beta <- params["beta"]                # gamma
  rho1 <- params["rho1"]                # normal
  rho2 <- params["rho2"]                # normal
  d_s <- params["d_s"]                  # log-normal
  d_l <- params["d_l"]                  # log-normal
  
  # Antibody priors from: data/ab_prior_parameters.csv
  
  # calculate log-prior
  ret <-  dlnorm(EIR, log(44.9/365), 0.5, log=TRUE) +  
    dunif(gamma_llin, min = 0, max = 10, log = TRUE) +  
    dgamma(r_clin, shape=5.230235, rate = 5.287793, log = TRUE) +
    dbeta(v_max, shape1=19.836, shape2=2.204, log = TRUE) +
    dgamma(alpha, shape = 4.129436, rate = 4.064203, log = TRUE) +
    dunif(beta, min = 0, max = 11159, log = TRUE) +
    dnorm(rho1, mean = 0.807167555599538, sd = 0.0825488238386109, log = TRUE) +
    dnorm(rho2, mean = 0.0714033701850845, sd = 0.0787771092892667, log = TRUE) +
    dlnorm(d_s, meanlog = 3.79850436548426, sdlog = 0.0468110331879774, log = TRUE) +
    dlnorm(d_l, meanlog = 6.27915084227264, sdlog = 0.0756914748978481, log = TRUE)
  
  # return
  return(ret)
}

validation_logprior <- function(params, misc){
  
  # extract parameter values
  EIR <- params["EIR"]                  # log-normal

  # calculate log-prior
  ret <-  dunif(EIR, min = 0, max = 100/365, log=TRUE)
  
  # return
  return(ret)
}
