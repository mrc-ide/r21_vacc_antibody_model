# Average clinical hazard
get_clinical_hazard <- function(eir, age_at_enrollment, gamma_llin, vx,
                                season, icm, age, n, zeta, weight,
                                rho, a0, ub, db, b0, b1, v1, ib0, kb, uc, dc,
                                phi0, phi1, ic0, kc){
  
  # Time 0 in the model now corresponds to data: 28 days after the 3rd dose,
  # which is also the start of vaccination effect coming into place and when
  # seasonality curve begins
  
  # Model run starts at age 0 - age is equivalent to time since birth.
  # Vaccine effect only comes into place at our assumed start/enrollment date. 
  # Clinical hazard is returned only for timepoints after enrollment.
  
  # Adapted to have flexible seasonality over the years
  # Repeats the same seasonality pattern irrespective of its length
  # Ensures we always start at the same point in the seasonality curve
  # for all ages at enrollment
  time <- round(age - age_at_enrollment)
  vx_shift <- c(rep(0, sum(time < 0)), vx)[age]
  seasonality <-  season[(time %% length(season)) + 1]
  
  clinical_hazard <- rep(0, 365 * 10)
  for(h in 1:n){
    # Exposure to infectious bites
    epsilon <- get_epsilon(age = age, eir = eir * zeta[h], pa = rho, a0 = a0,
                           gamma_llin = gamma_llin)
    
    # Immunity against malaria infection
    ib <- acquire_immunity(exposure = epsilon, u = ub, d = db)
    
    # Probability infection
    b <- get_b(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
    
    # Hazard of infection
    infection_hazard <- get_infection_hazard(epsilon = epsilon, b = b,
                                             seasonality = seasonality,
                                             vaccine_efficacy = vx_shift)
    
    # Immunity against clinical malaria
    ica <- acquire_immunity(exposure = infection_hazard, u = uc, d = dc)
    
    # Probability an infection is clinical
    phi <- get_phi(ica = ica, icm = icm, phi0 = phi0, phi1 = phi1,
                   ic0 = ic0, kc = kc)
    
    index_period <- which(time > 0 & time <= (10 * 365))
    # Hazard of clinical infection
    clinical_hazard <- clinical_hazard +
      get_clinical_hazard_no_het(infection_hazard = infection_hazard,
                                 clinical_probability = phi)[index_period] * weight[h]
  }
  return(clinical_hazard)
}
# Estimate the average level of immunity in a 20 year old women
# to inform the maternal immunity
get_maternal_start <- function(eir, gamma_llin, 
                         season, age,
                         rho, a0, ub, db, b0, b1, v1, ib0, kb, uc, dc){
  
  age <- 0:(365*20)
  vx <- rep(0, length(age))
  seasonality <-  season[(age %% length(season)) + 1]
  
  # Exposure to infectious bites
  epsilon <- get_epsilon(age = age, eir = eir, pa = rho, a0 = a0,
                         gamma_llin = gamma_llin)
  
  # Immunity against malaria infection
  ib <- acquire_immunity(exposure = epsilon, u = ub, d = db)
  
  # Probability infection
  b <- get_b(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
  
  # Hazard of infection
  infection_hazard <- get_infection_hazard(epsilon = epsilon, b = b,
                                           seasonality = seasonality,
                                           vaccine_efficacy = vx)
  
  # Immunity against clinical malaria
  ica <- acquire_immunity(exposure = infection_hazard, u = uc, d = dc)
  return(ica[length(age)])
}