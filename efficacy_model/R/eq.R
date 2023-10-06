
# Immunity function. Acquisition of immunity against infection (IB) and clinical
# malaria (ICA). White et al (S26).
cpp_function('doubles acquire_immunity(doubles exposure, double u, double d) {
  int n = exposure.size();
  writable::doubles immunity(n);
  immunity[0] = 0;
  for(int i = 1; i < n; ++i) {
    immunity[i] = immunity[i - 1] + (exposure[i] / (exposure[i] * u + 1)) - (immunity[i - 1] / d); 
  }
  return immunity;
}')

get_icm <- function(age, ica20, pm, dm){
  pm * ica20 * exp(- age / dm)
}

get_epsilon <- function(age, eir, pa, a0, gamma_llin = 1){
  gamma_llin * eir *  (1 - pa * exp(-age / a0))
}

# get_phi <- function(ica, icm, phi0, phi1, ic0, kc){
#   phi0 * (phi1 + ((1 - phi1) / (1 + ((ica + icm) / ic0)^ kc)))
# }

cpp_function('doubles get_phi(doubles ica, doubles icm, double phi0, double phi1, double ic0, double kc) { 
  int n = ica.size(); 
  writable::doubles out(n); 
  for(int i = 0; i < n; ++i) { 
    out[i] = phi0 * (phi1 + ((1 - phi1) / (1 + pow(((ica[i] + icm[i]) / ic0), kc)))); 
  } 
  return out; 
}')

# get_b <- function(ib, b0, b1, ib0, kb){
#   b0 * (b1 + ((1 - b1) / (1 + (ib / ib0)^ kb)))
# }

cpp_function('doubles get_b(doubles ib, double b0, double b1, double ib0, double kb) {
  int n = ib.size();
  writable::doubles out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = b0 * (b1 + ((1 - b1) / (1 + pow((ib[i] / ib0), kb))));
  }
  return out;
}')

get_infection_hazard <- function(epsilon, b, seasonality = NULL, vaccine_efficacy = NULL){
  ih <- epsilon * b
  
  if(!is.null(seasonality)){
    ih <- ih * seasonality
  }
  
  if(!is.null(vaccine_efficacy)){
    ih <- ih * (1 - vaccine_efficacy)
  }
  
  return(ih)
}

get_clinical_hazard_no_het <- function(infection_hazard, clinical_probability){
  infection_hazard * clinical_probability
}