//lognormal_partial_pooling.stan
functions {
  real abm(real t, real phase, real peak1, real peak2, real peak3, real duration1, real duration2, real rho1, real rho2) {
    
    real ab;
    real ab_low;
    
    real r1 = log(2) / duration1;
    real r2 = log(2) / duration2;
    
    real peak = peak1;
    real tm = t;
    real rho = rho1;
    
    
    if(phase == 2){
      // Estimate titre 1 timestep before boost (to check model doesn't boost lower)
      tm = 364;
      ab_low = peak * ((rho * exp(-r1 * tm)) + ((1 - rho) *  exp(-r2 * tm)));
      // Set time, peak and rho for current datum
      tm = t - 364;
      peak = peak2;
      rho = rho2;
    }
    if(phase == 3){
      // Estimate titre 1 timestep before boost (to check model doesn't boost lower)
      tm = 729;
      rho = rho2;
      peak = peak2;
      ab_low = peak * ((rho * exp(-r1 * tm)) + ((1 - rho) *  exp(-r2 * tm)));
      // Set time, peak and rho for current datum
      tm = t - 729;
      peak = peak3;
      rho = rho2;
    }
    // If boost is lower than previous timestep estimate adjust
    if(peak < ab_low){
      peak = ab_low;
    }

    ab = peak * ((rho * exp(-r1 * tm)) + ((1 - rho) *  exp(-r2 * tm)));
    return ab;
  }
}
data {
  int<lower=1> N; // number of individuals
  int<lower=1> R; // number of records
  int<lower=1,upper=N> individual[R]; // individual index
  real<lower=0> ab[R]; // antibody titre
  int<lower=0> t[R]; // time
  int<lower=1> phase[R]; // Phase (primary, booster 1, booster 2)
  real<lower=0> peak1[N]; // Primary peak
  real<lower=0> peak2[N]; // Booster 1 peak
  real<lower=0> peak3[N]; // Booster 2 peak
  
  // Information on priors
  real prior_duration1[2];
  real prior_duration2[2];
  real prior_logit_rho1[2];
  real prior_logit_rho2[2];
  real prior_sigma[2];
  real prior_re_sd[2];
}
parameters {
  vector<lower=0>[N] duration1;
  vector<lower=0>[N] duration2;
  vector[N] logit_rho1;
  vector[N] logit_rho2;
  real<lower=0> sigma;
  real<lower=0> global_duration1;
  real<lower=0> global_duration1_sd;
  real<lower=0> global_duration2;
  real<lower=0> global_duration2_sd;
  real global_rho1;
  real<lower=0> global_rho1_sd;
  real global_rho2;
  real<lower=0> global_rho2_sd;
}
transformed parameters {
  vector[N] rho1 = inv_logit(logit_rho1);
  vector[N] rho2 = inv_logit(logit_rho2);
}
model {
  vector[R] predicted_titre;
  for(i in 1:R){
    int p = individual[i];
    predicted_titre[i] = log(abm(t[i], phase[i], peak1[p], peak2[p], peak3[p], duration1[p], duration2[p], rho1[p], rho2[p]));
  } 
  
  // Priors
  duration1 ~ lognormal(log(prior_duration1[1]), prior_duration1[2]);
  duration2 ~ lognormal(log(prior_duration2[1]), prior_duration2[2]);
  logit_rho1 ~ normal(prior_logit_rho1[1], prior_logit_rho1[2]);
  logit_rho2 ~ normal(prior_logit_rho2[1], prior_logit_rho2[2]);
  sigma ~ lognormal(log(prior_sigma[1]), prior_sigma[2]);
  
  // Hyper-priors
  global_duration1_sd ~ gamma(prior_re_sd[1], prior_re_sd[2]);
  global_duration2_sd ~ gamma(prior_re_sd[1], prior_re_sd[2]);
  global_rho1_sd ~ gamma(prior_re_sd[1], prior_re_sd[2]);
  global_rho2_sd ~ gamma(prior_re_sd[1], prior_re_sd[2]);
  
  // Random effects
  duration1 ~ lognormal(log(global_duration1), global_duration1_sd);
  duration2 ~ lognormal(log(global_duration2), global_duration2_sd);
  logit_rho1 ~ normal(global_rho1, global_rho1_sd);
  logit_rho2 ~ normal(global_rho2, global_rho2_sd);
  
  // Likelihood
  ab ~ lognormal(predicted_titre, sigma);
}
