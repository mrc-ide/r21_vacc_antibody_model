# Note - due to checks on down boost this needs to be run with the full t timeseries
antibody_titre <- function(t, phase, peak1, peak2, peak3, duration1, duration2, rho1, rho2, t_boost1 = 364, t_boost2 = 729){
  r1 <- log(2) / duration1
  r2 <- log(2) / duration2
  
  t[phase == 2] <-  t[phase == 2] - t_boost1
  t[phase == 3] <-  t[phase == 3] - t_boost2
  
  rho <- rep(rho1, length(t))
  rho[phase == 2] <- rho2
  rho[phase == 3] <- rho2
  
  peak <- rep(peak1, length(t))
  peak[phase == 2] <- peak2
  peak[phase == 3] <- peak3

  ab <- peak * ((rho * exp(-r1 * t)) + ((1 - rho) *  exp(-r2 * t)))
  
  # Check for down-boosts
  last_phase_1 <- tail(ab[phase == 1], 1)
  if(last_phase_1 > peak2){
    peak[phase == 2] <- last_phase_1
  }
  last_phase_2 <- tail(ab[phase == 2], 1)
  if(last_phase_2 > peak3){
    peak[phase == 3] <- last_phase_2
  }
  ab <- peak * ((rho * exp(-r1 * t)) + ((1 - rho) *  exp(-r2 * t)))
  
  return(ab)
}

# cpp implementation of the above
# cpp_function('doubles antibody_titre_cpp(integers t, doubles peak, doubles rho, doubles t_boost, double d_s, double d_l) {
#   int n = t.size();
#   double r_s = log(2) / d_s;
#   double r_l = log(2) / d_l;
#   
#   writable::doubles ab(n);
#   
#   int boost_level = 0;
#   for(int i = 0; i < n; ++i) {
#       if(t[i] == t_boost[boost_level]){
#         boost_level ++;
#       }
#     int t_cur = t[i];
#     if(boost_level > 0){
#       t_cur = t[i] - t_boost[boost_level - 1];
#     }
#     
#     ab[i] = peak[boost_level] * (rho[boost_level] * exp(-r_s * t_cur) + (1 - rho[boost_level]) * exp(-r_l * t_cur));
#   }
#   return ab;
# }')
