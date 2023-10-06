#' Antibody titre model
#' 
#' Following White et al (2015). The model is generalised for N boosters, where,
#' for parameters peak and rho a vector can be supplied. The first value of the
#' vector corresponds to the parameter after the initial 3-dose schedule, with
#' subsequent values defining the value following each boost. Half lives of the
#' short- and long-lived components of the anti body response are assumed to be
#' constant across doses/boosters. Defaults are from White Table 3 (5-17 months)
#' 
#'
#' @param t Vector of times (days since dose 3)
#' @param peak Peak titres.
#' @param rho Proportion of response that is short-lived. 
#' @param t_boost timing of booster doses (days).
#' @param d_s Half-life of short-lived component of antibody response.
#' @param d_l Half-life of long-lived component of antibody response.
#'
#' @return A vector of antibody titres.
antibody_titre <- function(t, peak = c(589, 264), rho = c(0.88, 0.7),
                           t_boost = 548, d_s = 45, d_l = 591){
  r_s <- log(2) / d_s
  r_l <- log(2) / d_l
  t_boost <- c(0, t_boost)
  
  ab <- rep(NA, length(t))
  for(i in seq_along(t_boost)){
    t_cur <- t - t_boost[i]
    index <- t >= t_boost[i]
    ab[index] <- peak[i] * (rho[i] * exp(-r_s * t_cur[index]) + (1 - rho[i]) * exp(-r_l * t_cur[index]))
  }
  return(ab)
}

# cpp implementation of the above
cpp_function('doubles antibody_titre_cpp(integers t, doubles peak, doubles rho, doubles t_boost, double d_s, double d_l) {
  int n = t.size();
  double r_s = log(2) / d_s;
  double r_l = log(2) / d_l;
  
  writable::doubles ab(n);
  
  int boost_level = 0;
  for(int i = 0; i < n; ++i) {
      if(t[i] == t_boost[boost_level]){
        boost_level ++;
      }
    int t_cur = t[i];
    if(boost_level > 0){
      t_cur = t[i] - t_boost[boost_level - 1];
    }
    
    ab[i] = peak[boost_level] * (rho[boost_level] * exp(-r_s * t_cur) + (1 - rho[boost_level]) * exp(-r_l * t_cur));
  }
  return ab;
}')
