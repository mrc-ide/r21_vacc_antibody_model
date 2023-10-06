#' Vaccine efficacy model
#' 
#' Following White et al (2015) Defaults are from White Table 3 (5-17 months)
#'
#' @param ab Antibody titre
#' @param v_max Maximum efficacy against infection
#' @param alpha Shape parameter of dose–response curve
#' @param beta Scale parameter of dose–response curve
#'
#' @return Vector of vaccine efficacy.
vaccine_efficacy <- function(ab, v_max = 0.93, alpha = 0.74, beta = 99.2,
                             avidity = 0, t_boost = 304){
  
  v_max * (1 - (1 / (1 + ((ab / beta) ^ alpha))))
  
}
