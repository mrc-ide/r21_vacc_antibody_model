# PfPR EIR match ---------------------------------------------------------------

PRmatch <- function(x, y){ # x = scenario # , y = parameter draw #
  
  # read in selected scenario
  data <- readRDS('./malariasimulation_runs/03_output/baseline_parameters.rds')[x,] 
  
  # choose a parameter set from baseline scenarios
  p <- unlist(data$params, recursive = F)
  
  # calibration ref: https://mrc-ide.github.io/cali/articles/Basic_calibration.html
  # define target: PfPR2-10 value
  target <- data$pfpr
  
  # write function for time points to average PfPR - years 4-6
  # 6 years is enough to match by PfPR - 21 years is needed to match by c inc
  year <- 365
  p$timesteps <- 9 * year   # simulation run time = 9 years
  
  summary_mean_pfpr_2_10_6y9y <- function(x){
    
    x$year <- ceiling(x$timestep / year)
    x <- x |> filter(year >= 7)
    prev_2_10 <- mean(x$n_detect_730_3650 / x$n_730_3650)
    return(prev_2_10)
    
  }
  
  # run calibration model
  set.seed(123)
  out <- cali::calibrate(parameters = p,
                         target = target,
                         summary_function = summary_mean_pfpr_2_10_6y9y,
                         tolerance = 0.02,
                         low = 0.001,
                         high = 1500) 
  
  # store init_EIR results as an .rds file to be read in later
  PR <- data.frame(scenarioID = x, drawID = y)
  PR$starting_EIR <- out
  PR$ID <- data$ID
  
  saveRDS(PR, paste0('./malariasimulation_runs/03_output/PR_EIR/PRmatch_draws_', data$ID, '.rds'))
  
}

