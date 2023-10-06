# Cost effectiveness -----------------------------------------------------------

cost_effectiveness <- function(x, y){

  # create file index
  index <- seq(x, y, 1)

  # function for processing HPC data
  cost_effectiveness_b <- function(x){ # input = index of file to process

  # process HPC output: condense by age group over simulation time ----
  output <- HPC_processing(x)

  # run mortality and DALY functions ----
  output <- mortality_rate(output)
  output <- outcome_uncertainty(output)
  output <- daly_components(output)

  # condense into one line per run ----
  output <- output |>
    select(-inc, -sev, -mortality_rate) |> # get rid of rate vars

    # create vars for childhood cases
    mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '91.25-1825'), n, 0), # u5 denominator
           n_91.25_1825 = ifelse(age == '91.25-1825', n, 0), # SMC denominator
           u5_cases = ifelse(age %in% c('0-91.25', '91.25-1825'), cases, 0),
           u5_severe = ifelse(age %in% c('0-91.25', '91.25-1825'), severe_cases, 0),
           u5_deaths = ifelse(age %in% c('0-91.25', '91.25-1825'), deaths, 0),
           u5_dalys = ifelse(age %in% c('0-91.25', '91.25-1825'), daly, 0)) |>

    mutate_at(vars(n, n_0_1825, n_91.25_1825,
                   u5_cases, u5_severe, u5_deaths, u5_dalys,
                   inc_clinical, inc_severe,
                   cases, cases_lower, cases_upper, severe_cases,
                   deaths, deaths_lower, deaths_upper,
                   yll:daly_lower), sum, na.rm = T) |>  # condense outputs over all ages in population
    select(-age, -age_upper, -age_lower) |>
    distinct()

  # add costs (intervention-specific and total) ----
  output <- add_costs(output)

  print(x)

  return(output)

  }

  # run cost_effectiveness function
  dalyoutput <- map_dfr(index, cost_effectiveness_b)

  saveRDS(dalyoutput, paste0('./malariasimulation_runs/03_output/HPC_processing/run_', x, '_', y, '.rds'))



}

