# input: data frame
# process: reads in dataframe, calculates cost-effectiveness of interventions compared to baseline
# output: data frame with additional variables for cost-effectiveness per case or DALY averted

outcome_averted <- function(x # dataframe to read in and process
                            ){

  output <- x

  # separate out baseline scenarios
  none <- output |>
    filter(RTSS == "none" & R21 == "none" & # filter out interventions
             (SMC == 0 | (seasonality == "highly seasonal"))) |>
    rename(daly_baseline = daly,
           cases_baseline = cases,
           severe_baseline = severe_cases,
           deaths_baseline = deaths,

           u5_dalys_baseline = u5_dalys,
           u5_cases_baseline = u5_cases,
           u5_severe_baseline = u5_severe,
           u5_deaths_baseline = u5_deaths,

           cost_total_baseline = cost_total,
           cost_total_u5_baseline = cost_total_u5) |>

    select(file, ID, scenario, drawID, daly_baseline, cases_baseline,
           severe_baseline, deaths_baseline, u5_dalys_baseline,
           u5_cases_baseline, u5_severe_baseline, u5_deaths_baseline, 
           cost_total_baseline, cost_total_u5_baseline) |>
    
    distinct()
    

  # separate out non baseline scenarios and merge
  base_IDs <- none$file

  scenarios <- output |> filter(!(file %in% base_IDs)) |>
    left_join(none |> select(-file, -scenario) |> distinct(), by = c("ID", "drawID")) |>

    # calculate cost-effectiveness
    mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
           CE_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
           CE_case = (cost_total - cost_total_baseline) / (cases_baseline - cases),
           CE_u5_case = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases)
    ) |>

    # assign scenarios
    mutate(intervention = case_when(

      (SMC==0 | (seasonality=="highly seasonal")) & RTSS=="none" & R21=="none" ~ "none",
      (SMC==0 | (seasonality=="highly seasonal")) & RTSS=="EPI" ~ "RTS,S age-based",
      (SMC==0 | (seasonality=="highly seasonal")) & RTSS=="EPI12" ~ "RTS,S age-based (12m booster)",
      (SMC==0 | (seasonality=="highly seasonal")) & RTSS=="SV" ~ "RTS,S seasonal",
      (SMC==0 | (seasonality=="highly seasonal")) & RTSS=="hybrid" ~ "RTS,S hybrid",
      (SMC==0 | (seasonality=="highly seasonal")) & R21=="EPI" ~ "R21 age-based",
      (SMC==0 | (seasonality=="highly seasonal")) & R21=="SV" ~ "R21 seasonal",
      (SMC==0 | (seasonality=="highly seasonal")) & R21=="hybrid" ~ "R21 hybrid")) |>

    mutate(intervention_f = factor(intervention, levels = c("none", "RTS,S age-based", "RTS,S age-based (12m booster)", "RTS,S seasonal", "RTS,S hybrid", "R21 age-based", "R21 seasonal", "R21 hybrid"))) |>

    mutate(rank = as.numeric(intervention_f))

  return(scenarios)

}
