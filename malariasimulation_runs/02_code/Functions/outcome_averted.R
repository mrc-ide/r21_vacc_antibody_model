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
           
           cost_clinical_baseline = cost_clinical,
           cost_severe_baseline = cost_severe,
           cost_vax_baseline = cost_vax,
           
           cost_total_baseline = cost_total,
           cost_total_u5_baseline = cost_total_u5,
           
           # discounted outcomes
           
           daly_baseline_discounted = daly_discounted,
           cases_baseline_discounted = cases_discounted,
           severe_baseline_discounted = severe_cases_discounted,
           deaths_baseline_discounted = deaths_discounted,
           
           u5_dalys_baseline_discounted = u5_dalys_discounted,
           u5_cases_baseline_discounted = u5_cases_discounted,
           u5_severe_baseline_discounted = u5_severe_discounted,
           u5_deaths_baseline_discounted = u5_deaths_discounted,
           
           cost_clinical_baseline_discounted = cost_clinical_discounted,
           cost_severe_baseline_discounted = cost_severe_discounted,
           cost_vax_baseline_discounted = cost_vax_discounted,
           
           cost_total_baseline_discounted = cost_total_discounted,
           cost_total_u5_baseline_discounted = cost_total_u5_discounted) |>
    
    select(file, ID, scenario, drawID, daly_baseline, cases_baseline,
           severe_baseline, deaths_baseline, u5_dalys_baseline,
           u5_cases_baseline, u5_severe_baseline, u5_deaths_baseline, 
           cost_clinical_baseline, cost_severe_baseline, cost_vax_baseline,
           cost_total_baseline, cost_total_u5_baseline,
           daly_baseline_discounted, cases_baseline_discounted,
           severe_baseline_discounted, deaths_baseline_discounted, u5_dalys_baseline_discounted,
           u5_cases_baseline_discounted, u5_severe_baseline_discounted, u5_deaths_baseline_discounted, 
           cost_clinical_baseline_discounted, cost_severe_baseline_discounted, cost_vax_baseline_discounted,
           cost_total_baseline_discounted, cost_total_u5_baseline_discounted) |>
    
    distinct()
  
  
  # separate out non baseline scenarios and merge
  base_IDs <- none$file
  
  scenarios <- output |> filter(!(file %in% base_IDs)) |>
    left_join(none |> select(-file, -scenario) |> distinct(), by = c("ID", "drawID")) |>
    
    # calculate cost-effectiveness
    mutate(
      cases_averted = (cases_baseline - cases),
      dalys_averted = (daly_baseline - daly),
      incremental_cost = (cost_total - cost_total_baseline),
      cost_saved = (cost_clinical_baseline + cost_severe_baseline) - (cost_clinical + cost_severe),
      cost_vax_diff = (cost_vax - cost_vax_baseline),
      CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
      CE_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
      CE_case = (cost_total - cost_total_baseline) / (cases_baseline - cases),
      CE_u5_case = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases),
      
      # discounted outcomes
      
      cases_averted_discounted = (cases_baseline_discounted - cases_discounted),
      dalys_averted_discounted = (daly_baseline_discounted - daly_discounted),
      incremental_cost_discounted = (cost_total_discounted - cost_total_baseline_discounted),
      cost_saved_discounted = (cost_clinical_baseline_discounted + cost_severe_baseline_discounted) - 
        (cost_clinical_discounted + cost_severe_discounted),
      cost_vax_diff_discounted = (cost_vax_discounted - cost_vax_baseline_discounted),
      CE_discounted = (cost_total_discounted - cost_total_baseline_discounted) / 
        (daly_baseline_discounted - daly_discounted),
      CE_case_discounted = (cost_total_discounted - cost_total_baseline_discounted) / 
        (cases_baseline_discounted - cases_discounted)
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
