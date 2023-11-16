# input: data frame
# process: reads in dataframe, calculates cost of interventions and cost-effectiveness
# output: data frame with additional variables for cost and cost-effectiveness

add_costs <- function(x # dataframe to read in and process
                      ){

  SMCcost <- 0.9075                          # SMC $0.9075 per dose (no SMC in runs!)
  
  rtss_cost_per_dose <- c(3.79, 12.01)       # RTS,S per dose $3, $9.3 + consumables cost
  r21_cost_per_dose <- c(2.7, 3.97, 5.23)    # R21 per dose $2, $3, $4 + consumables cost
  
  delivery_cost_seasonal <- 3.75             # same for RTS,S and R21
  delivery_cost_epi <- 1.48
  delivery_cost_hybrid <- 2.36
  
  # create combinations of dose cost and delivery cost
  rtsscost_df_seasonal <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_seasonal)
  R21cost_df_seasonal <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_seasonal)
  rtsscost_df_epi <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_epi)
  R21cost_df_epi <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_epi)
  rtsscost_df_hybrid <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_hybrid)
  R21cost_df_hybrid <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_hybrid)
  novax_cost <- expand_grid(cost_per_dose = 0, delivery_cost = 0)
  
  RDT <- 0.48 + (0.48 * 0.15)            # RDT $0.48 unit cost + (unit cost * 15% delivery markup)
  AL_adult <- 0.45 * 24                  # $7.2 clinical treatment cost ($0.3 * 24 doses)
  AL_child <- 0.45 * 12                  # $3.6 clinical treatment cost ($0.3 * 12 doses)
  
  outpatient <- 2.08                     # (Median WHO Choice cost for SSA)
  inpatient <- 24.96                     # (Median WHO Choice cost for SSA, assuming average duration of stay of 3 days)
  
  # clinical: RDT cost + Drug course cost + facility cost (outpatient)
  TREATcost_adult <- RDT + AL_adult + outpatient
  TREATcost_child <- RDT + AL_child + outpatient
  
  # severe: RDT cost + Drug course cost + facility cost (inpatient)
  SEVcost_adult <- RDT + AL_adult + inpatient
  SEVcost_child <- RDT + AL_child + inpatient
  
  
  population <- x$population[1]
  sim_length <- x$sim_length[1]
  
  # prepare to add costs to dataset
  dalyoutput_cost <- x
  
  if(x$R21 == "none" & x$RTSS == "SV"){ # merge in RTSS costing dataframe
    dalyoutput_cost <- dalyoutput_cost |> merge(rtsscost_df_seasonal) |> ungroup() |> rowwise()
  }
  if(x$R21 == "none" & x$RTSS %in% c("EPI", "EPI12")){ # merge in RTSS costing dataframe
    dalyoutput_cost <- dalyoutput_cost |> merge(rtsscost_df_epi) |> ungroup() |> rowwise()
  }
  if(x$R21 == "none" & x$RTSS == "hybrid"){ # merge in RTSS costing dataframe
    dalyoutput_cost <- dalyoutput_cost |> merge(rtsscost_df_hybrid) |> ungroup() |> rowwise()
  }
  
  if(x$RTSS == "none" & x$R21 == "SV"){ # merge in R21 costing dataframe
    dalyoutput_cost <- dalyoutput_cost |>  merge(R21cost_df_seasonal) |> ungroup() |> rowwise()
  }
  if(x$RTSS == "none" & x$R21 == "EPI"){ # merge in R21 costing dataframe
    dalyoutput_cost <- dalyoutput_cost |>  merge(R21cost_df_epi) |> ungroup() |> rowwise()
  }
  if(x$RTSS == "none" & x$R21 == "hybrid"){ # merge in R21 costing dataframe
    dalyoutput_cost <- dalyoutput_cost |>  merge(R21cost_df_hybrid) |> ungroup() |> rowwise()
  }
  
  if(x$RTSS == "none" & x$R21 == "none"){ # merge in R21 costing dataframe
    dalyoutput_cost <- dalyoutput_cost |>  merge(novax_cost) |> ungroup() |> rowwise()
  }
    
  
  # create cost variables
  # 77% of treatment costs are from the public sector (DHS, SSA)
  dalyoutput_cost <- dalyoutput_cost |>
    
    mutate(smc_timesteps = case_when(SMC > 0 & seasonality == "highly seasonal" ~ 4 * (sim_length / 365),
                                     SMC > 0 & seasonality == "seasonal" ~ 5 * (sim_length / 365),
                                     SMC == 0 ~ 0)) |>
    
    mutate(# non-severe treatment
      cost_clinical = (((cases - severe_cases - u5_cases) * treatment * TREATcost_adult) +
                         ((u5_cases - u5_severe) * treatment * TREATcost_child)) * 0.77,
      # severe treatment
      cost_severe = ((severe_cases * treatment * SEVcost_adult) +
                       (u5_severe * treatment * SEVcost_child)) * 0.77,
      # SMC
      cost_SMC = n_91.25_1825 * SMC * SMCcost * smc_timesteps,
      # RTSS / R21
      cost_vax = (dose1 + dose2 + dose3 + dose4 + dose5) * (cost_per_dose + delivery_cost),
      
      # TOTAL
      cost_total = cost_clinical + cost_severe + cost_SMC + cost_vax,
      
      # cost just among children
      cost_total_u5 =
        # cost clinical
        ((u5_cases-u5_severe) * treatment * TREATcost_child)*.77 +
        # cost severe
        (u5_severe * treatment * SEVcost_child)*.77 +
        # cost SMC and cost VAX are all among children
        cost_SMC + cost_vax,
      
      # discounted costs
      
      cost_clinical_discounted = (((cases_discounted - severe_cases_discounted - u5_cases_discounted) * treatment * TREATcost_adult) +
                                    ((u5_cases_discounted - u5_severe_discounted) * treatment * TREATcost_child)) * 0.77,
      # severe treatment
      cost_severe_discounted = ((severe_cases_discounted * treatment * SEVcost_adult) +
                                  (u5_severe_discounted * treatment * SEVcost_child)) * 0.77,
      # SMC
      cost_SMC_discounted = n_91.25_1825_discounted * SMC * SMCcost * smc_timesteps,
      # RTSS / R21
      cost_vax_discounted = (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted) * 
        (cost_per_dose + delivery_cost),
      
      # TOTAL
      cost_total_discounted = cost_clinical_discounted + cost_severe_discounted + cost_SMC_discounted + 
        cost_vax_discounted,
      
      # cost just among children
      cost_total_u5_discounted =
        # cost clinical
        ((u5_cases_discounted-u5_severe_discounted) * treatment * TREATcost_child) * 0.77 +
        # cost severe
        (u5_severe_discounted * treatment * SEVcost_child) * 0.77 +
        # cost SMC and cost VAX are all among children
        cost_SMC_discounted + cost_vax_discounted)
  
  return(dalyoutput_cost)
}
