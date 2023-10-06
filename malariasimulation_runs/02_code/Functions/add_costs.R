# input: data frame
# process: reads in dataframe, calculates cost of interventions and cost-effectiveness
# output: data frame with additional variables for cost and cost-effectiveness

add_costs <- function(x # dataframe to read in and process
                      ){

  SMCcost <- 0.9075                     

  rtss_cost_per_dose <- c(3.79, 12.01)      
  r21_cost_per_dose <- c(2.69, 3.79, 5.24)  
  delivery_cost_seasonal <- 3.35            
  delivery_cost_epi <- 1.33
  delivery_cost_hybrid <- 2.09

  # create combinations of dose cost and delivery cost
  rtsscost_df_seasonal <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_seasonal)
  R21cost_df_seasonal <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_seasonal)
  rtsscost_df_epi <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_epi)
  R21cost_df_epi <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_epi)
  rtsscost_df_hybrid <- expand_grid(cost_per_dose = rtss_cost_per_dose, delivery_cost = delivery_cost_hybrid)
  R21cost_df_hybrid <- expand_grid(cost_per_dose = r21_cost_per_dose, delivery_cost = delivery_cost_hybrid)
  novax_cost <- expand_grid(cost_per_dose = 0, delivery_cost = 0)
  
  RDT <- 0.46 + (0.46 * 0.15)          
  AL_adult <- 0.3 * 24                  
  AL_child <- 0.3 * 12                  

  outpatient <- 1.87                    
  inpatient <- 8.71                     

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
           cost_clinical = ((cases - severe_cases - u5_cases) * treatment * TREATcost_adult) * .77 +
             ((u5_cases - u5_severe) * treatment * TREATcost_child) * .77,
           # severe treatment
           cost_severe = (severe_cases * treatment * SEVcost_adult) * .77 +
             (u5_severe * treatment * SEVcost_child) * .77,
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
             cost_SMC + cost_vax)

  return(dalyoutput_cost)

}
