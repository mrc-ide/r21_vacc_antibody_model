# input: index of malariasimulation HPC run
# process: reads in raw HPC output, adds R21 doses, condenses output over simulation length
# output: data frame with one row per age group

HPC_processing <- function(x){

  # read in specified rds file
    output <- readRDS(paste0("./malariasimulation_runs/03_output/HPC/general_", format(x, scientific = FALSE), ".rds"))

    # extract type of RTS,S / R21 intervention
    RTSS = output$RTSS[1]
    R21 = output$R21[1]
    fifth = output$fifth[1]
    
    # add vaccine doses
    if(RTSS == "none" | R21 == "none"){
      output <- output |> rowwise() |>
        mutate(dose1 = 0,
               dose2 = 0,
               dose3 = 0,
               dose4 = 0,
               dose5 = 0) |>
        ungroup()
    }
    
    if(RTSS == "EPI" | RTSS == "EPI12" | R21 == "EPI" | RTSS == "hybrid" | R21 == "hybrid"){
      output <- output |> rowwise() |>
        mutate(dose1 = n_pev_epi_dose_1,
               dose2 = n_pev_epi_dose_2,
               dose3 = n_pev_epi_dose_3,
               dose4 = n_pev_epi_booster_1) |>
        ungroup()
    }
    
    if(RTSS == "SV" | R21 == "SV"){
      output <- output |> rowwise() |>
        mutate(dose1 = n_pev_mass_dose_1,
               dose2 = n_pev_mass_dose_2,
               dose3 = n_pev_mass_dose_3,
               dose4 = n_pev_mass_booster_1) |>
        ungroup()
    }
    
    if(fifth == 1){
      if(RTSS == "EPI" | RTSS == "EPI12" | R21 == "EPI" | RTSS == "hybrid" | R21 == "hybrid"){
        output <- output |> rowwise() |>
          mutate(dose5 = n_pev_epi_booster_2) |>
          ungroup()
      }
      
      if(RTSS == "SV" | R21 == "SV"){
        output <- output |> rowwise() |>
          mutate(dose5 = n_pev_mass_booster_2) |>
          ungroup()
      }
    }
    
    # summarize data over the first 15 years (mult of 3 for ITNs)
    sim_length <- output$sim_length[1] / 365
    
    ## DISCOUNTING
    discount_function <- function(x, rate, year) {
      x * (1/(1+rate)^(year-1))
    }
    # Discount cases and severe cases, then take sum over 15 years
    
    discounted_output <- output |> filter(year <= sim_length) |> # first 15 years
      mutate(n_91.25_1825_discounted = discount_function(n_91.25_1825, rate = 0.03, year = year)) |> # for calculation of SMC cost
      mutate_at(vars(n_inc_severe_0_91.25:p_inc_severe_36500_73000), 
                ~ discount_function(.x, rate = 0.03, year = year)) |>  # discounting
      mutate_at(vars(n_inc_clinical_0_91.25:p_inc_clinical_36500_73000), 
                ~ discount_function(.x, rate = 0.03, year = year)) |>
      mutate_at(vars(dose1:dose5), 
                ~ discount_function(.x, rate = 0.03, year = year)) |>   # for calculation of vaccine cost
      mutate_at(vars(n_0_91.25:n_36500_73000), mean, na.rm = TRUE) |>   # mean of n in each age group
      mutate_at(vars(n_inc_severe_0_91.25:dose5), sum, na.rm = TRUE) |> # sum of cases and vax doses
      mutate_at(vars(n_91.25_1825_discounted), mean, na.rm = TRUE) |>   # mean of n in each age group
      mutate_at(vars(n_inc_severe_0_91.25:p_inc_severe_36500_73000), round , 0) |> 
      mutate_at(vars(n_inc_clinical_0_91.25:p_inc_clinical_36500_73000), round, 0) |>
      mutate_at(vars(dose1:dose5), round, 0) |>
      select(-month, -year) |>
      distinct() |>
      
      # calculate outputs by age
      dplyr::select(ID:n_36500_73000, n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                    n_inc_severe_0_91.25:n_inc_severe_36500_73000, n_detect_730_3650:n_730_3650,
                    n_treated, n_infections, dose1:dose5, n_91.25_1825_discounted) |>
      
      # moving from wide to long age groups
      pivot_longer(cols = c(n_0_91.25:n_36500_73000,
                            n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                            n_inc_severe_0_91.25:n_inc_severe_36500_73000),
                   names_to = c('age'), values_to = c('value')) |>
      mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
             inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
             inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
             age = gsub('n_inc_clinical_', '', age),                         # combining age vars
             age = gsub('n_inc_severe_', '', age),
             age = gsub('n_', '', age),
             age = gsub('_', '-', age)) |>
      group_by(age) |>
      select(-value) |>
      mutate_at(vars(n:inc_severe), sum, na.rm = TRUE) |> # consolidate
      distinct() |> ungroup()
    
    discounted_output <- discounted_output %>%
      select(ID:fifth, age, dose1:dose5, n_91.25_1825_discounted, inc_clinical, inc_severe) %>%
      rename(inc_clinical_discounted = inc_clinical, inc_severe_discounted = inc_severe,
             dose1_discounted = dose1, dose2_discounted = dose2, dose3_discounted = dose3,
             dose4_discounted = dose4, dose5_discounted = dose5)
    
    # Outputs without discounting
    output <- output |> filter(year <= sim_length) |> # first 15 years
      mutate_at(vars(n_0_91.25:n_36500_73000), mean, na.rm = TRUE) |>   # mean of n in each age group
      mutate_at(vars(n_inc_severe_0_91.25:dose5), sum, na.rm = TRUE) |> # sum of cases and vax doses
      select(-month, -year) |>
      distinct() |>
      
      # calculate outputs by age
      dplyr::select(ID:n_36500_73000, n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                    n_inc_severe_0_91.25:n_inc_severe_36500_73000, n_detect_730_3650:n_730_3650,
                    n_treated, n_infections, dose1:dose5) |>
      
      # moving from wide to long age groups
      pivot_longer(cols = c(n_0_91.25:n_36500_73000,
                            n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                            n_inc_severe_0_91.25:n_inc_severe_36500_73000),
                   names_to = c('age'), values_to = c('value')) |>
      mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
             inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
             inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
             age = gsub('n_inc_clinical_', '', age),                         # combining age vars
             age = gsub('n_inc_severe_', '', age),
             age = gsub('n_', '', age),
             age = gsub('_', '-', age)) |>
      group_by(age) |>
      select(-value) |>
      mutate_at(vars(n:inc_severe), sum, na.rm = TRUE) |> # consolidate
      distinct() |> ungroup()
    
    output <- left_join(output, discounted_output)
    
    output <- output |>
      separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) |>
      mutate(age_lower = as.numeric(age_lower)/365,
             age_upper = as.numeric(age_upper)/365,
             inc = inc_clinical / n,
             sev = inc_severe / n,
             cases = inc_clinical,
             severe_cases = inc_severe,
             cases_discounted = inc_clinical_discounted,
             severe_cases_discounted = inc_severe_discounted)
    


    output$file <- paste0("./malariasimulation_runs/03_output/HPC/general_", x, ".rds")


  return(output)

}


