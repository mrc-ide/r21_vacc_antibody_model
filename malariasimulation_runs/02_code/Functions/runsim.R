# Set interventions and run malariasimulation ----------------------------------

runsim <- function(x){ # x = scenario #
  
  year <- 365
  month <- year / 12
  
  # read in selected scenario
  data <- readRDS("./malariasimulation_runs/03_output/run_parameters.rds")[x,]
  match <- readRDS("./malariasimulation_runs/03_output/EIRestimates.rds") |> select(-scenarioID)
  
  # EIR / prev match from "PfPR_EIR_match.R"
  data <- data |> left_join(match, by = c("drawID", "ID"))

  # EIR equilibrium ----------
  params <- set_equilibrium(unlist(data$params, recursive = F), as.numeric(data$starting_EIR))

  # run simulation ----------
  set.seed(123)

  output <- run_simulation(
    timesteps = data$warmup + data$sim_length,
    # correlations = correlations,
    parameters = params) |>

    # add vars to output
    mutate(ID = data$ID,
           scenario = data$scenarioID,
           drawID = data$drawID,
           EIR = data$starting_EIR,
           warmup = data$warmup,
           sim_length = data$sim_length,
           population = data$population,
           pfpr = data$pfpr,
           timestep = timestep - data$warmup,
           seasonality = data$seas_name,
           speciesprop = paste(data$speciesprop, sep = ",", collapse = ""),
           treatment = data$treatment,
           SMC = data$SMC,
           RTSS = data$RTSS,
           RTSScov = data$RTSScov,
           R21 = data$R21, 
           R21cov = data$R21cov,
           fifth = data$fifth) |>
    ungroup() |>
    filter(timestep > 0) |> # remove warmup period

    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) |>

    # keep only necessary variables
    dplyr::select(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop, treatment, SMC, RTSS, RTSScov, R21, R21cov, fifth,
                  starts_with("n_inc_severe"), starts_with("p_inc_severe"),
                  starts_with("n_pev"),
                  starts_with("n_inc"), starts_with("p_inc"),
                  starts_with("n_detect"), starts_with("p_detect"),
                  starts_with("n_"), -n_bitten, n_treated, n_infections) |>

    # take means of populations and sums of cases by month
    group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop, treatment, SMC, RTSS, RTSScov, R21, R21cov, fifth) |>

    mutate_at(vars(n_0_91.25:n_36500_73000, n_730_3650,
                   n_detect_730_3650, p_detect_730_3650), mean, na.rm = TRUE) |>
    mutate_at(vars(n_inc_severe_0_91.25:p_inc_clinical_36500_73000,
                   n_treated, n_infections), sum, na.rm = TRUE) |>

    dplyr::select(n_0_91.25:n_36500_73000,
                  n_inc_severe_0_91.25:p_inc_clinical_36500_73000,
                  n_detect_730_3650, p_detect_730_3650,
                  n_730_3650,
                  n_treated, n_infections) |>
    distinct()


  # save output ----------
  saveRDS(output, paste0("./malariasimulation_runs/03_output/HPC/general_", x, ".rds"))
}
