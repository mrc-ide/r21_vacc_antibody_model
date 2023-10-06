# Generate parameters for malariasimulation ------------------------------------

generate_params <- function(inputpath,   # path to input scenarios
                            outputpath){ # path where output file will be stored
  
  # read in dataframe of all scenario combinations
  scenarios <- readRDS(inputpath)
  
  # generate parameters
  generate_params2 <- function(x){ # x = scenario number
    
    # pull one scenario at a time
    data <- scenarios[x, ]
    
    # assign values
    population = data$population
    seasonality = data$seasonality
    seas_name = data$seas_name
    pfpr = data$pfpr
    warmup = data$warmup
    sim_length = data$sim_length
    speciesprop = data$speciesprop
    treatment = data$treatment
    SMC = data$SMC
    RTSS = data$RTSS
    RTSScov = data$RTSScov
    R21 = data$R21
    R21cov = data$R21cov
    fifth = data$fifth
    ID = data$ID
    drawID = data$drawID
    drawID_R21 = data$drawID_R21
    
    year <- 365
    month <- year / 12
    
    # starting parameters ----------
    params <- get_parameters(list(
      human_population = population,
      model_seasonality = TRUE,
      # rainfall fourier parameters
      g0 = unlist(seasonality)[1],
      g = unlist(seasonality)[2:4],
      h = unlist(seasonality)[5:7],
      individual_mosquitoes = FALSE))
    
    # parameter draws ----------
    params <- set_parameter_draw(
      parameters = params,
      draw = drawID)
    
    # R21 vax and booster profiles ----------
    # < pull values from R21 parameter draws
    r21_params <- readRDS("./malariasimulation_runs/01_data/r21_draws.rds")[drawID_R21,]  
    
    r21_profile <- rtss_profile
    
    r21_profile$vmax = r21_params$v_max
    r21_profile$alpha = r21_params$alpha
    r21_profile$beta = r21_params$beta
    r21_profile$cs = c(r21_params$mu_r21_cs, r21_params$sd_r21_cs)
    r21_profile$rho = c(r21_params$mu_r21_rho, r21_params$sd_r21_rho)
    r21_profile$ds = c(r21_params$mu_r21_ds, r21_params$sd_r21_ds)
    r21_profile$dl = c(r21_params$mu_r21_dl, r21_params$sd_r21_dl)
    
    r21_booster_profile <- r21_profile
    r21_booster_profile$cs = c(r21_params$mu_r21_cs_boost, r21_params$sd_r21_cs_boost)
    r21_booster_profile$rho = c(r21_params$mu_r21_rho_boost, r21_params$sd_r21_rho_boost)
    
    r21_booster_profile2 <- r21_booster_profile
    r21_booster_profile2$cs = c(r21_params$mu_r21_cs_boost2, r21_params$sd_r21_cs_boost2)
    
    # outcome definitions ----------
    # incidence for every 5 year age group
    params$clinical_incidence_rendering_min_ages = c(0, 0.25, seq(5, 100, 5)) * year
    params$clinical_incidence_rendering_max_ages = c(0.25, seq(5, 100, 5), 200) * year
    params$severe_incidence_rendering_min_ages = c(0, 0.25, seq(5, 100, 5)) * year
    params$severe_incidence_rendering_max_ages = c(0.25, seq(5, 100, 5), 200) * year
    
    # prevalence 2-10 year olds
    params$prevalence_rendering_min_ages = 2 * year
    params$prevalence_rendering_max_ages = 10 * year
    
    # demography ----------
    africa_demog <- read.csv("./malariasimulation_runs/01_data/ssa_demography_2021.csv")
    ages <- round(africa_demog$age_upper * year) # top of age bracket
    deathrates <- africa_demog$mortality_rate / 365   # age-specific death rates
    
    params <- set_demography(
      params,
      agegroups = ages,
      timesteps = 0,
      deathrates = matrix(deathrates, nrow = 1))
    
    # vectors ----------
    params <- set_species(
      parameters = params,
      species = list(arab_params, fun_params, gamb_params),
      proportions = unlist(speciesprop))
    
    # proportion of bites taken in bed for each species
    # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
    # proportion of bites taken indoors for each species
    params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020
    
    # treatment ----------
    if (treatment > 0) {
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      # AL default, SP: https://doi.org/10.1016/S2214-109X(22)00416-8 supplement
      params$drug_prophylaxis_scale <- c(10.6, 39.34) 
      params$drug_prophylaxis_shape <- c(11.3, 3.40) 
      
      params <- set_clinical_treatment(
        parameters = params,
        drug = 1,
        timesteps = c(1),
        coverages = c(treatment)
      )  }
    
    # SMC ----------
    smc_timesteps <- 0
    
    if (SMC > 0 & seas_name == "seasonal") {
      peak <- peak_season_offset(params)
      # 5 doses, centered around peak
      first <- round(warmup + c(peak + c(-2, -1, 0, 1, 2) * month), 0)
      firststeps <- sort(rep(first, sim_length/year))
      yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
      timesteps <- yearsteps + firststeps
      
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      # AL default, SP: https://doi.org/10.1016/S2214-109X(22)00416-8 supplement
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_smc(
        parameters = params,
        drug = 2,
        timesteps = sort(timesteps),
        coverages = rep(SMC, length(timesteps)),
        min_age = round(0.25 * year),
        max_age = round(5 * year))
      
      smc_timesteps <- params$smc_timesteps - warmup
    }
    
    if (SMC > 0 & seas_name == "highly seasonal") {
      peak <- peak_season_offset(params)
      # 4 doses, centered around peak
      first <- round(c(peak + c(-1, 0, 1, 2) * month), 0)
      firststeps <- sort(rep(first, (warmup + sim_length)/year))
      yearsteps <- rep(c(0, seq(year, (warmup + sim_length) - year, year)), length(first))
      timesteps <- yearsteps + firststeps
      
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      # AL default, SP: https://doi.org/10.1016/S2214-109X(22)00416-8 supplement
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_smc(
        parameters = params,
        drug = 2,
        timesteps = sort(timesteps),
        coverages = rep(SMC, length(timesteps)),
        min_ages = rep((0.25 * year), length(timesteps)),
        max_ages = rep((5 * year), length(timesteps)))
      
      # var for outputting to check SMC timings are correct
      smc_timesteps <- params$smc_timesteps - warmup
    }
    
    # RTS,S EPI ----------
    if (RTSS == "EPI") {
      params$pev_doses <- round(c(0, 1.5 * month, 3 * month))
      boosters <- if(fifth == 0) round(c(18 * month)) else round(c(18 * month, 30 * month))
      
      params <- set_pev_epi(
        parameters = params,
        profile = rtss_profile,
        coverages = RTSScov,
        timesteps = warmup,
        age = round(6 * month), # 6, 7,5, 9 months
        min_wait = 0,
        booster_timestep = boosters,    # 27 months
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = rep(list(rtss_booster_profile), length(boosters)),
        seasonal_boosters = FALSE
      )  }
    
    # RTS,S EPI 12 month booster ----------
    if (RTSS == "EPI12") {
      params$pev_doses <- round(c(0, 1 * month, 2 * month))
      boosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month))
      
      # update booster to have same effect as dose 3 per Thompson et al. 2022
      rtss_booster_profile2 <- rtss_booster_profile
      rtss_booster_profile2$cs <- c(6.37008, 0.35)
      
      params <- set_pev_epi(
        parameters = params,
        profile = rtss_profile,
        coverages = RTSScov,
        timesteps = warmup,
        age = round(6 * month), # 6, 7, 8 months
        min_wait = 0,
        booster_timestep = boosters,    # 20 months
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = rep(list(rtss_booster_profile2), length(boosters)),
        seasonal_boosters = FALSE)  
      
      }    
    
    # RTS,S SV ----------
    rtss_mass_timesteps <- 0
    
    if (RTSS == "SV") {
      peak <- peak_season_offset(params)
      if(seas_name == "highly seasonal"){
        first <- round(warmup + (peak - month * 3.5), 0)}
      if(seas_name == "seasonal"){
        first <- round(warmup + (peak - month * 5.5), 0)}
      timesteps <- c(first, first+seq(year, sim_length, year))
      params$pev_doses <- round(c(0, 1 * month, 2 * month))
      
      boosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month))
      
      # update booster to have same effect as dose 3 per Thompson et al. 2022
      rtss_booster_profile2 <- rtss_booster_profile
      rtss_booster_profile2$cs <- c(6.37008, 0.35)
      
      params <- set_mass_pev(
        parameters = params,
        profile = rtss_profile,
        timesteps = timesteps,
        coverages = rep(RTSScov,length(timesteps)),
        min_ages = round(5 * month),
        max_ages = round(17 * month),
        min_wait = 0,
        booster_timestep = boosters,
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = rep(list(rtss_booster_profile2), length(boosters)))
      
      # var for outputting to check RTS,S timings are correct
      rtss_mass_timesteps <- params$rtss_mass_timesteps - warmup
      
    }
    
    # RTS,S hybrid ----------
    if (RTSS == "hybrid") {
      params$pev_doses <- round(c(0, 1.5 * month, 3 * month))
      
      peak <- peak_season_offset(params)
      
      if(seas_name == "highly seasonal"){
        boost <- round((peak - month * 1.5), 0)}
      if(seas_name == "seasonal"){
        boost <- round((peak - month * 3.5), 0)}
      
      boosters <- if(fifth == 0) boost else c(boost, boost + year)
      
      # update booster to have same effect as dose 3 per Thompson et al. 2022
      rtss_booster_profile2 <- rtss_booster_profile
      rtss_booster_profile2$cs <- c(6.37008, 0.35)
      
      params <- set_pev_epi(
        parameters = params,
        profile = rtss_profile,
        coverages = RTSScov,
        timesteps = warmup,
        age = round(6 * month),  # 6, 7,5, 9 months
        min_wait = round(6 * month),
        booster_timestep = boosters,
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = rep(list(rtss_booster_profile2), length(boosters)),
        seasonal_boosters = TRUE)

    }
    
    
    # R21 EPI ----------
    if (R21 == "EPI") {
      params$pev_doses <- round(c(0, 1 * month, 2 * month))
      boosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month))
      booster_profile <- if(fifth == 0) list(r21_booster_profile) else list(r21_booster_profile, r21_booster_profile2)
      
      
      params <- set_pev_epi(
        parameters = params,
        profile = r21_profile,
        coverages = R21cov,
        timesteps = warmup,
        age = round(6 * month), # 6, 7, 8 months
        min_wait = 0,
        booster_timestep = boosters,    # 20 months
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = booster_profile,
        seasonal_boosters = FALSE
        )  }
    
    # R21 SV ----------
    rtss_mass_timesteps <- 0
    
    if (R21 == "SV") {
      peak <- peak_season_offset(params)
      if(seas_name == "highly seasonal"){
        first <- round(warmup + (peak - month * 3.5), 0)}
      if(seas_name == "seasonal"){
        first <- round(warmup + (peak - month * 5.5), 0)}
      timesteps <- c(first, first+seq(year, sim_length, year))
      params$pev_doses <- round(c(0, 1 * month, 2 * month))
      
      boosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month))
      booster_profile <- if(fifth == 0) list(r21_booster_profile) else list(r21_booster_profile, r21_booster_profile2)
      
      params <- set_mass_pev(
        parameters = params,
        profile = r21_profile,
        timesteps = timesteps,
        coverages = rep(R21cov,length(timesteps)),
        min_ages = round(5 * month),
        max_ages = round(17 * month),
        min_wait = 0,
        booster_timestep = boosters,
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = booster_profile)
      
      # var for outputting to check RTS,S timings are correct
      rtss_mass_timesteps <- params$rtss_mass_timesteps - warmup
      
    }
    
    # R21 hybrid ----------
    if (R21 == "hybrid") {
      params$pev_doses <- round(c(0, 1 * month, 2 * month))
      
      peak <- peak_season_offset(params)
      
      if(seas_name == "highly seasonal"){
        boost <- round((peak - month * 1.5), 0)}
      if(seas_name == "seasonal"){
        boost <- round((peak - month * 3.5), 0)}
      
      boosters <- if(fifth == 0) boost else c(boost, boost + year)
      booster_profile <- if(fifth == 0) list(r21_booster_profile) else list(r21_booster_profile, r21_booster_profile2)
      
      params <- set_pev_epi(
        parameters = params,
        profile = r21_profile,
        coverages = R21cov,
        timesteps = warmup,
        age = round(6 * month),  # 6, 7, 8 months 
        min_wait = round(6 * month),
        booster_timestep = boosters,
        booster_coverage = rep(.80, length(boosters)),
        booster_profile = booster_profile,
        seasonal_boosters = TRUE)
  
    }
  
    # save as data.frame
    data$params <- list(params)
    data$scenarioID <- x
    
    # print count
    print(x)
    
    return(data)
    
  }
  
  # loop through function to generate parameters one by one
  output <- map_dfr(1:nrow(scenarios), generate_params2)
  
  
  # save output ----------
  saveRDS(output, outputpath)
  
}
