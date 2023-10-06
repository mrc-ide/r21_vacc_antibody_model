# Run scenarios through malariasimulation on HPC

source("./malariasimulation_runs/02_code/data_and_libraries.R")

# Set up your job --------------------------------------------------------------
# MODEL set-up ----
# year
year <- 365

# population
population <- 200000

# run time
warmup <- 21 * year       # needs to be multiple of 3 for ITN distribution
sim_length <- 15 * year   # value > 0 

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- readRDS("./malariasimulation_runs/03_output/parameter_draws.rds")$draw
drawID_R21 <- c(1:50) 

draws <- tibble(drawID, drawID_R21) 

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65, 0.01)

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
seas_name <- "highly seasonal"
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- tibble(seasonality, seas_name)

seas_name <- "seasonal"
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- tibble(seasonality, seas_name)

seas_name <- "perennial"
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
s3 <- tibble(seasonality, seas_name)

stable <- bind_rows(s1, s2, s3)

# vectors
# list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

# INTERVENTIONS ----
# treatment coverage
treatment <- c(0.45) 

# SMC: 0, 0.75
SMC <- c(0, 0.75) 

# RTS,S: none, EPI, SV, hybrid
RTSS <- c("none", "EPI", "SV", "hybrid", "EPI12") 

# RTS,S coverage
RTSScov <- c(0, 0.80) 

# R21: none, EPI, SV, hybrid
R21 <- c("none", "EPI", "SV", "hybrid") 

# R21 coverage
R21cov <- c(0, 0.80) 

# adding a fifth RTS,S dose: 0, 1
fifth <- c(0, 1)  

interventions <- crossing(treatment, SMC, RTSS, RTSScov, R21, R21cov, fifth)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, draws) |>
  mutate(ID = paste(pfpr, seas_name, drawID, sep = "_")) 

# remove non-applicable scenarios
combo <- combo |>
  filter(!(seas_name == "highly seasonal" & SMC == 0)) |>
  filter(!(seas_name %in% c("seasonal", "perennial") & SMC != 0)) |>
  filter(!(RTSS == "none" & RTSScov > 0)) |> 
  filter(!(RTSScov == 0 & RTSS %in% c("EPI", "SV", "hybrid", "EPI12"))) |> 
  filter(!(R21 == "none" & R21cov > 0)) |> 
  filter(!(R21cov == 0 & R21 %in% c("EPI", "SV", "hybrid"))) |> 
  filter(!(seas_name == "perennial" & (RTSS == "SV" | RTSS == "hybrid" | R21 == "SV" | R21 == "hybrid"))) |>
  filter(!(RTSS %in% c("EPI", "SV", "hybrid", "EPI12") & R21 %in% c("EPI", "SV", "hybrid"))) |>
  filter(!(seas_name == "highly seasonal"))

# put variables into the same order as function arguments
combo <- combo |> 
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         R21,               # R21 strategy
         R21cov,            # R21 coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID,            # parameter draw no.
         drawID_R21         # R21 parameter draw no.
  ) |> as.data.frame()

saveRDS(combo, "./malariasimulation_runs/03_output/runscenarios.rds")

# generate parameter list in malariasimulation format
source("./malariasimulation_runs/02_code/Functions/generate_params.R")

generate_params("./malariasimulation_runs/03_output/runscenarios.rds", # file path to pull
                "./malariasimulation_runs/03_output/run_parameters.rds")      # file path to push
# glance
paramlist <- readRDS("./malariasimulation_runs/03_output/run_parameters.rds")


# Run tasks --------------------------------------------------------------------
x <- c(1:nrow(combo)) # number of runs

# define all combinations of scenarios and draws
index <- tibble(x = x)

source("./malariasimulation_runs/02_code/Functions/runsim.R")

# Example run
runsim(x=index$x[1])
