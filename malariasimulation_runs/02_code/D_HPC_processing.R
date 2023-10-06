# Process malariasimulation runs and add costing output on HPC

source("./malariasimulation_runs/02_code/data_and_libraries.R")

source("./malariasimulation_runs/02_code/Functions/HPC_processing.R")     # one line per age group
source("./malariasimulation_runs/02_code/Functions/deaths_dalys.R")       # calc deaths & dalys
source("./malariasimulation_runs/02_code/Functions/add_costs.R")          # add cost estimates
source("./malariasimulation_runs/02_code/Functions/outcome_averted.R")    # cases & dalys averted
source("./malariasimulation_runs/02_code/Functions/cost_effectiveness.R") # calc CE

# Run tasks --------------------------------------------------------------------
combo <- readRDS("./malariasimulation_runs/03_output/runscenarios.rds")

index <- seq(1, nrow(combo), 1)

# break tasks to run in groups of 1,000
cost_effectiveness(1,1000)
cost_effectiveness(1001,2000)
cost_effectiveness(2001,3000)
cost_effectiveness(3001,4000)
cost_effectiveness(4001,5000)
cost_effectiveness(5001,6000)
cost_effectiveness(6001,7000)
cost_effectiveness(7001,8000)
cost_effectiveness(8001,9000)
cost_effectiveness(9001,10000)
cost_effectiveness(10001,11000)
cost_effectiveness(11001,12000)
cost_effectiveness(12001,13200)

# Results ----------------------------------------------------------------------
# read in all parameter draw runs and process
files <- list.files(path = "./malariasimulation_runs/03_output/HPC_processing/", pattern = "run*", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
dalyoutput <- rbindlist(dat_list, fill = TRUE, idcol = "identifier")

# save output
saveRDS(dalyoutput, "./malariasimulation_runs/03_output/dalyoutput_draws.rds")

# calculate cases / DALYs averted
source("./malariasimulation_runs/02_code/Functions/outcome_averted.R")
output <- outcome_averted(dalyoutput)

# save final output
saveRDS(output, paste0(path, "./malariasimulation_runs/03_output/scenarios_draws.rds"))

# check intervention classifications
table(output$intervention)
table(dalyoutput$RTSS); table(dalyoutput$R21)
table(dalyoutput$RTSS, dalyoutput$R21)
summary(output$CE)

