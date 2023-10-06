# Set up draws of malariasimulation parameters
# 50 draws are used for uncertainty runs

source("./malariasimulation_runs/02_code/data_and_libraries.R")

# get draws --------------------------------------------------------------------

# choose a random sample of 50 draws to use in malariasimulation::set_parameter_draw()
set.seed(123)
rsample <- sample(x = c(1:1000), size = 50, replace = FALSE, prob = NULL)
rsample <- tibble(draw = rsample)

saveRDS(rsample, "./malariasimulation_runs/03_output/parameter_draws.rds")

# read in R21 parameters - pulled from a private repo
r21 <- read_tsv("./malariasimulation_runs/01_data/R21_params.txt", col_names = F) # github .csv links in .txt file

ab_params <- read_csv(url(paste0(r21[2, 1]))) |> # antibody
  pivot_wider(names_from = par, values_from = c(mu, sd)) |>
  dplyr::select(-draw)

e_params <- read_csv(url(paste0(r21[4, 1]))) # efficacy

head(rsample); head(e_params); head(ab_params)

r21_params <- bind_cols(rsample, e_params, ab_params)

saveRDS(r21_params, "./malariasimulation_runs/03_output/r21_draws.rds")

