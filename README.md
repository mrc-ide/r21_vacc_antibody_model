# R21/Matrix-M vaccine antibody and vaccine efficacy model

Code to 
- estimate the association between antibody dynamics and vaccine efficacy against clinical malaria
of the R21/Matrix-M vaccine using Phase 2 trial data
- project the public health impact and cost-effectiveness of the R21/Matrix-M vaccine

As shown in: Schmit, Topazian et al. The public health impact and cost-effectiveness of the R21/Matrix-M malaria vaccine: a mathematical modelling study.

## Directory

```
.
├── ab_model                           # Scripts to estimate the dynamics of anti-CSP antibody titres
|   ├── 1_fit_ab_model.R                 # Fit antibody model to immunogenicity data
|   ├── 2_visualise_ab_model.R           # Visualise and save parameter estimates from antibody model fit
|   ├── ab_model.stan                    # Antibody model in STAN
|   ├── output                           # Model outputs
|   ├── R                                # Functions to run and fit model
|   |   ├── antibody_model.R               # Antibody titre dynamics model
|   |   ├── draw_priors.R                  # Function to draw prior samples
|   |   ├── stan_utils.R                   # Helper functions
├── efficacy_model                     # Scripts to estimate the relationship between antibody titres and vaccine efficacy
|   ├── 3_fit_efficacy_model.R           # Fit efficacy model to case incidence data
|   ├── 4_model_validation.R             # Validate model against Phase 3 trial efficacy
|   ├── data                             # Data files
|   ├── output                           # Model outputs
|   ├── R                                # Functions to run and fit model
|   |   ├── antibody_model.R               # Antibody titre dynamics model
|   |   ├── efficacy_fitting.R             # Functions for model fitting
|   |   ├── eq.R                           # Function for immunity and infection hazard
|   |   ├── hazard.R                       # Case incidence model
|   |   ├── season.R                       # Seasonality function
|   |   ├── utils.R                        # Helper functions
|   |   ├── vaccine_efficacy_model.R       # Model for relationship between antibody titres and vaccine efficacy
├── malariasimulation_runs             # Scripts to projecte public health impact and cost-effectiveness 
|   ├── 01_data                          # Data files
|   ├── 02_code                          # Coding scripts
|   |   ├── data_and_libraries.R           # Libraries
|   |   ├── A_parameter_draws.R            # List parameters for each draw
|   |   ├── B_PfPR_EIR_match.R             # Match PfPR to EIR
|   |   ├── C_HPC_runs.R                   # Model runs
|   |   ├── D_HPC_processing.R             # Processing model runs
|   |   ├── Functions                      # Functions to simulate and calculate model outcomes
|   |   |   ├── add_costs.R                  # Helper-functionadd costs
|   |   |   ├── cost_effectiveness.R         # Process CE
|   |   |   ├── deaths_dalys.R               # Helper-function calculate deaths & DALYs
|   |   |   ├── generate_params.R            # Generate dataframe of parameter lists
|   |   |   ├── HPC_processing.R             # Helper-function process model runs
|   |   |   ├── outcome_averted.R            # Helper-function calculate DALYs and cases averted
|   |   |   ├── PRmatch.R                    # Calibrate model via EIR and PfPR
|   |   |   ├── runsim.R                     # Run malariasimulation
|   ├── 03_output                         # Model outputs
├── r21_vacc_antibody_model.Rproj      # R.Studio project file
└── README.md                          # Project overview

```
