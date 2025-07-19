Code from: No one-size-fits-all: trait-dependent effects of local plant diversity on pollinators and pollination service in a densifying city
Contact: Merin Reji Chacko (merin.rejichacko@gmail.com)

## Usage

To reproduce the models:

1. Go to Envidat. NOTE: Most scripts in `02_analyses/` will not run without data being added to the "raw_data/RejiChacko_etal_2025_EnviDat" directory.
Simply copy and unzip all the files from here: Reji Chacko, M., Moretti, M., Frey, D. (2025). Comprehensive dataset of pollinator diversity and visitation rates with individual-based traits and pollination success across four urban garden plant species. EnviDat. https://www.doi.org/10.16904/envidat.676.
2. Open `github.Rproj` in RStudio.
3. Sequently run scripts from "scripts/02_analyses/".
4. OPTIONAL: to reproduce the files in the "cleaned_data/" directory, the scripts from "scripts/01_data_prep" can be run. But the output files have already been provided in the "cleaned_data" directory.

## Directory structure
.
├── cleaned_data # these datasets are created by the scripts in the directory: 01_data_preparation
│   ├── pollinators_sum_stat.txt
│   ├── pollinator_abundance_aggregated_garden_phytometer.txt
│   ├── pollinator_SR_aggregated_garden_phytometer.txt
│   ├── sampling_effort_gardens_plants_aggregated.txt
│   └── sampling_effort_gardens_plants_aggregated_per_date.txt
│
├── environments # files in here will be created in the process of running the scripts
│   ├── 1a_1b_environment.RData # example name
│   └── 3a_environment.RData # example name
│
├── raw_data
│   ├── explanatory_variables.txt # Data from: Frey, D. & Moretti, M. A comprehensive dataset on cultivated and spontaneously growing vascular plants in urban gardens. Data in Brief 25, 103982 (2019).
│   ├── plant_floristic_data.txt # Data from: Frey, D. & Moretti, M. A comprehensive dataset on cultivated and spontaneously growing vascular plants in urban gardens. Data in Brief 25, 103982 (2019).
│   └── RejiChacko_etal_2025_EnviDat/   # Download contents from: Reji Chacko, M., Moretti, M. & Frey, D. Comprehensive dataset of pollinator diversity and visitation rates with individual-based traits and pollination success across four urban garden plant species. (2025). doi:10.16904/envidat.676
│
├── results # outputs of the scripts, provided here for ease of use
│   ├── abundance_models_coefficents_all_scales.txt
│   ├── abundance_models_model_fit_all_scales.txt
│   ├── beta_coefficients_all_plants.txt
│   ├── beta_coefficients_all_plants_for_paper.txt
│   ├── beta_coefficients_reproductive_success_all.txt
│   ├── beta_coefficients_reproductive_success_garden_landscape.txt
│   ├── Beta_coefficients_TraitFiltering_ITD.txt
│   └── Beta_coefficients_TraitFiltering_tongue.txt
│
├── scripts
│   ├── 01_data_preparation # this is provided for reproducibility, but the outputs are provided in "cleaned_data"
│   │   ├── 01_summarising_fieldwork_data.R
│   │   ├── 02_summarising_pollinators.R
│   │   ├── 03_calculating_pollinator_abundances.R
│   │   └── 04_calculating_pollinator_species_richness.R
│   │
│   └── 02_analyses # these are the codes needed to reproduce the analyses
│       ├── 1a_abundance_models_figs_S1_S2.R
│       ├── 1b_species_richness_models.R
│       ├── 1c_effect_plots_abundance_fig_2.R
│       ├── 1d_effect_plots_richness_fig_S3.R
│       ├── 1e_coef_plots_abundance_and_richness_figs_3_S4.R
│       ├── 2_trait_models__figs_4_5.R
│       ├── 3a_reproductive_succ_1_data_prep.R
│       ├── 3b_reproductive_succ_2_carrot~abundance.R
│       ├── 3c_reproductive_succ_3_radish_fruitset~abundance.R
│       ├── 3d_reproductive_succ_4_radish_seedset~abundance.R
│       ├── 3e_reproductive_succ_5_onobrychis~abundance.R
│       ├── 3f_reproductive_succ_6_comfreyfruitset~abundance.R
│       ├── 3g_reproductive_succ_7_comfreyseedset~abundance.R
│       ├── 3h_reproductive_succ_8_bayesian_inference~abundance.R
│       ├── 3i_reproductive_succ_9_effect_plots~abundance_fig_S5.R
│       ├── 4a_reproductive_succ_1_carrot~habitat.R
│       ├── 4b_reproductive_succ_2_radish_fruitset~habitat.R
│       ├── 4c_reproductive_succ_3_radish_seedset~habitat.R
│       ├── 4d_reproductive_succ_4_sainfoin~habitat.R
│       ├── 4e_reproductive_succ_5_comfrey_fruitset~habitat.R
│       ├── 4f_reproductive_succ_6_comfrey_seedset~habitat.R
│       ├── 4g_reproductive_succ_7_bayesian_inference~habitat.R
│       └── 4h_reproductive_succ_8_effect_plots~habitat_figure_6.R
│
├── github.Rproj # Open this R Project for ease of use
├── R_requirements.txt # list of R packages needed to run codes
└── README.txt # YOU ARE HERE
