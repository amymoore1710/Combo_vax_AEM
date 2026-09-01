# libraries ----
pacman::p_load(dplyr,
               tidyr,
               future,
               furrr,
               purrr,
               scales,
               stringr)

options(future.globals.maxSize = 3 * 1024^3)

# Functions (also defines endpoint vectors) ----
source(here::here("code","4a-sens-spec-analysis.R"))
source(here::here("code","4b-ve-analysis.R"))

################################################################################
# Define dataset names and file paths ----
dataset_names <- c(
  "Main", "VE.inf0.0", "VE.inf0.2", "VE.dis0.3", "VE.dis0.5",
  "VE.sev0.0", "VE.sev0.3", "VE.ct0.0", "VE.ct0.2"
)

file_names <- c("ve_0.10_0.40_0.15_0.10", "ve_0.00_0.40_0.15_0.10", "ve_0.20_0.40_0.15_0.10",
                "ve_0.10_0.30_0.15_0.10", "ve_0.10_0.50_0.15_0.10", "ve_0.10_0.40_0.00_0.10",
                "ve_0.10_0.40_0.30_0.10", "ve_0.10_0.40_0.15_0.00", "ve_0.10_0.40_0.15_0.20")

dataset_paths <- paste0("sim data/sim_results_", file_names, ".RDS")

### Step 1: Estimate Sensitivity & Specificity ----
sens_spec_results <- list()

for (i in seq_along(dataset_names)) {
  scenario     <- dataset_names[i]
  dataset_path <- here::here(dataset_paths[i])

  data      <- readRDS(dataset_path)
  ss_result <- f.proc.ss(data) %>% mutate(scenario = scenario)

  sens_spec_results[[scenario]] <- ss_result
  rm(data, ss_result); gc()
}

ss_combined <- bind_rows(sens_spec_results)
saveRDS(ss_combined, here::here("clean data", "sens_spec_results.RDS"))
rm(sens_spec_results, ss_combined)

### Step 2: Estimate Incidence, Risk, and VE ----
inc_results  <- list()
risk_results <- list()
ve_results   <- list()

for (i in seq_along(dataset_names)) {
  scenario     <- dataset_names[i]
  dataset_path <- here::here(dataset_paths[i])
  file_name    <- file_names[i]

  super_true_ve <- sub("^[^_]*_[^_]*_([^_]*).*", "\\1", file_name)

  data      <- readRDS(dataset_path)
  ve_result <- f.proc.ve(data, super_true_ve)
  ve_result <- map(ve_result, ~ mutate(.x, scenario = scenario))

  saveRDS(ve_result, here::here("clean data", "scenario ve results",
                                paste0(scenario, "_v2.rds")))

  inc_results[[scenario]]  <- ve_result$incidence
  risk_results[[scenario]] <- ve_result$risk
  ve_results[[scenario]]   <- ve_result$ve

  rm(data, ve_result); gc()
}

inc_df  <- bind_rows(inc_results)
risk_df <- bind_rows(risk_results)
ve_df   <- bind_rows(ve_results)

saveRDS(inc_df,  here::here("clean data", "inc_results3.RDS"))
saveRDS(risk_df, here::here("clean data", "risk_results3.RDS"))
saveRDS(ve_df,   here::here("clean data", "ve_results3.RDS"))