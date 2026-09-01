pacman::p_load(future.apply, Rcpp, dplyr)

# parameters ----
source(here::here("code", "1-parameters.R"))

# compile C++ functions once, before any parallelism ----
message("Compiling C++ functions...")
sourceCpp(here::here("code", "2-functions.cpp"))
cpp_lib <- sourceCpp(here::here("code", "2-functions.cpp"), rebuild = FALSE)$output_file
message("Done.\n")

# set number of children and simulations ----
n_simulations <- 1000
n_children    <- 20000

# set up parallel workers once, before the scenario loop ----
message(sprintf("Setting up parallel workers (%d simulations x %d children)...",
                n_simulations, n_children))
plan(multisession)
message("Done.\n")

# quick single-run smoke test ----
message("Running single test simulation...")
params_test <- f.param(ve_infection = 0.10, ve_disease = 0.40,
                       ve_severity  = 0.15, ve_ct      = 0.10)
one_sim_cpp <- simulate_trial(
  n_children, n_months, countries,
  params_test[["diarrhea_ir"]], params_test[["shigella_sub"]], other_sub,
  params_test[["sev_params"]], params_test[["other_sev_params"]],
  params_test[["shigella_ct"]], other_ct
)
message(sprintf("Test simulation OK — %d rows, %d columns.\n",
                nrow(one_sim_cpp), ncol(one_sim_cpp)))

# VE parameter sets ----
ve_params_list <- list(
  list(ve_infection = 0.10, ve_disease = 0.40, ve_severity = 0.15, ve_ct = 0.10),
  list(ve_infection = 0.00, ve_disease = 0.40, ve_severity = 0.15, ve_ct = 0.10),
  list(ve_infection = 0.20, ve_disease = 0.40, ve_severity = 0.15, ve_ct = 0.10),
  list(ve_infection = 0.10, ve_disease = 0.30, ve_severity = 0.15, ve_ct = 0.10),
  list(ve_infection = 0.10, ve_disease = 0.50, ve_severity = 0.15, ve_ct = 0.10),
  list(ve_infection = 0.10, ve_disease = 0.40, ve_severity = 0.00, ve_ct = 0.10),
  list(ve_infection = 0.10, ve_disease = 0.40, ve_severity = 0.30, ve_ct = 0.10),
  list(ve_infection = 0.10, ve_disease = 0.40, ve_severity = 0.15, ve_ct = 0.00),
  list(ve_infection = 0.10, ve_disease = 0.40, ve_severity = 0.15, ve_ct = 0.20)
)

n_scenarios   <- length(ve_params_list)
total_start   <- proc.time()

# scenario loop ----
for (i in seq_along(ve_params_list)) {

  ve  <- ve_params_list[[i]]
  tag <- sprintf("VEinf=%.2f / VEdis=%.2f / VEsev=%.2f / VEct=%.2f",
                 ve$ve_infection, ve$ve_disease, ve$ve_severity, ve$ve_ct)

  message(sprintf("[%d/%d] Starting scenario: %s", i, n_scenarios, tag))
  scenario_start <- proc.time()

  # build parameters for this scenario
  params <- f.param(
    ve_infection = ve$ve_infection,
    ve_disease   = ve$ve_disease,
    ve_severity  = ve$ve_severity,
    ve_ct        = ve$ve_ct
  )
  diarrhea_ir2      <- params[["diarrhea_ir"]]
  shigella_sub2     <- params[["shigella_sub"]]
  sev_params2       <- params[["sev_params"]]
  other_sev_params2 <- params[["other_sev_params"]]
  shigella_ct2      <- params[["shigella_ct"]]

  # run simulations in parallel
  message(sprintf("  Running %d simulations in parallel...", n_simulations))

  simulations <- future_lapply(
    seq_len(n_simulations),
    future.seed = 2024L,
    function(x) {
      Rcpp::sourceCpp(file = here::here("code", "2-functions.cpp"), rebuild = FALSE)
      simulate_trial(n_children, n_months, countries,
                     diarrhea_ir2, shigella_sub2, other_sub,
                     sev_params2, other_sev_params2,
                     shigella_ct2, other_ct)
    }
  )

  sim_elapsed <- (proc.time() - scenario_start)[["elapsed"]]
  message(sprintf("  Simulations complete (%.1f sec). Post-processing...", sim_elapsed))

  combined_df <- bind_rows(simulations, .id = "sim")

  df <- combined_df %>%
    mutate(
      sim              = as.numeric(sim),
      diarrhea         = ifelse(shigella_diarrhea == 1 | other_diarrhea == 1, 1, 0),
      shigella_ct      = ifelse(is.na(shigella_ct), 35, shigella_ct),
      other_ct         = ifelse(is.na(other_ct), 35, other_ct),
      shigella_score   = ifelse(is.na(shigella_score), 0, shigella_score),
      other_score      = ifelse(is.na(other_score), 0, other_score),
      shigella_gemsmsd = ifelse(is.na(shigella_gemsmsd), 0, shigella_gemsmsd),
      other_gemsmsd    = ifelse(is.na(other_gemsmsd), 0, other_gemsmsd),
      shigella_culture = ifelse(is.na(shigella_culture), 0, shigella_culture),
      # severity flags
      sev.score.diarrhea          = ifelse(diarrhea == 1 &
                                             (shigella_score >= 6 | other_score >= 6), 1, 0),
      sev.gems.diarrhea           = ifelse(diarrhea == 1 &
                                             (shigella_gemsmsd == 1 | other_gemsmsd == 1), 1, 0),
      sev.score.shigella_diarrhea = ifelse(shigella_diarrhea == 1 & shigella_score >= 6, 1, 0),
      sev.gems.shigella_diarrhea  = ifelse(shigella_diarrhea == 1 & shigella_gemsmsd == 1, 1, 0)
    ) %>%
    mutate(
      # PCR-based endpoints
      endpt1        = ifelse(diarrhea == 1 & shigella_ct < 35, 1, 0),
      endpt2.1      = ifelse(diarrhea == 1 & shigella_ct < 28.8, 1, 0),
      endpt2.2      = ifelse(diarrhea == 1 & shigella_ct < 30.4, 1, 0),
      endpt3        = ifelse(diarrhea == 1 & shigella_ct < 35 & other_inf == 0, 1, 0),
      endpt4.1      = ifelse(diarrhea == 1 & shigella_ct < 28.8 & other_ct >= 30, 1, 0),
      endpt4.2      = ifelse(diarrhea == 1 & shigella_ct < 30.4 & other_ct >= 30, 1, 0),
      # Culture endpoint
      endpt.culture = ifelse(diarrhea == 1 & shigella_culture == 1, 1, 0),
      # Severity x PCR
      sev1.endpt1        = ifelse(sev.score.diarrhea == 1 & shigella_ct < 35, 1, 0),
      sev1.endpt2.1      = ifelse(sev.score.diarrhea == 1 & shigella_ct < 28.8, 1, 0),
      sev1.endpt2.2      = ifelse(sev.score.diarrhea == 1 & shigella_ct < 30.4, 1, 0),
      sev1.endpt3        = ifelse(sev.score.diarrhea == 1 & shigella_ct < 35 & other_inf == 0, 1, 0),
      sev1.endpt4.1      = ifelse(sev.score.diarrhea == 1 & shigella_ct < 28.8 & other_ct >= 30, 1, 0),
      sev1.endpt4.2      = ifelse(sev.score.diarrhea == 1 & shigella_ct < 30.4 & other_ct >= 30, 1, 0),
      sev2.endpt1        = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 35, 1, 0),
      sev2.endpt2.1      = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 28.8, 1, 0),
      sev2.endpt2.2      = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 30.4, 1, 0),
      sev2.endpt3        = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 35 & other_inf == 0, 1, 0),
      sev2.endpt4.1      = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 28.8 & other_ct >= 30, 1, 0),
      sev2.endpt4.2      = ifelse(sev.gems.diarrhea == 1 & shigella_ct < 30.4 & other_ct >= 30, 1, 0),
      # Severity x culture
      sev1.endpt.culture = ifelse(sev.score.diarrhea == 1 & shigella_culture == 1, 1, 0),
      sev2.endpt.culture = ifelse(sev.gems.diarrhea  == 1 & shigella_culture == 1, 1, 0)
    )

  filename <- sprintf("sim_results_ve_%.2f_%.2f_%.2f_%.2f.RDS",
                      ve$ve_infection, ve$ve_disease,
                      ve$ve_severity,  ve$ve_ct)

  saveRDS(df, here::here("sim data", filename))

  total_elapsed    <- (proc.time() - scenario_start)[["elapsed"]]
  overall_elapsed  <- (proc.time() - total_start)[["elapsed"]]
  remaining        <- (overall_elapsed / i) * (n_scenarios - i)

  message(sprintf("  Saved: %s", filename))
  message(sprintf("  Scenario time:  %.1f sec (sim: %.1f sec, post-process: %.1f sec)",
                  total_elapsed, sim_elapsed, total_elapsed - sim_elapsed))
  message(sprintf("  Overall: %d/%d scenarios complete — est. %.1f min remaining.\n",
                  i, n_scenarios, remaining / 60))
}

total_time <- (proc.time() - total_start)[["elapsed"]]
message(sprintf("All %d scenarios complete. Total time: %.1f min.",
                n_scenarios, total_time / 60))