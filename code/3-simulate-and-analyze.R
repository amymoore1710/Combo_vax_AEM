pacman::p_load(future.apply, 
               Rcpp, 
               dplyr, 
               tidyr, 
               purrr,
               parallel,
               foreach,
               doParallel)

options(future.globals.maxSize = 2 * 1024^3)

# Parameters and analysis functions ----
source(here::here("code", "1-parameters.R"))
source(here::here("code", "functions", "sens-spec-analysis.R"))     # defines endpoints, label_endpoints, sens_spec_fun
source(here::here("code", "functions", "sim-analysis-functions.R"))    # defines within_sim_summarise, f.aggregate_ss, f.aggregate_ve

# Compile C++ once in main process ----
message("Compiling simulation functions...")
sourceCpp(here::here("code", "2-functions.cpp"))
message("Compiling person-time functions...")
sourceCpp(here::here("code", "functions", "person-time.cpp"))
message("Done.\n")

# chunk_size controls how many sims run before workers are recycled ----
# reduce if memory is still tight; increase to reduce worker startup overhead
  # ***Cannot use simualtion "chunks" for fitting CoxPH model so gonna have to reduce sample size
# n_simulations <- 1000
# n_children    <- 20000
# chunk_size    <- 100

n_simulations <- 10
n_children    <- 100
chunk_size    <- 0 


message(sprintf("Settings: %d sims x %d children, chunk size %d.\n",
                n_simulations, n_children, chunk_size))

# ---------------------------------------------------------------
# Smoke test: run ONE simulation at small scale to gut-check the output
# before committing to the full scenario loop. Uses the "Main" scenario's
# VE parameters but a small n_children so it runs in seconds. Inspect
# `test_raw` (row-level simulate_trial output) and the printed summaries
# below to sanity-check things like score/Ct distributions, event rates,
# and that the severity/culture chain looks reasonable.
# ---------------------------------------------------------------
run_smoke_test <- FALSE   # set FALSE to skip and go straight to the full run

if (run_smoke_test) {
  message("Running smoke test (1 simulation, small n_children)...")

  test_params <- f.param(ve_infection_shig = 0.10, ve_disease_shig = 0.40,
                         ve_severity_shig  = 0.60, ve_ct_shig      = 0.10,
                         ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40,
                         ve_severity_ETEC  = 0.60, ve_ct_ETEC      = 0.10)

  test_raw <- simulate_trial(
    n_children = 100, 
    n_months = n_months, 
    countries = countries,
    diarrhea_ir = test_params[["diarrhea_ir"]], 
    shigella_sub = test_params[["shigella_sub"]], 
    ETEC_sub = test_params[["ETEC_sub"]],
    other_sub = other_sub,
    shigella_sev_params = test_params[["shigella_sev_params"]], 
    ETEC_sev_params = test_params[["ETEC_sev_params"]],
    other_sev_params = test_params[["other_sev_params"]],
    shigella_quantity = test_params[["shigella_quantity"]], 
    ETEC_quantity = test_params[["ETEC_quantity"]],
    other_quantity = other_quantity
  )

  message(sprintf("  Smoke test produced %d rows, %d columns.", nrow(test_raw), ncol(test_raw)))

  message("\n--- Column summary ---")
  print(summary(test_raw))

  message("\n--- Event rates ---")
  print(test_raw %>%
    summarise(
      n_rows                = n(),
      pct_shigella_diarrhea = mean(shigella_diarrhea, na.rm = TRUE),
      pct_other_diarrhea    = mean(other_diarrhea, na.rm = TRUE),
      pct_shigella_subclin  = mean(shigella_sub, na.rm = TRUE),
      pct_vax               = mean(vax)
    ))

  message("\n--- Shigella score / gemsdef / culture, among Shigella diarrhea cases ---")
  print(test_raw %>%
    filter(shigella_diarrhea == 1) %>%
    mutate(shigella_ct = 35 - 3.322 * shigella_quantity) %>%
    summarise(
      n             = n(),
      mean_score    = mean(shigella_score, na.rm = TRUE),
      pct_gemsmsd   = mean(shigella_gemsmsd, na.rm = TRUE),
      mean_quantity = mean(shigella_quantity, na.rm = TRUE),
      mean_ct       = mean(shigella_ct, na.rm = TRUE),
      sd_ct         = sd(shigella_ct, na.rm = TRUE),
      pct_culture   = mean(shigella_culture, na.rm = TRUE)
    ))

  message("\n--- Shigella quantity / culture, among subclinical detections ---")
  print(test_raw %>%
    filter(shigella_diarrhea != 1, shigella_sub == 1) %>%
    mutate(shigella_ct = 35 - 3.322 * shigella_quantity) %>%
    summarise(
      n             = n(),
      mean_quantity = mean(shigella_quantity, na.rm = TRUE),
      mean_ct       = mean(shigella_ct, na.rm = TRUE),
      sd_ct         = sd(shigella_ct, na.rm = TRUE),
      pct_culture   = mean(shigella_culture, na.rm = TRUE)
    ))

  message("\n--- Shigella Ct by severity group (converted from simulated quantity) ---")
  print(test_raw %>%
    filter(!is.na(shigella_quantity)) %>%
    mutate(
      shigella_ct = 35 - 3.322 * shigella_quantity,
      group = case_when(
        shigella_diarrhea == 1 & (shigella_score >= 6 | shigella_gemsmsd == 1) ~ "Diarrhea - Severe",
        shigella_diarrhea == 1                                                  ~ "Diarrhea - Mild",
        TRUE                                                                    ~ "Subclinical"
      )
    ) %>%
    group_by(group) %>%
    summarise(
      n             = n(),
      mean_quantity = mean(shigella_quantity), sd_quantity = sd(shigella_quantity),
      mean_ct       = mean(shigella_ct),        sd_ct       = sd(shigella_ct),
      min_ct        = min(shigella_ct),         max_ct      = max(shigella_ct),
      .groups = "drop"
    ))

  message("\n--- Other pathogen detection / quantity / Ct ---")
  print(test_raw %>%
    filter(other_inf == 1) %>%
    mutate(other_ct = 35 - 3.322 * other_quantity) %>%
    summarise(
      n                   = n(),
      pct_other_inf       = mean(other_inf, na.rm = TRUE),
      mean_quantity_other = mean(other_quantity, na.rm = TRUE),
      mean_ct_other       = mean(other_ct, na.rm = TRUE),
      sd_ct_other         = sd(other_ct, na.rm = TRUE),
      min_ct_other        = min(other_ct, na.rm = TRUE),
      max_ct_other        = max(other_ct, na.rm = TRUE)
    ))

  message("\n--- Sens/spec confusion matrix: all Shigella diarrhea vs. Shigella MSD (score>=6), Bangladesh, endpoint 1 ---")
  # Reuses the exact TP/FN/TN/FP logic from within_sim_summarise()/sens_spec_fun()
  # (endpt1 = any Shigella detected among diarrhea; sev1.endpt1 = same, restricted
  # to score>=6 "MSD" diarrhea) so these numbers match what the full pipeline
  # would produce for a single simulation.
  smoke_ss <- within_sim_summarise(test_raw, sim_id = 1)$ss %>%
    filter(country_id == "ALL", endpoint %in% c("endpt1", "sev1.endpt1")) %>%
    mutate(group = ifelse(endpoint == "endpt1", "All Shigella diarrhea", "Shigella MSD (score>=6)"))

  print(smoke_ss %>% select(group, TP, FN, TN, FP, sensitivity, specificity))

  all_row <- smoke_ss %>% filter(group == "All Shigella diarrhea")
  msd_row <- smoke_ss %>% filter(group == "Shigella MSD (score>=6)")

  message(sprintf(
    paste0(
      "\nSentence check: using endpoint 1 (no false negatives) in a single simulation, ",
      "there were %s true positives, %s true negatives, and %s false positives for all Shigella diarrhea ",
      "(specificity = %.1f%%), compared to %s true positives, %s true negatives and %s false positives for ",
      "Shigella MSD with severity score >=6 (specificity = %.1f%%)."
    ),
    format(all_row$TP, big.mark = ","), format(all_row$TN, big.mark = ","), format(all_row$FP, big.mark = ","),
    all_row$specificity * 100,
    format(msd_row$TP, big.mark = ","), format(msd_row$TN, big.mark = ","), format(msd_row$FP, big.mark = ","),
    msd_row$specificity * 100
  ))

  message("\nSmoke test complete — inspect test_raw / printed summaries above before proceeding.\n")

  # Clear the ve_max calibration cache populated by the smoke test's
  # f.param() call above. Without this, the real scenario loop's first
  # encounter with ve_severity = 0.60 (the "Main" scenario) would silently
  # reuse the smoke test's calibration instead of running its own — the
  # cache is meant to guarantee identical ve_max ACROSS scenarios that
  # share a target, not to be seeded by an unrelated test run.
  rm(list = ls(.ve_max_cache), envir = .ve_max_cache)
  message("Cleared ve_max cache populated during smoke test.\n")
}

rm(test_raw)

# ---------------------------------------------------------------
# VE scenario definitions
# ---------------------------------------------------------------

ve_params_list <- list(
  list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
       ve_severity_shig = 0.60, ve_ct_shig = 0.10,
       ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
       ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
  list(ve_infection_shig = 0.00, ve_disease_shig = 0.40, 
       ve_severity_shig = 0.60, ve_ct_shig = 0.10,
       ve_infection_ETEC = 0.00, ve_disease_ETEC = 0.40, 
       ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10)
)

# ve_params_list <- list(
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.00, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.00, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.20, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.20, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.30, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.30, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.50, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.50, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.40, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.40, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.80, ve_ct_shig = 0.10,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.80, ve_ct_ETEC = 0.10),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.00,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.00),
#   list(ve_infection_shig = 0.10, ve_disease_shig = 0.40, 
#        ve_severity_shig = 0.60, ve_ct_shig = 0.20,
#        ve_infection_ETEC = 0.10, ve_disease_ETEC = 0.40, 
#        ve_severity_ETEC = 0.60, ve_ct_ETEC = 0.20)
# )

dataset_names <- c(
  "Main", "VE.inf0.0", "VE.inf0.2", "VE.dis0.3", "VE.dis0.5",
  "VE.sev0.0", "VE.sev0.3", "VE.ct0.0", "VE.ct0.2"
)

n_scenarios  <- length(ve_params_list)
total_start  <- proc.time()

# # Accumulators for cross-scenario outputs
# ss_all   <- list()
# inc_all  <- list()
# risk_all <- list()
# ve_all   <- list()

# ---------------------------------------------------------------
# Main scenario loop
# ---------------------------------------------------------------
for (i in seq_along(ve_params_list)) {

  ve       <- ve_params_list[[i]]
  scenario <- dataset_names[i]
  tag      <- sprintf("VEinf=%.2f / VEdis=%.2f / VEsev=%.2f / VEct=%.2f",
                      ve$ve_infection_shig, ve$ve_disease_shig, ve$ve_severity_shig, ve$ve_ct_shig)
  filename <- sprintf("sVEinf_%.2f_sVEdis=%.2f_sVEsev=%.2f_sVEct=%.2f_eVEinf_%.2f_eVEdis=%.2f_eVEsev=%.2f_eVEct=%.2f",
                      ve$ve_infection_shig, ve$ve_disease_shig, ve$ve_severity_shig, ve$ve_ct_shig,
                      ve$ve_infection_ETEC, ve$ve_disease_ETEC, ve$ve_severity_ETEC, ve$ve_ct_ETEC)

  message(sprintf("[%d/%d] %s", i, n_scenarios, tag))
  scenario_start <- proc.time()

  params <- f.param(ve_infection_shig = ve$ve_infection_shig, ve_disease_shig = ve$ve_disease_shig,
                    ve_severity_shig  = ve$ve_severity_shig,  ve_ct_shig      = ve$ve_ct_shig,
                    ve_infection_ETEC = ve$ve_infection_ETEC, ve_disease_ETEC = ve$ve_disease_ETEC,
                    ve_severity_ETEC  = ve$ve_severity_ETEC,  ve_ct_ETEC      = ve$ve_ct_ETEC)

  diarrhea_ir2          <- params[["diarrhea_ir"]]
  shigella_sub2         <- params[["shigella_sub"]]
  ETEC_sub2             <- params[["ETEC_sub"]]
  shigella_sev_params2  <- params[["shigella_sev_params"]]
  ETEC_sev_params2      <- params[["ETEC_sev_params"]]
  other_sev_params2     <- params[["other_sev_params"]]
  shigella_quantity2    <- params[["shigella_quantity"]]
  ETEC_quantity2        <- params[["ETEC_quantity"]]

  super_true_ve <- as.character(ve$ve_disease_shig)

  # Low-noise analytic true VE (per site + pooled, score and gems MSD
  # separately) computed directly from this scenario's calibrated ve_max
  # and DGP parameters — see compute_true_ve_by_site() in 1-parameters.R.
  # Passed into f.aggregate_ve() to use as the unsuffixed "True VE"/
  # "VE Bias" benchmark instead of the noisier simulated-trial average.
  analytic_true_ve_shig <- compute_true_ve_by_site_shig(
    score_params = shigella_sev_params2$score_params,
    gems_coefs   = shigella_sev_params2$gems_params$coefs,
    diarrhea_ir  = diarrhea_ir2,
    ve_max       = shigella_sev_params2$ve_max,
    ve_disease   = ve$ve_disease_shig
  )
  
  analytic_true_ve_ETEC <- compute_true_ve_by_site_ETEC(
    score_params = ETEC_sev_params2$score_params,
    gems_coefs   = ETEC_sev_params2$gems_params$coefs,
    diarrhea_ir  = diarrhea_ir2,
    ve_max       = ETEC_sev_params2$ve_max,
    ve_disease   = ve$ve_disease_ETEC
  )

  # --- chunked simulation ---
  # Run n_simulations in chunks of chunk_size. Each chunk is completed and
  # its raw summaries discarded before the next chunk starts, so memory
  # never accumulates beyond one chunk at a time.
  # Only the tiny per-sim intermediate rows (~1 row per sim x country x vax
  # x endpoint) are kept across chunks for the final quantile step.
  
  # ***Cannot chunk for CoxPH analysis...
  # chunks   <- split(seq_len(n_simulations),
  #                   ceiling(seq_len(n_simulations) / chunk_size))
  # n_chunks <- length(chunks)
  # 
  # ss_rows <- vector("list", n_simulations)
  # pt_rows <- vector("list", n_simulations)

  # for (ch in seq_along(chunks)) {
  #   sim_ids     <- chunks[[ch]]
  #   chunk_start <- proc.time()
  #   message(sprintf("  Chunk %d/%d (sims %d-%d)...",
  #                   ch, n_chunks, min(sim_ids), max(sim_ids)))
  # 
  #   # Recycle workers at the start of every chunk: kills old worker processes
  #   # and spawns fresh ones so worker-held memory is fully released.
  #   plan(sequential); gc(); plan(multisession)
  # 
  #   chunk_summaries <- future_lapply(
  #     sim_ids,
  #     future.seed = 2024L,
  #     function(x) {
  #       Rcpp::sourceCpp(here::here("code", "2-functions.cpp"),              rebuild = FALSE)
  #       Rcpp::sourceCpp(here::here("code", "functions", "person-time.cpp"), rebuild = FALSE)
  #       source(here::here("code", "functions", "sens-spec-analysis.R"))
  #       source(here::here("code", "functions", "sim-analysis-functions.R"))
  #       raw <- simulate_trial(n_children, n_months, countries,
  #                             diarrhea_ir2, shigella_sub2, other_sub,
  #                             sev_params2, other_sev_params2,
  #                             shigella_quantity2, other_quantity)
  #       within_sim_summarise(raw, sim_id = x)
  #     }
  #   )
  # 
  #   for (j in seq_along(sim_ids)) {
  #     sid            <- sim_ids[j]
  #     ss_rows[[sid]] <- chunk_summaries[[j]]$ss
  #     pt_rows[[sid]] <- chunk_summaries[[j]]$pt
  #   }
  #   # Kill workers immediately after extracting results
  #   rm(chunk_summaries); plan(sequential); gc()
  # 
  #   n_stored <- sum(!sapply(pt_rows, is.null))
  #   message(sprintf("    Done (%.1f sec). Stored %d/%d sim summaries.",
  #                   (proc.time() - chunk_start)[["elapsed"]],
  #                   n_stored, n_simulations))
  # }
  
  
  # Set up Parallel Computing
  detectCores()
  num_cores <- detectCores() - 2
  cl <- makeCluster(num_cores)
  
  simulations <- foreach(i = 1:n_simulations, .combine = rbind) %dopar% {
    simulate_trial(
      n_children = n_children, 
      n_months = n_months, 
      countries = countries,
      diarrhea_ir = diarrhea_ir2, 
      shigella_sub = shigella_sub2, 
      ETEC_sub = ETEC_sub2,
      other_sub = other_sub,
      shigella_sev_params = shigella_sev_params2, 
      ETEC_sev_params = ETEC_sev_params2,
      other_sev_params = other_sev_params2,
      shigella_quantity = shigella_quantity2, 
      ETEC_quantity = ETEC_quantity2,
      other_quantity = other_quantity
    )
  }
  
  simulations$sim <- rep(x = 1:n_simulations, each = n_children*4*12)

  sim_elapsed <- (proc.time() - scenario_start)[["elapsed"]]
  message(sprintf("  All chunks complete (%.1f sec). Aggregating...", sim_elapsed))
  
  simulations <- endpoint_creation(simulations)
  
  # Final aggregation across all sims
  # ss_result <- f.aggregate_ss(ss_rows) %>% mutate(scenario = scenario)
  # ss_all[[scenario]] <- ss_result
  # 
  # ve_result <- f.aggregate_ve(pt_rows, super_true_ve, analytic_true_ve)
  # ve_result <- map(ve_result, ~ mutate(.x, scenario = scenario))
  # 
  # saveRDS(ve_result, here::here("clean data", "scenario ve results",
  #                               paste0(scenario, "_v2.rds")))
  # 
  # inc_all[[scenario]]  <- ve_result$incidence
  # risk_all[[scenario]] <- ve_result$risk
  # ve_all[[scenario]]   <- ve_result$ve
  # 
  # rm(ss_rows, pt_rows, ss_result, ve_result); gc()

  total_elapsed   <- (proc.time() - scenario_start)[["elapsed"]]
  overall_elapsed <- (proc.time() - total_start)[["elapsed"]]
  remaining       <- (overall_elapsed / i) * (n_scenarios - i)

  message(sprintf("  Scenario time: %.1f sec (sim: %.1f, aggregate: %.1f)",
                  total_elapsed, sim_elapsed, total_elapsed - sim_elapsed))
  message(sprintf("  %d/%d complete — est. %.1f min remaining.\n",
                  i, n_scenarios, remaining / 60))
}

# # Save combined outputs ----
saveRDS(simulations, here::here("sim data", paste0(filename, ".RDS")))

# saveRDS(bind_rows(ss_all),   here::here("clean data", "sens_spec_results.RDS"))
# saveRDS(bind_rows(inc_all),  here::here("clean data", "inc_results3.RDS"))
# saveRDS(bind_rows(risk_all), here::here("clean data", "risk_results3.RDS"))
# saveRDS(bind_rows(ve_all),   here::here("clean data", "ve_results3.RDS"))

total_time <- (proc.time() - total_start)[["elapsed"]]
message(sprintf("\nAll %d scenarios complete. Total time: %.1f min.",
                n_scenarios, total_time / 60))
