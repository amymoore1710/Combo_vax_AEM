# ---------------------------------------------------------------
# sim-analysis-functions.R
#
# Functions for summarising a single simulation run and aggregating
# per-sim summaries into final sensitivity/specificity, incidence,
# risk, and VE results across simulations.
#
# Sourced by 3-simulate-and-analyze.R, both in the main process and
# inside each parallel worker. Depends on 4a-sens-spec-analysis.R
# being sourced first (for endpoints, sevscore_endpoints,
# sevgems_endpoints, label_endpoints(), sens_spec_fun()) and on
# person-time.cpp being compiled (for calculate_person_time_all_cpp()).
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# within_sim_summarise()
#
# Called inside each parallel worker. Takes the raw simulate_trial
# output for ONE simulation run and returns two small summaries:
#
#   $ss  — TP/FN/TN/FP per (country_id, endpoint)   [~4 countries x 21 endpoints]
#   $pt  — person-time per (country_id, child, endpoint, vax)
#          [~n_children x 21 endpoints x 4 countries]
#
# The raw row-level data is never returned to the main process.
# ---------------------------------------------------------------
within_sim_summarise <- function(raw, sim_id) {

  # derived here so they're always available regardless of calling environment
  all_endpoints <- c(endpoints, sevscore_endpoints, sevgems_endpoints)
  truth_cols    <- c(rep("shigella_diarrhea",           length(endpoints)),
                     rep("sev.score.shigella_diarrhea", length(sevscore_endpoints)),
                     rep("sev.gems.shigella_diarrhea",  length(sevgems_endpoints)))

  # --- derive endpoint columns ---
  df <- raw %>%
    mutate(
      sim              = sim_id,
      diarrhea         = ifelse(shigella_diarrhea == 1 | ETEC_diarrhea == 1 | other_diarrhea == 1, 1, 0),
      # convert simulated pathogen quantity to Ct for endpoint thresholds
      # (quantity is simulated directly to avoid Ct truncation — see
      # quantity_to_ct() in 2-functions.cpp/2-functions.R for rationale)
      shigella_ct      = ifelse(is.na(shigella_quantity), 35, 35 - 3.322 * shigella_quantity),
      ETEC_ct          = ifelse(is.na(ETEC_quantity),     35, 35 - 3.322 * ETEC_quantity),
      other_ct         = ifelse(is.na(other_quantity),    35, 35 - 3.322 * other_quantity),
      shigella_score   = ifelse(is.na(shigella_score),    0, shigella_score),
      ETEC_score       = ifelse(is.na(ETEC_score),        0, ETEC_score),
      other_score      = ifelse(is.na(other_score),       0, other_score),
      shigella_gemsmsd = ifelse(is.na(shigella_gemsmsd),  0, shigella_gemsmsd),
      ETEC_gemsmsd     = ifelse(is.na(ETEC_gemsmsd),      0, ETEC_gemsmsd),
      other_gemsmsd    = ifelse(is.na(other_gemsmsd),     0, other_gemsmsd),
      shigella_culture = ifelse(is.na(shigella_culture),  0, shigella_culture),
      ETEC_culture     = ifelse(is.na(ETEC_culture),      0, ETEC_culture),
      sev.score.diarrhea          = ifelse(diarrhea == 1 &
                                             (shigella_score >= 6 | other_score >= 6), 1, 0),
      sev.gems.diarrhea           = ifelse(diarrhea == 1 &
                                             (shigella_gemsmsd == 1 | other_gemsmsd == 1), 1, 0),
      sev.score.shigella_diarrhea = ifelse(shigella_diarrhea == 1 & shigella_score >= 6, 1, 0),
      sev.gems.shigella_diarrhea  = ifelse(shigella_diarrhea == 1 & shigella_gemsmsd == 1, 1, 0),
      sev.score.ETEC_diarrhea     = ifelse(ETEC_diarrhea == 1 & ETEC_score >= 6, 1, 0),
      sev.gems.ETEC_diarrhea      = ifelse(ETEC_diarrhea == 1 & ETEC_gemsmsd == 1, 1, 0),
      endpt1             = ifelse(diarrhea == 1 & shigella_ct < 35, 1, 0),
      endpt2.1           = ifelse(diarrhea == 1 & shigella_ct < 28.8, 1, 0),
      endpt2.2           = ifelse(diarrhea == 1 & shigella_ct < 30.4, 1, 0),
      endpt3             = ifelse(diarrhea == 1 & shigella_ct < 35 & other_inf == 0, 1, 0),
      endpt4.1           = ifelse(diarrhea == 1 & shigella_ct < 28.8 & other_ct >= 30, 1, 0),
      endpt4.2           = ifelse(diarrhea == 1 & shigella_ct < 30.4 & other_ct >= 30, 1, 0),
      endpt.culture      = ifelse(diarrhea == 1 & shigella_culture == 1, 1, 0),
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
      sev1.endpt.culture = ifelse(sev.score.diarrhea == 1 & shigella_culture == 1, 1, 0),
      sev2.endpt.culture = ifelse(sev.gems.diarrhea  == 1 & shigella_culture == 1, 1, 0)
    )

  # --- sens/spec summary: TP/FN/TN/FP per (country_id, endpoint) ---
  df_diar           <- df %>% filter(diarrhea == 1)
  df_sev_score_diar <- df_diar %>% filter(sev.score.diarrhea == 1)
  df_sev_gems_diar  <- df_diar %>% filter(sev.gems.diarrhea  == 1)

  ss <- bind_rows(
    map(endpoints,          ~ sens_spec_fun(df_diar,           "shigella_diarrhea",           .x)),
    map(sevscore_endpoints, ~ sens_spec_fun(df_sev_score_diar, "sev.score.shigella_diarrhea", .x)),
    map(sevgems_endpoints,  ~ sens_spec_fun(df_sev_gems_diar,  "sev.gems.shigella_diarrhea",  .x))
  ) %>%
    mutate(sim = sim_id) %>%
    select(sim, country_id, endpoint, TP, FN, TN, FP, sensitivity, specificity)

  # pooled ("All Sites") counts: sum TP/FN/TN/FP across countries within this
  # sim before recomputing sensitivity/specificity, rather than averaging
  # site-level rates, so sites with more events contribute proportionally
  # more to the pooled estimate.
  ss_pooled <- ss %>%
    group_by(sim, endpoint) %>%
    summarise(
      TP = sum(TP), FN = sum(FN), TN = sum(TN), FP = sum(FP),
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      country_id  = "ALL"
    ) %>%
    select(sim, country_id, endpoint, TP, FN, TN, FP, sensitivity, specificity)

  ss <- bind_rows(ss, ss_pooled)

  # --- per-sim incidence summary: aggregate to country x vax x endpoint ---
  # calculate_person_time_all_cpp returns one row per child x country x endpoint.
  # We immediately reduce to 168 rows (4 countries x 2 vax x 21 endpoints)
  # so that only this tiny summary is returned to the main process.
  pt_by_country <- calculate_person_time_all_cpp(df, truth_cols, all_endpoints) %>%
    left_join(df %>% select(sim, country_id, child, vax) %>% distinct(),
              by = c("sim", "country_id", "child"))

  pt <- pt_by_country %>%
    group_by(sim, country_id, vax, endpoint) %>%
    summarise(
      risk_true      = mean(true_event),
      risk_obs       = mean(obs_event),
      incidence_true = sum(true_event)  / sum(person_time_true),
      incidence_obs  = sum(obs_event)   / sum(person_time_obs),
      incidence_bias = incidence_obs - incidence_true,
      .groups = "drop"
    )

  # pooled ("All Sites") incidence/risk: sum events and person-time across
  # countries within this sim before taking ratios, mirroring the ss
  # pooling above, so sites contribute in proportion to their event/
  # person-time totals rather than being averaged as equally-weighted rates.
  pt_pooled <- pt_by_country %>%
    group_by(sim, vax, endpoint) %>%
    summarise(
      risk_true      = mean(true_event),
      risk_obs       = mean(obs_event),
      incidence_true = sum(true_event)  / sum(person_time_true),
      incidence_obs  = sum(obs_event)   / sum(person_time_obs),
      incidence_bias = incidence_obs - incidence_true,
      .groups = "drop"
    ) %>%
    mutate(country_id = "ALL")

  pt <- bind_rows(pt, pt_pooled)

  list(ss = ss, pt = pt)
}

# ---------------------------------------------------------------
# f.aggregate_ss()  — final sens/spec CI from list of per-sim ss summaries
# ---------------------------------------------------------------
f.aggregate_ss <- function(ss_list) {
  bind_rows(ss_list) %>%
    group_by(country_id, endpoint) %>%
    summarise(
      sens_mean  = mean(sensitivity, na.rm = TRUE),
      sens_lower = quantile(sensitivity, 0.025, na.rm = TRUE),
      sens_upper = quantile(sensitivity, 0.975, na.rm = TRUE),
      spec_mean  = mean(specificity,  na.rm = TRUE),
      spec_lower = quantile(specificity, 0.025, na.rm = TRUE),
      spec_upper = quantile(specificity, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    label_endpoints() %>%
    pivot_longer(c(starts_with("sens_"), starts_with("spec_")),
                 names_to  = c("measure_raw", ".value"),
                 names_sep = "_") %>%
    mutate(measure = ifelse(measure_raw == "sens", "Sensitivity", "Specificity")) %>%
    select(-measure_raw) %>%
    rename(mean = mean, lower = lower, upper = upper) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK", "ALL"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan", "All Sites")))
}

# ---------------------------------------------------------------
# f.aggregate_ve()  — final incidence/risk/VE from list of per-sim pt summaries
# (mirrors f.proc.ve but accepts already-computed person-time rows)
# ---------------------------------------------------------------
f.aggregate_ve <- function(pt_list, super_true_ve, analytic_true_ve) {

  # pt_list entries are already aggregated to (sim x country x vax x endpoint)
  # inside within_sim_summarise — just bind them
  incidence_by_country_sim <- bind_rows(pt_list)

  # CI summary
  incidence_by_country_ci <- incidence_by_country_sim %>%
    group_by(country_id, vax, endpoint) %>%
    summarise(
      true_risk_mean       = mean(risk_true),
      true_risk_lower      = quantile(risk_true, 0.025),
      true_risk_upper      = quantile(risk_true, 0.975),
      obs_risk_mean        = mean(risk_obs),
      obs_risk_lower       = quantile(risk_obs, 0.025),
      obs_risk_upper       = quantile(risk_obs, 0.975),
      true_incidence_mean  = mean(incidence_true),
      true_incidence_lower = quantile(incidence_true, 0.025),
      true_incidence_upper = quantile(incidence_true, 0.975),
      obs_incidence_mean   = mean(incidence_obs),
      obs_incidence_lower  = quantile(incidence_obs, 0.025),
      obs_incidence_upper  = quantile(incidence_obs, 0.975),
      bias_incidence_mean  = mean(incidence_bias),
      bias_incidence_lower = quantile(incidence_bias, 0.025),
      bias_incidence_upper = quantile(incidence_bias, 0.975),
      .groups = "drop"
    ) %>%
    label_endpoints()

  # Format incidence
  inc_df <- bind_rows(
    incidence_by_country_ci %>%
      select(-starts_with("obs_"), -starts_with("bias_")) %>%
      mutate(measure = "True",
             mean  = true_incidence_mean  * 12 * 100,
             lower = true_incidence_lower * 12 * 100,
             upper = true_incidence_upper * 12 * 100) %>%
      select(-starts_with("true_")),
    incidence_by_country_ci %>%
      select(-starts_with("true_"), -starts_with("bias_")) %>%
      mutate(measure = "Observed",
             mean  = obs_incidence_mean  * 12 * 100,
             lower = obs_incidence_lower * 12 * 100,
             upper = obs_incidence_upper * 12 * 100) %>%
      select(-starts_with("obs_")),
    incidence_by_country_ci %>%
      select(-starts_with("true_"), -starts_with("obs_")) %>%
      mutate(measure = "Bias",
             mean  = bias_incidence_mean  * 12 * 100,
             lower = bias_incidence_lower * 12 * 100,
             upper = bias_incidence_upper * 12 * 100) %>%
      select(-starts_with("bias_"))
  ) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK", "ALL"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan", "All Sites")))

  # Format risk
  risk_df <- bind_rows(
    incidence_by_country_ci %>%
      select(-starts_with("obs_"), -starts_with("bias_")) %>%
      mutate(measure = "True",
             mean = true_risk_mean, lower = true_risk_lower, upper = true_risk_upper) %>%
      select(-starts_with("true_")),
    incidence_by_country_ci %>%
      select(-starts_with("true_"), -starts_with("bias_")) %>%
      mutate(measure = "Observed",
             mean = obs_risk_mean, lower = obs_risk_lower, upper = obs_risk_upper) %>%
      select(-starts_with("obs_"))
  ) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK", "ALL"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan", "All Sites")))

  # VE
  ve_df <- incidence_by_country_sim %>%
    select(sim, country_id, vax, endpoint, incidence_true, incidence_obs) %>%
    mutate(incidence_true = ifelse(incidence_true == 0, 1e-20, incidence_true),
           incidence_obs  = ifelse(incidence_obs  == 0, 1e-20, incidence_obs)) %>%
    pivot_wider(names_from  = vax,
                values_from = c(incidence_true, incidence_obs),
                names_glue  = "{.value}_{vax}") %>%
    mutate(
      ve_sim_true = pmax(0, 1 - (incidence_true_1 / incidence_true_0)),
      ve_obs      = pmax(0, 1 - (incidence_obs_1  / incidence_obs_0)),
      # "(sim)" bias: PAIRED per sim — this trial's observed VE vs. this
      # SAME trial's own true VE. Reflects measurement/misclassification
      # bias only, netting out between-sim sampling noise in ve_sim_true
      # itself (since both sides of the subtraction come from the same
      # replicate).
      ve_sim_bias = ve_obs - ve_sim_true
    ) %>%
    # unsuffixed "true VE": looked up from analytic_true_ve (see
    # compute_true_ve_by_site() in 1-parameters.R) rather than computed from
    # simulated trial data — a low-noise, per-site value (score/gems MSD)
    # or the exact-by-construction super_true_ve (All diarrhea), instead of
    # the previous simple unweighted average of the four sites' own noisy
    # simulated true VEs. country_id == "ALL" gets analytic_true_ve's own
    # IR_shigella-weighted pooled row, matching the same "ALL" convention
    # used throughout this function. For "All diarrhea" endpoints this is
    # moot either way: ve_disease is applied as the same multiplicative
    # factor to every site's IR_shigella (vax_IR_shigella = IR_shigella *
    # (1 - ve_disease); see f.param() in 1-parameters.R), so the true VE
    # against any Shigella diarrhea is identical across sites by
    # construction and super_true_ve already IS that value everywhere. For
    # MSD endpoints (score and gems), ve_max is a single scalar too, but
    # works through a nonlinear/discrete mechanism, so the realized true VE
    # genuinely differs by site (see calibrate_ve_max() comments in
    # 1-parameters.R) — analytic_true_ve captures that site-specific value
    # directly from the DGP parameters instead of from noisy simulated
    # event counts. Endpoint family (score vs. gems vs. plain) is
    # determined by the "sev1."/"sev2."/no-prefix naming convention set up
    # in within_sim_summarise() above.
    left_join(analytic_true_ve, by = "country_id") %>%
    mutate(
      ve_true = case_when(
        grepl("^sev1\\.", endpoint) ~ true_ve_score,
        grepl("^sev2\\.", endpoint) ~ true_ve_gems,
        TRUE                        ~ as.numeric(super_true_ve)
      ),
      ve_bias = ve_obs - ve_true
    ) %>%
    select(-true_ve_score, -true_ve_gems) %>%
    # "(site-specific)": each sim's observed VE compared to a FIXED
    # benchmark — the average true VE for THIS site across all its sims —
    # rather than paired per sim as in "(sim)" above. Point estimate ends
    # up numerically equal to "(sim)"'s (mean(A-B) = mean(A) - mean(B)
    # either way), but the bias distribution here also carries the
    # between-sim sampling noise in ve_sim_true that the paired "(sim)"
    # version nets out, since every sim in this site is compared against
    # the same fixed average rather than its own trial's true VE.
    group_by(country_id, endpoint) %>%
    mutate(
      ve_true_site = mean(ve_sim_true),
      ve_bias_site = ve_obs - ve_true_site
    ) %>%
    ungroup() %>%
    group_by(country_id, endpoint) %>%
    summarise(
      true_ve            = mean(ve_true),
      true_ve_sim_mean   = mean(ve_sim_true),
      true_ve_sim_lower  = quantile(ve_sim_true, 0.025),
      true_ve_sim_upper  = quantile(ve_sim_true, 0.975),
      true_ve_site       = mean(ve_true_site),
      obs_ve_mean        = mean(ve_obs),
      obs_ve_lower       = quantile(ve_obs, 0.025),
      obs_ve_upper       = quantile(ve_obs, 0.975),
      bias_ve_mean       = mean(ve_bias),
      bias_ve_lower      = quantile(ve_bias, 0.025),
      bias_ve_upper      = quantile(ve_bias, 0.975),
      bias_ve_sim_mean   = mean(ve_sim_bias),
      bias_ve_sim_lower  = quantile(ve_sim_bias, 0.025),
      bias_ve_sim_upper  = quantile(ve_sim_bias, 0.975),
      bias_ve_site_mean  = mean(ve_bias_site),
      bias_ve_site_lower = quantile(ve_bias_site, 0.025),
      bias_ve_site_upper = quantile(ve_bias_site, 0.975),
      .groups = "drop"
    ) %>%
    label_endpoints()

  ve_out <- bind_rows(
    ve_df %>% select(country_id, endpoint, severity, true_ve) %>%
      mutate(measure = "True VE", mean = true_ve, lower = true_ve, upper = true_ve) %>%
      select(-true_ve),
    ve_df %>% select(country_id, endpoint, severity, starts_with("true_ve_sim_")) %>%
      mutate(measure = "True VE (sim)") %>%
      rename(mean = true_ve_sim_mean, lower = true_ve_sim_lower, upper = true_ve_sim_upper),
    ve_df %>% select(country_id, endpoint, severity, true_ve_site) %>%
      mutate(measure = "True VE (site-specific)", mean = true_ve_site, lower = true_ve_site, upper = true_ve_site) %>%
      select(-true_ve_site),
    ve_df %>% select(country_id, endpoint, severity, starts_with("obs_ve_")) %>%
      mutate(measure = "Observed VE") %>%
      rename(mean = obs_ve_mean, lower = obs_ve_lower, upper = obs_ve_upper),
    ve_df %>% select(country_id, endpoint, severity, bias_ve_mean, bias_ve_lower, bias_ve_upper) %>%
      mutate(measure = "VE Bias") %>%
      rename(mean = bias_ve_mean, lower = bias_ve_lower, upper = bias_ve_upper),
    ve_df %>% select(country_id, endpoint, severity, starts_with("bias_ve_sim_")) %>%
      mutate(measure = "VE Bias (sim)") %>%
      rename(mean = bias_ve_sim_mean, lower = bias_ve_sim_lower, upper = bias_ve_sim_upper),
    ve_df %>% select(country_id, endpoint, severity, starts_with("bias_ve_site_")) %>%
      mutate(measure = "VE Bias (site-specific)") %>%
      rename(mean = bias_ve_site_mean, lower = bias_ve_site_lower, upper = bias_ve_site_upper)
  ) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK", "ALL"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan", "All Sites")))

  list(incidence = inc_df, risk = risk_df, ve = ve_out)
}

endpoint_creation <- function(raw) {
  # --- derive endpoint columns ---
  df <- raw %>%
    mutate(
      sim              = sim,
      diarrhea         = ifelse(shigella_diarrhea == 1 | ETEC_diarrhea == 1 | other_diarrhea == 1, 1, 0),
      # convert simulated pathogen quantity to Ct for endpoint thresholds
      # (quantity is simulated directly to avoid Ct truncation — see
      # quantity_to_ct() in 2-functions.cpp/2-functions.R for rationale)
      shigella_ct      = ifelse(is.na(shigella_quantity), 35, 35 - 3.322 * shigella_quantity),
      ETEC_ct          = ifelse(is.na(ETEC_quantity),     35, 35 - 3.322 * ETEC_quantity),
      other_ct         = ifelse(is.na(other_quantity),    35, 35 - 3.322 * other_quantity),
      shigella_score   = ifelse(is.na(shigella_score),    0, shigella_score),
      ETEC_score       = ifelse(is.na(ETEC_score),        0, ETEC_score),
      other_score      = ifelse(is.na(other_score),       0, other_score),
      shigella_gemsmsd = ifelse(is.na(shigella_gemsmsd),  0, shigella_gemsmsd),
      ETEC_gemsmsd     = ifelse(is.na(ETEC_gemsmsd),      0, ETEC_gemsmsd),
      other_gemsmsd    = ifelse(is.na(other_gemsmsd),     0, other_gemsmsd),
      shigella_culture = ifelse(is.na(shigella_culture),  0, shigella_culture),
      ETEC_culture     = ifelse(is.na(ETEC_culture),      0, ETEC_culture),
      sev.score.diarrhea          = ifelse(diarrhea == 1 &
                                             (shigella_score >= 6 | ETEC_score >= 6 | other_score >= 6), 1, 0),
      sev.gems.diarrhea           = ifelse(diarrhea == 1 &
                                             (shigella_gemsmsd == 1 | ETEC_gemsmsd == 1| other_gemsmsd == 1), 1, 0),
      sev.score.shigella_diarrhea = ifelse(shigella_diarrhea == 1 & shigella_score >= 6, 1, 0),
      sev.gems.shigella_diarrhea  = ifelse(shigella_diarrhea == 1 & shigella_gemsmsd == 1, 1, 0),
      sev.score.ETEC_diarrhea     = ifelse(ETEC_diarrhea == 1 & ETEC_score >= 6, 1, 0),
      sev.gems.ETEC_diarrhea      = ifelse(ETEC_diarrhea == 1 & ETEC_gemsmsd == 1, 1, 0),
      shigella_afe_beta_param = case_when(country_id == "BG" ~ shigella_afe_beta[1],
                                          country_id == "IN" ~ shigella_afe_beta[2],
                                          country_id == "PE" ~ shigella_afe_beta[3],
                                          country_id == "PK" ~ shigella_afe_beta[4]),
      ETEC_afe_beta_param = case_when(country_id == "BG" ~ ETEC_afe_beta[1],
                                      country_id == "IN" ~ ETEC_afe_beta[2],
                                      country_id == "PE" ~ ETEC_afe_beta[3],
                                      country_id == "PK" ~ ETEC_afe_beta[4]),
      other_afe_beta_param = case_when(country_id == "BG" ~ other_afe_beta[1],
                                       country_id == "IN" ~ other_afe_beta[2],
                                       country_id == "PE" ~ other_afe_beta[3],
                                       country_id == "PK" ~ other_afe_beta[4]),
      shigella_afe = 1 - exp(shigella_afe_beta_param * (35 - shigella_ct)),
      ETEC_afe = 1 - exp(ETEC_afe_beta_param * (35 - ETEC_ct)),
      other_afe = 1 - exp(other_afe_beta_param * (35 - other_ct))
    )
  
  return(df)
}
  
  
  
  
  