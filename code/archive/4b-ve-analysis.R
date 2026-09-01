library(Rcpp)
sourceCpp(here::here("code", "functions", "person-time.cpp"))

# Function to process incidence and VE ----
f.proc.ve <- function(df, super_true_ve) {

  # Person-time for all endpoint groups — single pass over data ----
  all_endpoints <- c(endpoints, sevscore_endpoints, sevgems_endpoints)
  truth_cols    <- c(rep("shigella_diarrhea",           length(endpoints)),
                     rep("sev.score.shigella_diarrhea", length(sevscore_endpoints)),
                     rep("sev.gems.shigella_diarrhea",  length(sevgems_endpoints)))

  person_time_all <- calculate_person_time_all_cpp(df, truth_cols, all_endpoints)

  # Incidence by country and simulation ----
  incidence_by_country_sim <- person_time_all |>
    left_join(df |> select(sim, country_id, child, vax) |> distinct(),
              by = c("sim", "country_id", "child")) |>
    group_by(sim, country_id, vax, endpoint) |>
    summarise(
      risk_true       = mean(true_event),
      risk_obs        = mean(obs_event),
      incidence_true  = sum(true_event)  / sum(person_time_true),
      incidence_obs   = sum(obs_event)   / sum(person_time_obs),
      incidence_bias  = incidence_obs - incidence_true,
      .groups = "drop"
    )

  # CI summary ----
  incidence_by_country_ci <- incidence_by_country_sim |>
    group_by(country_id, vax, endpoint) |>
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
    ) |>
    label_endpoints()

  # Format incidence ----
  true_df <- incidence_by_country_ci %>%
    select(-starts_with("obs_"), -starts_with("bias_")) %>%
    mutate(measure = "True",
           mean    = true_incidence_mean  * 12 * 100,
           lower   = true_incidence_lower * 12 * 100,
           upper   = true_incidence_upper * 12 * 100) %>%
    select(-starts_with("true_"))

  obs_df <- incidence_by_country_ci %>%
    select(-starts_with("true_"), -starts_with("bias_")) %>%
    mutate(measure = "Observed",
           mean    = obs_incidence_mean  * 12 * 100,
           lower   = obs_incidence_lower * 12 * 100,
           upper   = obs_incidence_upper * 12 * 100) %>%
    select(-starts_with("obs_"))

  bias_df <- incidence_by_country_ci %>%
    select(-starts_with("true_"), -starts_with("obs_")) %>%
    mutate(measure = "Bias",
           mean    = bias_incidence_mean  * 12 * 100,
           lower   = bias_incidence_lower * 12 * 100,
           upper   = bias_incidence_upper * 12 * 100) %>%
    select(-starts_with("bias_"))

  inc_df <- true_df %>%
    add_row(obs_df) %>%
    add_row(bias_df) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan")))

  # Format risk ----
  true_risk_df <- incidence_by_country_ci %>%
    select(-starts_with("obs_"), -starts_with("bias_")) %>%
    mutate(measure = "True",
           mean    = true_risk_mean,
           lower   = true_risk_lower,
           upper   = true_risk_upper) %>%
    select(-starts_with("true_"))

  obs_risk_df <- incidence_by_country_ci %>%
    select(-starts_with("true_"), -starts_with("bias_")) %>%
    mutate(measure = "Observed",
           mean    = obs_risk_mean,
           lower   = obs_risk_lower,
           upper   = obs_risk_upper) %>%
    select(-starts_with("obs_"))

  risk_df <- true_risk_df %>%
    add_row(obs_risk_df) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan")))

  # VE ----
  ve_df <- incidence_by_country_sim |>
    select(sim, country_id, vax, endpoint, incidence_true, incidence_obs) |>
    mutate(incidence_true = ifelse(incidence_true == 0, 1e-20, incidence_true),
           incidence_obs  = ifelse(incidence_obs  == 0, 1e-20, incidence_obs)) |>
    pivot_wider(
      names_from  = vax,
      values_from = c(incidence_true, incidence_obs),
      names_glue  = "{.value}_{vax}"
    ) |>
    mutate(
      ve_sim_true = pmax(0, 1 - (incidence_true_1 / incidence_true_0)),
      ve_obs      = pmax(0, 1 - (incidence_obs_1  / incidence_obs_0)),
      ve_sim_bias = ve_obs - ve_sim_true
    ) |>
    group_by(endpoint) |>
    mutate(
      ve_true = ifelse(grepl("sev", endpoint),
                       mean(ve_sim_true),
                       as.numeric(super_true_ve)),
      ve_bias = ve_obs - ve_true
    ) |>
    ungroup() |>
    group_by(country_id, endpoint) |>
    summarise(
      true_ve           = mean(ve_true),
      true_ve_sim_mean  = mean(ve_sim_true),
      true_ve_sim_lower = quantile(ve_sim_true, 0.025),
      true_ve_sim_upper = quantile(ve_sim_true, 0.975),
      obs_ve_mean       = mean(ve_obs),
      obs_ve_lower      = quantile(ve_obs, 0.025),
      obs_ve_upper      = quantile(ve_obs, 0.975),
      bias_ve_sim_mean  = mean(ve_sim_bias),
      bias_ve_sim_lower = quantile(ve_sim_bias, 0.025),
      bias_ve_sim_upper = quantile(ve_sim_bias, 0.975),
      bias_ve_mean      = mean(ve_bias),
      bias_ve_lower     = quantile(ve_bias, 0.025),
      bias_ve_upper     = quantile(ve_bias, 0.975),
      .groups = "drop"
    ) |>
    label_endpoints()

  true_ve_df <- ve_df |>
    select(country_id, endpoint, severity, true_ve) |>
    mutate(measure = "True VE", mean = true_ve, lower = true_ve, upper = true_ve) |>
    select(-true_ve)

  true_ve_sim_df <- ve_df |>
    select(country_id, endpoint, severity, starts_with("true_ve_sim_")) |>
    mutate(measure = "True VE (sim)") |>
    rename(mean = true_ve_sim_mean, lower = true_ve_sim_lower, upper = true_ve_sim_upper)

  obs_ve_df <- ve_df |>
    select(country_id, endpoint, severity, starts_with("obs_ve_")) |>
    mutate(measure = "Observed VE") |>
    rename(mean = obs_ve_mean, lower = obs_ve_lower, upper = obs_ve_upper)

  bias_ve_df <- ve_df |>
    select(country_id, endpoint, severity, bias_ve_mean, bias_ve_lower, bias_ve_upper) |>
    mutate(measure = "VE Bias") |>
    rename(mean = bias_ve_mean, lower = bias_ve_lower, upper = bias_ve_upper)

  bias_ve_sim_df <- ve_df |>
    select(country_id, endpoint, severity, starts_with("bias_ve_sim_")) |>
    mutate(measure = "VE Bias (sim)") |>
    rename(mean = bias_ve_sim_mean, lower = bias_ve_sim_lower, upper = bias_ve_sim_upper)

  ve_out <- true_ve_df |>
    add_row(true_ve_sim_df) |>
    add_row(obs_ve_df) |>
    add_row(bias_ve_df) |>
    add_row(bias_ve_sim_df) |>
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan")))

  return(list(incidence = inc_df, risk = risk_df, ve = ve_out))
}