# Endpoint column name vectors (used in simulation, analysis, and workers) ----
endpoints <- c("endpt1", "endpt2.1", "endpt2.2", "endpt3",
               "endpt4.1", "endpt4.2", "endpt.culture")

sevscore_endpoints <- c("sev1.endpt1", "sev1.endpt2.1", "sev1.endpt2.2", "sev1.endpt3",
                        "sev1.endpt4.1", "sev1.endpt4.2", "sev1.endpt.culture")

sevgems_endpoints  <- c("sev2.endpt1", "sev2.endpt2.1", "sev2.endpt2.2", "sev2.endpt3",
                        "sev2.endpt4.1", "sev2.endpt4.2", "sev2.endpt.culture")

# Endpoint factor levels and labels (shared with 4b-ve-analysis.R) ----
endpoint_levels <- c("endpt1", "endpt2.1", "endpt2.2", "endpt3",
                     "endpt4.1", "endpt4.2", "endpt.culture")

endpoint_labels <- c(
  "Endpoint 1:\nAny Shigella",
  "Endpoint 2A:\nShigella Ct<28.8",
  "Endpoint 2B:\nShigella Ct<30.4",
  "Endpoint 3:\nAny Shigella +\nno other pathogen",
  "Endpoint 4A:\nShigella Ct<28.8 +\nno other pathogen Ct<30",
  "Endpoint 4B:\nShigella Ct<30.4 +\nno other pathogen Ct<30",
  "Endpoint Culture:\nCulture positive"
)

# Helper to strip sev prefix and apply consistent factor coding ----
label_endpoints <- function(df) {
  df %>%
    mutate(
      severity = case_when(
        grepl("sev1.", endpoint) ~ "MAL-ED Score >=6",
        grepl("sev2.", endpoint) ~ "GEMS MSD",
        TRUE                     ~ "All diarrhea"
      ),
      endpoint = ifelse(grepl("sev1.", endpoint) | grepl("sev2.", endpoint),
                        substr(endpoint, 6, nchar(endpoint)),
                        endpoint)
    ) %>%
    mutate(endpoint = factor(endpoint,
                             levels = endpoint_levels,
                             labels = endpoint_labels))
}

# Function to calculate sensitivity and specificity for one endpoint ----
sens_spec_fun <- function(df, truth, endpoint) {
  df %>%
    group_by(country_id, sim) %>%
    summarise(
      TP = sum(!!sym(truth) == 1 & !!sym(endpoint) == 1),
      FN = sum(!!sym(truth) == 1 & !!sym(endpoint) == 0),
      TN = sum(!!sym(truth) == 0 & !!sym(endpoint) == 0),
      FP = sum(!!sym(truth) == 0 & !!sym(endpoint) == 1),
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      endpoint    = endpoint
    )
}

# Function to process all sensitivity and specificity calculations ----
f.proc.ss <- function(df) {

  df.diar           <- df %>% filter(diarrhea == 1)
  df.sev.score.diar <- df.diar %>% filter(sev.score.diarrhea == 1)
  df.sev.gems.diar  <- df.diar %>% filter(sev.gems.diarrhea  == 1)

  sens_spec_by_sim <- map(endpoints, function(ep) {
    sens_spec_fun(df.diar, "shigella_diarrhea", ep)
  }) %>% list_rbind()

  sevscore_sens_spec_by_sim <- map(sevscore_endpoints, function(ep) {
    sens_spec_fun(df.sev.score.diar, "sev.score.shigella_diarrhea", ep)
  }) %>% list_rbind()

  sevgems_sens_spec_by_sim <- map(sevgems_endpoints, function(ep) {
    sens_spec_fun(df.sev.gems.diar, "sev.gems.shigella_diarrhea", ep)
  }) %>% list_rbind()

  # Summarise across simulations ----
  metrics_by_country_ci <- sens_spec_by_sim %>%
    add_row(sevscore_sens_spec_by_sim) %>%
    add_row(sevgems_sens_spec_by_sim) %>%
    group_by(country_id, endpoint) %>%
    summarise(
      sens_mean  = mean(sensitivity, na.rm = TRUE),
      sens_lower = quantile(sensitivity, 0.025, na.rm = TRUE),
      sens_upper = quantile(sensitivity, 0.975, na.rm = TRUE),
      spec_mean  = mean(specificity, na.rm = TRUE),
      spec_lower = quantile(specificity, 0.025, na.rm = TRUE),
      spec_upper = quantile(specificity, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    label_endpoints()

  # Format ----
  sens <- metrics_by_country_ci %>%
    select(-starts_with("spec")) %>%
    mutate(measure = "Sensitivity") %>%
    rename(mean = sens_mean, lower = sens_lower, upper = sens_upper)

  spec <- metrics_by_country_ci %>%
    select(-starts_with("sens")) %>%
    mutate(measure = "Specificity") %>%
    rename(mean = spec_mean, lower = spec_lower, upper = spec_upper)

  sens %>%
    add_row(spec) %>%
    mutate(country_id = factor(country_id,
                               levels = c("BG", "IN", "PE", "PK"),
                               labels = c("Bangladesh", "India", "Peru", "Pakistan")))
}