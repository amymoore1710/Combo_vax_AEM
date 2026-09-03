pacman::p_load(dplyr, tidyr)

# Trial structure parameters ----
n_months   <- 12
age_groups <- c("12-17 months", "18-24 months")
countries  <- c("BG", "PE", "PK", "IN")

# MAL-ED parameters ----
diarrhea_ir  <- readRDS(here::here("sim param data", "diarrhea_incidence.RDS"))
shigella_sub <- readRDS(here::here("sim param data", "shigella_subclinical_prob.RDS"))
ETEC_sub <- readRDS(here::here("sim param data", "ETEC_subclinical_prob.RDS"))
other_sub    <- readRDS(here::here("sim param data", "other_subclinical_prob.RDS"))
sev_params       <- readRDS(here::here("sim param data", "shigella_severity.RDS"))
ETEC_sev_params       <- readRDS(here::here("sim param data", "ETEC_severity.RDS"))
other_sev_params <- readRDS(here::here("sim param data", "other_severity.RDS"))
shigella_quantity <- readRDS(here::here("sim param data", "shigella_quantity.RDS"))
ETEC_quantity <- readRDS(here::here("sim param data", "ETEC_quantity.RDS"))
other_quantity <- readRDS(here::here("sim param data", "other_quantity.RDS"))

# incidence rates ----
ir.tbl <- diarrhea_ir %>%
  pivot_longer(cols = c("IR_other", "IR_shigella", "IR_ETEC"),
               names_to = "parm", values_to = "val") %>%
  mutate(parm = ifelse(parm == "IR_shigella",
                       "Incidence rate of Shigella diarrhea",
                       ifelse(parm == "IR_ETEC", 
                              "Incidence rate of ETEC diarrhea",
                              "Incidence rate of other diarrhea")),
         val = sprintf("%.3f", val))

# shigella subclinical ----
shig.sub.tbl <- shigella_sub %>%
  mutate(parm = ifelse(prev.diarrhea == 1,
                       "Probability of detection of subclinical Shigella per month, diarrhea in prior month",
                       "Probability of detection of subclinical Shigella per month, no diarrhea in prior"),
         val = sprintf("%.3f", predicted_prob)) %>%
  select(country_id, agegrp, parm, val)

# ETEC subclinical ----
ETEC.sub.tbl <- ETEC_sub %>%
  mutate(parm = ifelse(prev.diarrhea == 1,
                       "Probability of detection of subclinical ETEC per month, diarrhea in prior month",
                       "Probability of detection of subclinical ETEC per month, no diarrhea in prior"),
         val = sprintf("%.3f", predicted_prob)) %>%
  select(country_id, agegrp, parm, val)

# other subclinical ----
other.sub.tbl <- other_sub %>%
  mutate(parm = ifelse(other.diar == 1,
                       "Probability of detection of non-Shigella pathogen per month, other diarrhea this month",
                       "Probability of detection of non-Shigella pathogen per month, no other diarrhea this month"),
         val  = sprintf("%.3f", predicted_prob)) %>%
  select(country_id, agegrp, parm, val)

# NOTE: subclinical culture probability is now part of the unified culture
# model below (culture.tbl), which uses Ct value and shig.diar status as
# predictors rather than a separate agegrp-only subclinical model.

# severity: score mean/SD from lookup table ----
shig.score.tbl <- sev_params$score_params %>%
  mutate(parm = "Modified Vesikari score of Shigella diarrhea (mean, sd)",
         val  = paste0(sprintf("%.2f", mean_score),
                       " (", sprintf("%.2f", sd_score), ")")) %>%
  select(country_id, agegrp, parm, val)

ETEC.score.tbl <- ETEC_sev_params$score_params %>%
  mutate(parm = "Modified Vesikari score of ETEC diarrhea (mean, sd)",
         val  = paste0(sprintf("%.2f", mean_score),
                       " (", sprintf("%.2f", sd_score), ")")) %>%
  select(country_id, agegrp, parm, val)

other.score.tbl <- other_sev_params$score_params %>%
  mutate(parm = "Modified Vesikari score of other diarrhea (mean, sd)",
         val  = paste0(sprintf("%.2f", mean_score),
                       " (", sprintf("%.2f", sd_score), ")")) %>%
  select(country_id, agegrp, parm, val)

# gemsdef and culture: predicted probabilities from GLM coefs
# computed at each country x agegrp combination using plogis()
# helper: get predicted probability from a coef vector for a given cell.
# Uses named lookup (robust to coefficient order / reference levels),
# matching the fix applied in 2-functions.R / 2-functions.cpp.
# extra_name/extra_predictor are vectorized: pass equal-length vectors for
# models with more than one continuous/binary predictor beyond agegrp/country
# (e.g. the unified culture model uses c("shigella_eiec", "shig.diar")).
glm_pred <- function(coefs, age_group, country, extra_name = NULL, extra_predictor = 0) {
  lp <- coefs["(Intercept)"]

  age_key <- paste0("agegrp", age_group)
  if (age_key %in% names(coefs)) lp <- lp + coefs[age_key]

  country_key <- paste0("country_id", country)
  if (country_key %in% names(coefs)) lp <- lp + coefs[country_key]

  if (!is.null(extra_name)) {
    for (i in seq_along(extra_name)) {
      if (extra_name[i] %in% names(coefs))
        lp <- lp + coefs[extra_name[i]] * extra_predictor[i]
    }
  }

  unname(plogis(lp))
}

# P(gemsdef) at mean score for each cell
shig.gems.tbl <- sev_params$score_params %>%
  rowwise() %>%
  mutate(
    parm = "Probability of Shigella diarrhea meeting GEMS MSD criteria (at mean score)",
    val  = sprintf("%.3f",
                   glm_pred(sev_params$gems_params$coefs,
                            agegrp, as.character(country_id),
                            extra_name = "score", extra_predictor = mean_score))
  ) %>%
  ungroup() %>%
  select(country_id, agegrp, parm, val)

ETEC.gems.tbl <- ETEC_sev_params$score_params %>%
  rowwise() %>%
  mutate(
    parm = "Probability of ETEC diarrhea meeting GEMS MSD criteria (at mean score)",
    val  = sprintf("%.3f",
                   glm_pred(ETEC_sev_params$gems_params$coefs,
                            agegrp, as.character(country_id),
                            extra_name = "score", extra_predictor = mean_score))
  ) %>%
  ungroup() %>%
  select(country_id, agegrp, parm, val)

other.gems.tbl <- other_sev_params$score_params %>%
  rowwise() %>%
  mutate(
    parm = "Probability of other diarrhea meeting GEMS MSD criteria (at mean score)",
    val  = sprintf("%.3f",
                   glm_pred(other_sev_params$gems_params$coefs,
                            agegrp, as.character(country_id),
                            extra_name = "score", extra_predictor = mean_score))
  ) %>%
  ungroup() %>%
  select(country_id, agegrp, parm, val)

# P(culture+) at representative Ct values, crossed with shig.diar status.
# Uses the unified culture model (agegrp + country_id + shigella_eiec +
# shig.diar), evaluated at a "low Ct" (high pathogen quantity, Ct = 20) and
# a "high Ct" (low pathogen quantity, Ct = 30) representative value for
# each combination of diarrhea status, so the table shows how culture
# positivity varies with both Ct and Shigella-attributable diarrhea status.
culture.tbl <- expand.grid(
  country_id = countries,
  agegrp     = age_groups,
  shig.diar  = c(0, 1),
  ct_value   = c(20, 30),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    parm = sprintf(
      "Probability of positive culture, %s, Ct=%d",
      ifelse(shig.diar == 1, "Shigella diarrhea", "subclinical Shigella"),
      ct_value
    ),
    val  = sprintf("%.3f",
                   glm_pred(sev_params$culture_params$coefs,
                            agegrp, country_id,
                            extra_name = c("shigella_eiec", "shig.diar"),
                            extra_predictor = c(ct_value, shig.diar)))
  ) %>%
  ungroup() %>%
  select(country_id, agegrp, parm, val)

ETEC.culture.tbl <- expand.grid(
  country_id = countries,
  agegrp     = age_groups,
  ETEC.diar  = c(0, 1),
  ct_value   = c(20, 30),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    parm = sprintf(
      "Probability of positive culture, %s, Ct=%d",
      ifelse(ETEC.diar == 1, "ETEC diarrhea", "subclinical ETEC"),
      ct_value
    ),
    val  = sprintf("%.3f",
                   glm_pred(ETEC_sev_params$culture_params$coefs,
                            agegrp, country_id,
                            extra_name = c("ST_ETEC", "ETEC.diar"),
                            extra_predictor = c(ct_value, ETEC.diar)))
  ) %>%
  ungroup() %>%
  select(country_id, agegrp, parm, val)

# shigella pathogen quantity ----
# Gamma GLM (log link) fit on country x agegrp x severity. shape_quantity is
# the (constant-across-cells) Gamma shape parameter; reported alongside the
# mean since, unlike a Normal model, the implied SD is mean-dependent
# (SD = mean / sqrt(shape)) rather than fixed.
shig.quantity.tbl <- shigella_quantity %>%
  mutate(parm = case_when(
           severity == "Subclinical" ~ "Subclinical Shigella quantity (mean, shape)",
           severity == "Mild"        ~ "Mild Shigella diarrhea quantity (mean, shape)",
           severity == "Severe"      ~ "Severe Shigella diarrhea quantity (mean, shape)"),
         parm = factor(parm, levels = c("Subclinical Shigella quantity (mean, shape)",
                                        "Mild Shigella diarrhea quantity (mean, shape)",
                                        "Severe Shigella diarrhea quantity (mean, shape)")),
         val  = paste0(sprintf("%.2f", mean_quantity),
                       " (", sprintf("%.2f", shape_quantity), ")")) %>%
  arrange(country_id, agegrp, parm) %>%
  select(country_id, agegrp, parm, val)

# ETEC pathogen quantity ----
# Same Gamma GLM structure as Shigella above.
ETEC.quantity.tbl <- ETEC_quantity %>%
  mutate(parm = case_when(
    severity == "Subclinical" ~ "Subclinical ETEC quantity (mean, shape)",
    severity == "Mild"        ~ "Mild ETEC diarrhea quantity (mean, shape)",
    severity == "Severe"      ~ "Severe ETEC diarrhea quantity (mean, shape)"),
    parm = factor(parm, levels = c("Subclinical ETEC quantity (mean, shape)",
                                   "Mild ETEC diarrhea quantity (mean, shape)",
                                   "Severe ETEC diarrhea quantity (mean, shape)")),
    val  = paste0(sprintf("%.2f", mean_quantity),
                  " (", sprintf("%.2f", shape_quantity), ")")) %>%
  arrange(country_id, agegrp, parm) %>%
  select(country_id, agegrp, parm, val)

# other pathogen quantity ----
# Same Gamma GLM structure as Shigella above.
other.quantity.tbl <- other_quantity %>%
  mutate(parm = case_when(
           severity == "Subclinical" ~ "Subclinical other-pathogen quantity (mean, shape)",
           severity == "Mild"        ~ "Mild other diarrhea quantity (mean, shape)",
           severity == "Severe"      ~ "Severe other diarrhea quantity (mean, shape)"),
         parm = factor(parm, levels = c("Subclinical other-pathogen quantity (mean, shape)",
                                        "Mild other diarrhea quantity (mean, shape)",
                                        "Severe other diarrhea quantity (mean, shape)")),
         val  = paste0(sprintf("%.2f", mean_quantity),
                       " (", sprintf("%.2f", shape_quantity), ")")) %>%
  arrange(country_id, agegrp, parm) %>%
  select(country_id, agegrp, parm, val)

# join all ----
all.parms <- ir.tbl %>%
  add_row(shig.sub.tbl) %>%
  add_row(ETEC.sub.tbl) %>%
  add_row(other.sub.tbl) %>%
  add_row(shig.score.tbl) %>%
  add_row(ETEC.score.tbl) %>%
  add_row(other.score.tbl) %>%
  add_row(shig.gems.tbl) %>%
  add_row(ETEC.gems.tbl) %>%
  add_row(other.gems.tbl) %>%
  add_row(culture.tbl) %>%
  add_row(ETEC.culture.tbl) %>%
  add_row(shig.quantity.tbl) %>%
  add_row(ETEC.quantity.tbl) %>%
  add_row(other.quantity.tbl) %>%
  filter(agegrp %in% age_groups & country_id %in% countries) %>%
  arrange(country_id, agegrp) %>%
  pivot_wider(names_from = c("country_id", "agegrp"), values_from = "val")

# save ----
write.csv(all.parms, file = here::here("tables", "parms.csv"), row.names = FALSE)
