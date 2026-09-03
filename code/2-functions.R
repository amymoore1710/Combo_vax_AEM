# determine whether diarrhea occurs on a given month
diarrhea_event <- function(lambda) {
  if (rexp(1, lambda) <= 1) return(1) else return(0)
}

# pull the corresponding lambdas for shigella and other diarrhea
diarrhea_lambda <- function(data, vax_status, age_group, country) {
  row <- data[data$country_id == country & data$agegrp == age_group, ]
  ir_shigella <- if (vax_status == 1) row$vax_IR_shigella else row$IR_shigella
  ir_ETEC <- if (vax_status == 1) row$vax_IR_ETEC else row$IR_ETEC
  c(ir_shigella = ir_shigella, ir_other = row$IR_other)
}

# apply severity-dependent VE shrinkage to a drawn score for vaccinated
# children. Linear ramp: no effect at score = 1, maximal effect (ve_max)
# at score = 12. ve_max is calibrated in 1-parameters.R (calibrate_ve_max())
# so the resulting VE against GEMS MSD matches the target ve_severity.
shrink_score <- function(score, ve_max) {
  ve_at_score <- ve_max * (score - 1) / 11
  round(pmin(pmax(score * (1 - ve_at_score), 1), 12))
}

# simulate severity score from observed mean/SD table, then apply
# severity-dependent shrinkage for vaccinated children. The draw itself is
# identical for vax/unvax — the VE effect is applied after, conditional on
# the realized score.
simulate_score <- function(score_params, vax_status, age_group, country, ve_max) {
  row <- score_params[score_params$country_id == country &
                        score_params$agegrp == age_group, ]
  drawn_score <- round(pmin(pmax(rnorm(1, row$mean_score, row$sd_score), 1), 12))
  if (vax_status == 1) shrink_score(drawn_score, ve_max) else drawn_score
}

# simulate gemsdef from logistic regression GLM coefs
simulate_gemsdef <- function(gems_params, score, age_group, country) {
  lp <- logodds_from_coefs(gems_params$coefs, age_group, country,
                            extra_name = "score", extra_predictor = score)
  rbinom(1, 1, plogis(lp))
}

# Upper bound on simulated pathogen quantity: a Ct of 0 is the
# (theoretical) detection ceiling, so quantity should not exceed
# (35 - 0) / 3.322. The Gamma quantity models have no natural upper bound
# (unlike Ct, which is physically constrained to be >= 0), so without this
# cap a small fraction of draws — concentrated in higher-severity Shigella
# cells (~1-3%) — exceed it and convert to a negative, physically
# implausible Ct. Mirrors the existing lower-bound floor (quantity > 0).
MAX_QUANTITY <- 35 / 3.322

# convert simulated pathogen quantity to Ct for use as a predictor in the
# culture model (which was fit on Ct directly: shigella_micro ~ ... +
# shigella_eiec + shig.diar). Quantity is simulated rather than Ct itself
# so that no artificial truncation near the Ct=35 detection ceiling is
# needed — quantity only needs a natural lower bound at 0.
# Conversion: Ct = 35 - 3.322 * quantity (inverse of quantity = (35-Ct)/3.322
# used when estimating shigella_quantity.RDS in 0-estimate-params.R)
quantity_to_ct <- function(quantity) {
  35 - 3.322 * quantity
}

# simulate Shigella culture positivity from the unified logistic regression
# GLM coefs. Used for ANY PCR-positive Shigella specimen (clinical diarrhea
# or subclinical detection), predicted from agegrp, country, the drawn
# pathogen quantity (converted to Ct), and whether the episode is
# Shigella-attributable diarrhea (shig.diar). Replaces the old two-path
# approach (gemsdef-based for clinical cases, agegrp-only for subclinical).
simulate_culture <- function(culture_params, quantity, shig_diar, age_group, country) {
  ct_value <- quantity_to_ct(quantity)
  lp <- logodds_from_coefs(culture_params$coefs, age_group, country,
                            extra_name = c("shigella_eiec", "shig.diar"),
                            extra_predictor = c(ct_value, shig_diar))
  rbinom(1, 1, plogis(lp))
}

# simulate other-diarrhea score (no vax effect)
simulate_other_score <- function(score_params, age_group, country) {
  row <- score_params[score_params$country_id == country &
                        score_params$agegrp == age_group, ]
  round(pmin(pmax(rnorm(1, row$mean_score, row$sd_score), 1), 12))
}

# simulate other-diarrhea gemsdef (no vax effect)
simulate_other_gemsdef <- function(gems_params, score, age_group, country) {
  lp <- logodds_from_coefs(gems_params$coefs, age_group, country,
                            extra_name = "score", extra_predictor = score)
  rbinom(1, 1, plogis(lp))
}

# shared helper: build log-odds from a named GLM coef vector.
# Looks up all terms by name so it is robust to any number of coefficients,
# any reference level, and any coefficient ordering.
# extra_name: name(s) of continuous/binary predictor(s) beyond agegrp/country
#             (e.g. "score" for gemsdef model, or c("shigella_eiec",
#             "shig.diar") for the unified culture model). Vectorized:
#             extra_name and extra_predictor must be the same length.
logodds_from_coefs <- function(coefs, age_group, country,
                                extra_name, extra_predictor) {
  lp <- coefs["(Intercept)"]

  # age group dummy — construct key from value, add only if present in model
  age_key <- paste0("agegrp", age_group)
  if (age_key %in% names(coefs)) lp <- lp + coefs[age_key]

  # country dummy — construct key from value, add only if present in model
  country_key <- paste0("country_id", country)
  if (country_key %in% names(coefs)) lp <- lp + coefs[country_key]

  # continuous/binary predictors — looked up by explicit name, not position
  for (i in seq_along(extra_name)) {
    if (extra_name[i] %in% names(coefs)) {
      lp <- lp + coefs[extra_name[i]] * extra_predictor[i]
    }
  }

  unname(lp)
}

# simulate pathogen quantity for a Shigella infection. quantity ~ country x
# agegrp x severity (Subclinical/Mild/Severe), fit as a Gamma GLM (log link)
# with a constant shape parameter across cells (see 0-estimate-params.R).
# Drawn from rgamma() using shape/scale parameterization (scale = mean/shape
# so that the draw's mean matches the GLM-predicted mean_quantity). Quantity
# is clamped only at a strict lower bound (>0, enforced via a small positive
# floor) since simulating an actual infection implies quantity > 0; no
# upper bound is needed (unlike the old Ct-based approach which needed an
# artificial cap at 34.999 near the Ct=35 detection ceiling).
assign_shig_quantity <- function(data, vax_status, age_group, country, severity) {
  row <- data[data$country_id == country & data$agegrp == age_group &
                data$severity == severity, ]
  mu  <- if (vax_status == 1) row$vax_mean_quantity else row$mean_quantity
  pmin(pmax(rgamma(1, shape = row$shape_quantity, scale = mu / row$shape_quantity), 1e-6), MAX_QUANTITY)
}

# simulate pathogen quantity for an other-pathogen infection. Same Gamma GLM
# structure as Shigella (quantity ~ country x agegrp x severity), no vax
# effect. Quantity is clamped only at a strict lower bound (>0, enforced via
# a small positive floor) since simulating an actual infection implies
# quantity > 0; no upper bound is needed (unlike the old Ct-based approach
# which needed an artificial cap at 34.999 near the Ct=35 detection
# ceiling). Same conversion formula (Ct = 35 - 3.322*quantity) is used
# across all pathogens, including Shigella.
assign_other_quantity <- function(data, age_group, country, severity) {
  row <- data[data$country_id == country & data$agegrp == age_group &
                data$severity == severity, ]
  pmin(pmax(rgamma(1, shape = row$shape_quantity, scale = row$mean_quantity / row$shape_quantity), 1e-6), MAX_QUANTITY)
}

# detecting subclinical shigella
shigella_det <- function(data, vax_status, age_group, country, prev.diarrhea) {
  row  <- data[data$country_id == country & data$agegrp == age_group &
                 data$prev.diarrhea == prev.diarrhea, ]
  prob <- if (vax_status == 1) row$vax_predicted_prob else row$predicted_prob
  rbinom(1, 1, prob)
}

# detecting other subclinical infections
other_det <- function(data, age_group, country, other.diar) {
  row <- data[data$country_id == country & data$agegrp == age_group &
                data$other.diar == other.diar, ]
  rbinom(1, 1, row$predicted_prob)
}

# simulate trial
simulate_trial <- function(n_children, n_months, countries,
                           diarrhea_ir,
                           shigella_sub, other_sub,
                           sev_params,
                           other_sev_params,
                           shigella_quantity, other_quantity) {

  # unpack severity parameter objects
  shig_score_params  <- sev_params$score_params
  gems_params        <- sev_params$gems_params
  culture_params     <- sev_params$culture_params
  other_score_params <- other_sev_params$score_params
  other_gems_params  <- other_sev_params$gems_params

  # calibrated severity-dependent VE shrinkage parameter (single scalar,
  # computed once per scenario in f.param() / calibrate_ve_max()).
  # Defaults to 0 (no shrinkage) if not present.
  ve_max <- if (!is.null(sev_params$ve_max)) sev_params$ve_max else 0

  # base trial structure
  trial <- data.frame(
    country_id = rep(countries, each = n_months * n_children),
    child      = rep(1:n_children, each = n_months),
    month      = rep(1:n_months, times = n_children * length(countries))
  ) %>%
    mutate(agegrp = ifelse(month <= 6, "12-17 months", "18-24 months")) %>%
    group_by(country_id, child) %>%
    mutate(vax = rbinom(1, 1, 0.5)) %>%
    ungroup()

  # diarrhea events
  trial <- trial %>%
    rowwise() %>%
    mutate(
      lambdas           = list(diarrhea_lambda(diarrhea_ir, vax, agegrp, country_id)),
      shigella_diarrhea = diarrhea_event(lambdas[[1]]["ir_shigella"]),
      other_diarrhea    = diarrhea_event(lambdas[[1]]["ir_other"])
    ) %>%
    ungroup()

  # shigella severity chain: score -> gemsdef (culture moved below, after
  # Ct is known, since the new unified culture model depends on the drawn
  # Ct value rather than gemsdef status)
  trial <- trial %>%
    rowwise() %>%
    mutate(
      shigella_score = if (shigella_diarrhea == 1)
        simulate_score(shig_score_params, vax, agegrp, country_id, ve_max)
      else NA_real_,

      shigella_gemsmsd = if (shigella_diarrhea == 1)
        simulate_gemsdef(gems_params, shigella_score, agegrp, country_id)
      else NA_integer_,

      shigella_any_severity = if (shigella_diarrhea == 1)
        ifelse(shigella_score >= 6 | shigella_gemsmsd == 1, "Severe", "Mild")
      else NA_character_
    ) %>%
    ungroup()

  # other diarrhea severity: score -> gemsdef (no culture)
  trial <- trial %>%
    rowwise() %>%
    mutate(
      other_score = if (other_diarrhea == 1)
        simulate_other_score(other_score_params, agegrp, country_id)
      else NA_real_,

      other_gemsmsd = if (other_diarrhea == 1)
        simulate_other_gemsdef(other_gems_params, other_score, agegrp, country_id)
      else NA_integer_,

      other_any_severity = if (other_diarrhea == 1)
        ifelse(other_score >= 6 | other_gemsmsd == 1, "Severe", "Mild")
      else NA_character_
    ) %>%
    ungroup()

  # subclinical shigella detection (PCR detection, no diarrhea)
  trial <- trial %>%
    group_by(country_id, child) %>%
    mutate(prev.diarrhea = ifelse(lag(shigella_diarrhea, default = 0) == 1, 1, 0)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(shigella_sub = shigella_det(shigella_sub, vax, agegrp, country_id,
                                       prev.diarrhea)) %>%
    ungroup()

  # shigella pathogen quantity: drawn for ANY PCR-detectable Shigella this
  # month (clinical diarrhea OR subclinical detection) — must come before
  # culture status since the unified culture model uses Ct (converted from
  # quantity) directly.
  trial <- trial %>%
    rowwise() %>%
    mutate(shigella_quantity = case_when(
      shigella_diarrhea == 1 ~ assign_shig_quantity(shigella_quantity, vax, agegrp, country_id,
                                               shigella_any_severity),
      shigella_sub      == 1 ~ assign_shig_quantity(shigella_quantity, vax, agegrp, country_id,
                                               "Subclinical"),
      TRUE                   ~ NA_real_
    )) %>%
    ungroup()

  # shigella culture: unified model across ALL PCR-positive specimens
  # (clinical or subclinical), predicted from agegrp, country, the drawn
  # Ct value (converted from quantity), and shig.diar status. Replaces the
  # old two-path approach (gemsdef-based for clinical, agegrp-only for
  # subclinical).
  trial <- trial %>%
    rowwise() %>%
    mutate(
      shig_pcr_positive = (shigella_diarrhea == 1) | (shigella_sub == 1),
      shig_diar_flag    = ifelse(shigella_diarrhea == 1, 1, 0),
      shigella_culture  = if (shig_pcr_positive)
        simulate_culture(culture_params, shigella_quantity, shig_diar_flag, agegrp, country_id)
      else NA_integer_
    ) %>%
    select(-shig_pcr_positive, -shig_diar_flag) %>%
    ungroup()

  # other pathogen detection and quantity: detection probability now
  # depends on other_diarrhea (see 0-estimate-params.R), mirroring how
  # Shigella subclinical detection depends on prev.diarrhea. Quantity uses
  # the same severity structure as Shigella above — "Subclinical" when
  # detected outside of an other-diarrhea episode, otherwise the drawn
  # other-diarrhea severity (Mild/Severe).
  trial <- trial %>%
    rowwise() %>%
    mutate(
      other_inf = other_det(other_sub, agegrp, country_id, other_diarrhea),
      other_quantity = case_when(
        is.na(other_inf) || other_inf != 1 ~ NA_real_,
        other_diarrhea == 1 ~ assign_other_quantity(other_quantity, agegrp, country_id,
                                                     other_any_severity),
        TRUE ~ assign_other_quantity(other_quantity, agegrp, country_id, "Subclinical")
      )
    ) %>%
    ungroup()

  return(trial)
}