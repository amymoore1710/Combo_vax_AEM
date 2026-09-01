pacman::p_load(dplyr)

# Simulation settings ----
months_recruitment <- 12 #recruitment period for study

# Trial structure parameters ----
n_months <- 12  # Trial duration in months
age_groups <- c("12-17 months", "18-24 months")
countries <- c("BG", "PE", "PK", "IN")

# MAL-ED parameters ----
## IR of shigella and other diarrhea per child month
diarrhea_ir <- readRDS(here::here("sim param data","diarrhea_incidence.RDS"))
## Probability of shigella subclinical infection per child month
shigella_sub <- readRDS(here::here("sim param data","shigella_subclinical_prob.RDS"))
## Probability of other pathogen infection per child month
other_sub <- readRDS(here::here("sim param data","other_subclinical_prob.RDS"))
## Shigella severity: score mean/SD table + GLM coefs for gemsdef and culture
sev_params <- readRDS(here::here("sim param data","shigella_severity.RDS"))
## Other diarrhea severity: score mean/SD table + GLM coefs for gemsdef
other_sev_params <- readRDS(here::here("sim param data","other_severity.RDS"))
## Shigella pathogen quantity: Gamma GLM (log link) fit on country x agegrp
## x severity (Subclinical/Mild/Severe), with a constant shape parameter
## (quantity is simulated directly, then converted to Ct only where Ct is
## needed, e.g. as a predictor in the culture model — see quantity_to_ct()
## in 2-functions.cpp/2-functions.R for the conversion formula)
shigella_quantity <- readRDS(here::here("sim param data","shigella_quantity.RDS"))
## Other pathogen quantity: same Gamma GLM structure as Shigella, fit on
## country x agegrp x severity (Subclinical/Mild/Severe)
## (quantity is simulated directly, then converted to Ct for output —
## see quantity_to_ct() conversion formula in 2-functions.cpp/2-functions.R)
other_quantity <- readRDS(here::here("sim param data","other_quantity.RDS"))

# Vaccine efficacy assumptions ----
ve_infection <- 0.10  # 10% VE against infection
ve_disease   <- 0.40  # 40% VE against disease
ve_severity  <- 0.60  # target AVERAGE marginal VE against MSD, across both
                       # definitions ("score >= 6" and GEMS MSD) — achieved
                       # via a severity-dependent score shrinkage calibrated
                       # below (see calibrate_ve_max())
ve_ct        <- 0.10  # reduction in pathogen quantity among breakthrough cases

# ---------------------------------------------------------------
# Severity-dependent VE on score
#
# Rather than shifting the mean of the whole score distribution uniformly,
# the vaccine's effect on severity scales linearly with how severe the score
# would have been: no effect at score = 1, maximal effect (ve_max) at
# score = 12. This is applied AFTER the score is drawn (the shrinkage factor
# depends on the realized draw), not before, so there is no vax_mean_score
# lookup anymore — vax_status is decided downstream of the random draw.
#
#   ve_at_score(s) = ve_max * (s - 1) / 11
#   vax_score      = round(clamp(s * (1 - ve_at_score(s)), 1, 12))
#
# ve_max is NOT the same as the target average marginal VE against MSD
# (ve_severity) because (a) the shrinkage acts proportionally on the score,
# not as a flat probability shift, and (b) MSD truth (score-based or
# GEMS-based) requires BOTH a Shigella diarrhea episode AND severity, so
# the case-severity effect this shrinkage produces compounds multiplicatively
# with the separate ve_disease effect on overall Shigella diarrhea incidence
# — see calibrate_ve_max() for the exact marginal_ve = 1-(1-ve_disease)*
# (1-conditional_ve) relationship. We calibrate ve_max numerically so that
# the AVERAGE of the two resulting marginal VEs (score >= 6 and GEMS MSD),
# each incidence-weighted across country/agegrp cells, equals ve_severity.
#
# The underlying shrinkage mechanism still only acts on the score directly
# (not on the GEMS MSD probability) — GEMS MSD's conditional VE is whatever
# falls out of feeding the shrunk score into the gemsdef GLM, capped by that
# GLM's own intercept/country/age effects, which score shrinkage can't
# touch. So while ve_max is searched over the score mechanism, the
# calibration objective now folds in the GEMS side's achieved marginal VE
# too — averaging the two rather than targeting "score >= 6" alone as before.
#
# Known limitation — score >= 6 conditional VE bias (worst for Pakistan)
# -----------------------------------------------------------------------
# NOTE: the numbers below (ve_max ~= 0.39, the score-crossing table, and
# the per-site frac_rescuable/VE table) were derived under the OLD
# case-conditional-only calibration (target_ve_severity applied directly
# to the score >= 6 conditional VE). Under the current marginal-averaged
# calibration, the conditional VE this mechanism needs to hit is lower —
# solving 1-(1-ve_disease)*(1-cond) = ve_severity with ve_disease = 0.40,
# ve_severity = 0.60 gives cond ~= 0.33, not 0.60 — so the calibrated
# ve_max is expected to come out smaller than 0.39, and the exact
# crossing pattern/numbers below will shift accordingly. Treat this
# section as illustrative of the MECHANISM (integer rounding means only
# scores near the threshold get rescued, and how far each site's score
# distribution sits above 6 drives how much conditional VE is achievable)
# rather than as current numbers — re-derive after running the calibration
# once to get updated figures.
#
# With ve_max calibrated to ~0.39, integer rounding means only one score
# value crosses the threshold: score 6 → vax_score 5. Scores 7–12 all
# shrink but remain >= 6:
#
#   score  vax_score (ve_max ≈ 0.39)   crosses threshold?
#   6      5                             YES
#   7      6                             no
#   8      6                             no
#   9–12   6–7                           no
#
# This applies to ALL sites identically (ve_max is a single scalar). The
# consequence is that the true VE against score >= 6 reduces to an identity:
#
#   VE ≈ P(score = 6) / P(score >= 6)  =  "frac_rescuable"
#
# i.e., the fraction of above-threshold cases sitting exactly at score = 6.
# That fraction — and therefore the achievable VE — varies considerably
# across sites based on how far each score distribution extends above 6:
#
#   site          mean_score  frac_rescuable  true VE vs score >= 6
#   BG 12-17mo    3.67        ~67%            ~66%
#   BG 18-24mo    3.42        ~70%            ~70%
#   IN 12-17mo    4.58        ~42%            ~42%   <- also notably low
#   IN 18-24mo    2.87        ~69%            ~69%
#   PE 12-17mo    3.72        ~60%            ~60%
#   PE 18-24mo    3.77        ~52%            ~52%
#   PK 12-17mo    6.27        ~25%            ~25%   <- worst affected
#   PK 18-24mo    4.98        ~35%            ~34%
#
# Pakistan is the most affected because its score distributions are centered
# near or above the threshold (mean 6.27 at 12-17mo, 4.98 at 18-24mo) —
# about 75% of its score >= 6 cases have scores 7–12 and cannot be rescued.
# India (12-17mo) is also notably low at ~42%. Bangladesh achieves the
# highest rescue rates (~67–70%) because its low mean scores (~3.4–3.7)
# concentrate above-threshold cases right at score = 6.
#
# The VE bias seen in output reflects this: f.aggregate_ve()'s unsuffixed
# "True VE"/"VE Bias" measures use the average marginal VE across the four
# real sites (excluding the pooled "All Sites" row) as the benchmark for
# severity endpoints (to avoid circularity with country-specific true VEs
# that differ by design) — see sim-analysis-functions.R. Sites with low
# frac_rescuable (PK, IN 12-17mo) will tend to show negative VE bias for
# the score >= 6 family against that benchmark; sites with high
# frac_rescuable (BG) will tend to show positive VE bias. This is a
# structural modeling constraint, NOT trial inefficiency or diagnostic
# misclassification. Switching to equal site weights in the calibration
# does not resolve it: there is no single ve_max that brings all sites'
# conditional VE near the same value, because raising ve_max enough for
# Pakistan would push Bangladesh, Peru, and India much higher still.
#
# The GEMS MSD endpoint is less affected at all sites because GEMS MSD is a
# continuous logit: any score reduction (even one that keeps score >= 6)
# lowers the predicted GEMS MSD probability, so the vaccine's effect
# accumulates across the full score distribution rather than only at the
# 6→5 boundary.
# ---------------------------------------------------------------

# Apply the linear severity-dependent shrinkage to a vector of scores
shrink_score <- function(score, ve_max) {
  ve_at_score <- ve_max * (score - 1) / 11
  round(pmin(pmax(score * (1 - ve_at_score), 1), 12))
}

# Simulate the score -> shrink chain (vectorized, no full trial needed) to
# estimate the CONDITIONAL (case-severity) VE against "score >= 6" diarrhea
# for a candidate ve_max, for one country/agegrp cell — i.e. P(score >= 6)
# among children who already have Shigella diarrhea, vax vs. unvax. This
# targets the score threshold directly, rather than going through the
# gemsdef GLM — score >= 6 is fully determined by the (post-shrinkage)
# score itself, so there's no GLM-intercept floor limiting how high the
# achievable conditional VE can go (unlike the GEMS MSD mechanism below).
simulate_score_ve_cell <- function(ve_max, mean_score, sd_score, n_draws = 200000) {

  unvax_score <- round(pmin(pmax(rnorm(n_draws, mean_score, sd_score), 1), 12))
  vax_score   <- shrink_score(unvax_score, ve_max)

  p_unvax <- mean(unvax_score >= 6)
  p_vax   <- mean(vax_score   >= 6)

  list(p_unvax = p_unvax, p_vax = p_vax)
}

# Log-odds from a named GLM coef vector — mirrors logodds_from_coefs() in
# 2-functions.R / 2-functions.cpp exactly. Duplicated here (rather than
# sourced from 2-functions.R) because 1-parameters.R runs first in the
# pipeline and calibrate_ve_max() needs it to simulate gemsdef probabilities
# during calibration. Vectorized over `score_vec` (unlike the original,
# which handles one child at a time) since calibration works with a whole
# vector of Monte Carlo score draws per country/agegrp cell at once.
gemsdef_logodds_vec <- function(gems_coefs, age_group, country, score_vec) {
  lp <- gems_coefs["(Intercept)"]

  age_key <- paste0("agegrp", age_group)
  if (age_key %in% names(gems_coefs)) lp <- lp + gems_coefs[age_key]

  country_key <- paste0("country_id", country)
  if (country_key %in% names(gems_coefs)) lp <- lp + gems_coefs[country_key]

  score_coef <- if ("score" %in% names(gems_coefs)) gems_coefs["score"] else 0

  unname(lp + score_coef * score_vec)
}

# Simulate the score -> gemsdef chain (vectorized) to estimate the
# CONDITIONAL (case-severity) VE against GEMS MSD for a candidate ve_max,
# for one country/agegrp cell — i.e. P(gemsdef == 1) among children who
# already have Shigella diarrhea, vax vs. unvax. Uses the analytic mean
# P(gemsdef==1) = mean(plogis(logodds)) over the drawn score distribution
# rather than an extra layer of Bernoulli draws, since only the expectation
# is needed here (lower Monte Carlo noise for the same n_draws than
# rbinom(n_draws, 1, plogis(lp)) would give).
simulate_gems_ve_cell <- function(ve_max, mean_score, sd_score, gems_coefs,
                                  age_group, country, n_draws = 200000) {

  unvax_score <- round(pmin(pmax(rnorm(n_draws, mean_score, sd_score), 1), 12))
  vax_score   <- shrink_score(unvax_score, ve_max)

  lp_unvax <- gemsdef_logodds_vec(gems_coefs, age_group, country, unvax_score)
  lp_vax   <- gemsdef_logodds_vec(gems_coefs, age_group, country, vax_score)

  p_unvax <- mean(plogis(lp_unvax))
  p_vax   <- mean(plogis(lp_vax))

  list(p_unvax = p_unvax, p_vax = p_vax)
}

# ---------------------------------------------------------------
# Calibrate ve_max so that the AVERAGE of the two MARGINAL (population-
# level) VEs against MSD — "score >= 6" and GEMS MSD — equals
# target_ve_severity.
#
# "Marginal" here means VE against the truth columns actually used
# downstream (sev.score.shigella_diarrhea / sev.gems.shigella_diarrhea =
# shigella_diarrhea==1 & <severity criterion>), which require BOTH a
# Shigella diarrhea episode occurring AND that episode being severe. Since
# vaccination reduces both independently — overall Shigella diarrhea
# incidence via ve_disease (vax_IR_shigella = IR_shigella*(1-ve_disease)),
# and conditional case-severity via this shrinkage mechanism — the two
# effects compound multiplicatively rather than target_ve_severity being
# the achieved marginal VE directly:
#   marginal_ve = 1 - (1 - ve_disease) * (1 - conditional_ve)
# This function inverts that relationship (per severity definition, then
# averaged) to solve for the ve_max that makes the AVERAGE marginal VE hit
# target_ve_severity, folding ve_disease into the objective rather than
# treating target_ve_severity as a case-conditional target as before.
# ---------------------------------------------------------------
calibrate_ve_max <- function(target_ve_severity, score_params, gems_coefs,
                             diarrhea_ir, ve_disease, n_draws = 200000) {

  weights <- score_params %>%
    left_join(diarrhea_ir, by = c("country_id", "agegrp")) %>%
    mutate(weight = IR_shigella / sum(IR_shigella))

  achieved_ve_fn <- function(ve_max) {
    cell_results <- weights %>%
      rowwise() %>%
      mutate(
        score_res     = list(simulate_score_ve_cell(ve_max, mean_score, sd_score, n_draws)),
        p_unvax_score = score_res$p_unvax,
        p_vax_score   = score_res$p_vax,
        gems_res      = list(simulate_gems_ve_cell(ve_max, mean_score, sd_score,
                                                    gems_coefs, agegrp, country_id, n_draws)),
        p_unvax_gems  = gems_res$p_unvax,
        p_vax_gems    = gems_res$p_vax
      ) %>%
      ungroup()

    weighted_p_unvax_score <- sum(cell_results$p_unvax_score * cell_results$weight)
    weighted_p_vax_score   <- sum(cell_results$p_vax_score   * cell_results$weight)
    weighted_p_unvax_gems  <- sum(cell_results$p_unvax_gems  * cell_results$weight)
    weighted_p_vax_gems    <- sum(cell_results$p_vax_gems    * cell_results$weight)

    # conditional (case-severity) VE, per definition
    cond_ve_score <- 1 - (weighted_p_vax_score / weighted_p_unvax_score)
    cond_ve_gems  <- 1 - (weighted_p_vax_gems   / weighted_p_unvax_gems)

    # marginal (population-level) VE, per definition — compounds the
    # case-severity effect with the separate ve_disease reduction
    marginal_ve_score <- 1 - (1 - ve_disease) * (1 - cond_ve_score)
    marginal_ve_gems  <- 1 - (1 - ve_disease) * (1 - cond_ve_gems)

    mean(c(marginal_ve_score, marginal_ve_gems))
  }

  ve_objective <- function(ve_max) achieved_ve_fn(ve_max) - target_ve_severity

  # Sanity-check the achievable range before calling uniroot. At ve_max = 0
  # (no case-severity shrinkage at all), the achieved marginal VE is not 0
  # under this formulation — it's exactly ve_disease, since vaccination
  # still reduces overall Shigella diarrhea incidence even with no severity
  # effect. The achievable ceiling as ve_max -> 1 is capped below 100% for
  # GEMS MSD by its GLM intercept (unlike "score >= 6", which has no such
  # floor), so the ceiling here is whatever the average of the two
  # definitions reaches, not necessarily near 100%.
  ve_at_0    <- achieved_ve_fn(0)
  ve_at_ceil <- achieved_ve_fn(0.999)

  if (ve_at_ceil < target_ve_severity) {
    stop(sprintf(
      paste0(
        "Target average marginal VE against MSD (%.1f%%) is unreachable.\n",
        "  Achieved average marginal VE at ve_max = 0     : %.1f%%\n",
        "  Achieved average marginal VE at ve_max = 0.999 : %.1f%%  <-- ceiling for this mechanism\n",
        "Lower ve_severity below %.1f%%."
      ),
      target_ve_severity * 100, ve_at_0 * 100, ve_at_ceil * 100, ve_at_ceil * 100
    ))
  }

  if (ve_at_0 > target_ve_severity) {
    stop(sprintf(
      paste0(
        "Target average marginal VE against MSD (%.1f%%) is already exceeded at ve_max = 0 ",
        "(%.1f%%), i.e. by ve_disease alone. Raise ve_severity above %.1f%%, or lower ve_disease."
      ),
      target_ve_severity * 100, ve_at_0 * 100, ve_at_0 * 100
    ))
  }

  uniroot(ve_objective, interval = c(0, 0.999), tol = 0.005)$root
}

# ---------------------------------------------------------------
# Cache for calibrated ve_max values, keyed by (target ve_severity,
# ve_disease).
#
# calibrate_ve_max() involves Monte Carlo draws (rnorm with no fixed seed)
# inside uniroot()'s root-finding, so two independent calls with the same
# inputs will converge to slightly different ve_max values due to sampling
# noise. Scenarios that share the same ve_severity target AND ve_disease
# (e.g. "Main", "VE.inf0.0", both at ve_severity = 0.60 / ve_disease = 0.40)
# should use an IDENTICAL ve_max so that comparisons across those scenarios
# isolate the effect of the parameter that actually differs between them.
#
# NOTE: unlike the previous (case-conditional) calibration, ve_disease is
# now folded into the calibration objective itself (marginal VE compounds
# both effects), so scenarios that vary ve_disease at a fixed ve_severity
# (e.g. "VE.dis0.3") will now legitimately get a DIFFERENT ve_max than
# "Main" — that's expected, not a caching bug: a lower ve_disease needs
# more case-severity shrinkage to reach the same average marginal target.
# Caching by (ve_severity, ve_disease) together guarantees identical ve_max
# only for genuinely identical calibration inputs, and still avoids
# re-running the (relatively expensive) calibration search for repeats.
# ---------------------------------------------------------------
.ve_max_cache <- new.env(parent = emptyenv())

get_calibrated_ve_max <- function(ve_severity, score_params, gems_coefs,
                                  diarrhea_ir, ve_disease, n_draws = 200000) {
  cache_key <- paste0("sev", ve_severity, "_dis", ve_disease)

  if (!exists(cache_key, envir = .ve_max_cache, inherits = FALSE)) {
    message(sprintf(
      "Calibrating ve_max for average marginal VE against MSD (score/gems) = %.2f (ve_disease = %.2f)...",
      ve_severity, ve_disease))
    calibrated <- calibrate_ve_max(
      target_ve_severity = ve_severity,
      score_params        = score_params,
      gems_coefs          = gems_coefs,
      diarrhea_ir         = diarrhea_ir,
      ve_disease          = ve_disease,
      n_draws             = n_draws
    )
    assign(cache_key, calibrated, envir = .ve_max_cache)
    message(sprintf("  -> calibrated ve_max = %.4f (cached for %s)",
                    calibrated, cache_key))
  } else {
    message(sprintf("Using cached ve_max for %s -> ve_max = %.4f",
                    cache_key, get(cache_key, envir = .ve_max_cache)))
  }

  get(cache_key, envir = .ve_max_cache, inherits = FALSE)
}

# ---------------------------------------------------------------
# Analytic (very-low-noise) true VE against MSD, computed directly from the
# underlying DGP parameters rather than by simulating a trial-sized sample.
# Reuses the same score -> shrink -> gemsdef mechanism as calibrate_ve_max()
# (simulate_score_ve_cell()/simulate_gems_ve_cell()), but:
#   - run with a much larger n_draws (default 1e6) than the calibration
#     search uses, since this is a one-time calculation rather than
#     something re-run at every uniroot() step — no meaningful cost to
#     making it much more precise
#   - reported per country AND pooled ("ALL", IR_shigella-weighted across
#     every country x agegrp cell — relabeled "All Sites" downstream, same
#     convention as f.aggregate_ve()'s internal country_id), not just an
#     overall average
#   - reported separately for score >= 6 and GEMS MSD, not averaged
#     together the way calibrate_ve_max()'s objective is
#
# f.aggregate_ve()'s "True VE (sim)"/"True VE (site-specific)" measures are
# noisy because their VE ratio is built from actual trial-sized event
# counts (a few hundred children x 12 months, repeated nsim times) — a
# small-sample ratio estimate. This function's ratio is instead built from
# n_draws Monte Carlo draws PER CELL (default 1e6), so its noise floor is
# roughly sqrt(p(1-p)/1e6) ~= 0.0005 per probability, several orders of
# magnitude tighter than trial-sized noise. Use this as a clean reference
# "designed true VE" to sanity-check or replace the noisier simulated
# benchmarks.
#
# Usage (after running f.param() to get the calibrated ve_max):
#   params <- f.param(ve_infection, ve_disease, ve_severity, ve_ct)
#   compute_true_ve_by_site(
#     score_params = params$sev_params$score_params,
#     gems_coefs   = params$sev_params$gems_params$coefs,
#     diarrhea_ir  = params$diarrhea_ir,
#     ve_max       = params$sev_params$ve_max,
#     ve_disease   = ve_disease
#   )
# ---------------------------------------------------------------
compute_true_ve_by_site <- function(score_params, gems_coefs, diarrhea_ir,
                                    ve_max, ve_disease, n_draws = 1e6) {

  cells <- score_params %>%
    left_join(diarrhea_ir, by = c("country_id", "agegrp"))

  cell_probs <- cells %>%
    rowwise() %>%
    mutate(
      score_res     = list(simulate_score_ve_cell(ve_max, mean_score, sd_score, n_draws)),
      p_unvax_score = score_res$p_unvax,
      p_vax_score   = score_res$p_vax,
      gems_res      = list(simulate_gems_ve_cell(ve_max, mean_score, sd_score,
                                                  gems_coefs, agegrp, country_id, n_draws)),
      p_unvax_gems  = gems_res$p_unvax,
      p_vax_gems    = gems_res$p_vax
    ) %>%
    ungroup() %>%
    select(country_id, agegrp, IR_shigella,
           p_unvax_score, p_vax_score, p_unvax_gems, p_vax_gems)

  # incidence-weighted conditional -> marginal VE, for whatever set of
  # country x agegrp cells is passed in (a single country's 2 cells for
  # per-site rows, or all 8 cells for the pooled "All Sites" row)
  summarise_ve <- function(df) {
    cond_ve_score <- 1 - weighted.mean(df$p_vax_score, df$IR_shigella) /
                          weighted.mean(df$p_unvax_score, df$IR_shigella)
    cond_ve_gems  <- 1 - weighted.mean(df$p_vax_gems, df$IR_shigella) /
                          weighted.mean(df$p_unvax_gems, df$IR_shigella)
    tibble(
      true_ve_score = 1 - (1 - ve_disease) * (1 - cond_ve_score),
      true_ve_gems  = 1 - (1 - ve_disease) * (1 - cond_ve_gems)
    )
  }

  per_site <- cell_probs %>%
    group_by(country_id) %>%
    group_modify(~ summarise_ve(.x)) %>%
    ungroup()

  pooled <- summarise_ve(cell_probs) %>%
    mutate(country_id = "ALL")

  bind_rows(per_site, pooled) %>%
    select(country_id, true_ve_score, true_ve_gems)
}

f.param <- function(ve_infection,
                    ve_disease,
                    ve_severity,
                    ve_ct) {
  
  # Diarrhea incidence: apply VE against disease to shigella IR
  diarrhea_ir_out <- diarrhea_ir %>% 
    mutate(vax_IR_shigella = IR_shigella * (1 - ve_disease))
  
  # Subclinical probability: apply VE against infection
  shigella_sub_out <- shigella_sub %>% 
    left_join(diarrhea_ir_out, by = c("country_id", "agegrp")) %>% 
    mutate(vax_predicted_prob = (predicted_prob + (IR_shigella - vax_IR_shigella)) * (1 - ve_infection))
  
  # Shigella severity: calibrate the linear severity-dependent shrinkage
  # parameter (ve_max) so the AVERAGE of the marginal VEs against
  # "score >= 6" MSD and GEMS MSD equals ve_severity (see calibrate_ve_max()
  # for why "marginal" — it compounds with ve_disease, since both MSD truth
  # definitions require a Shigella diarrhea episode AND severity). The score
  # draw itself is unaffected by vax status — shrinkage is applied post-draw
  # inside simulate_trial. Cached so that scenarios sharing the same
  # (ve_severity, ve_disease) pair get an identical ve_max.
  ve_max_calibrated <- get_calibrated_ve_max(
    ve_severity   = ve_severity,
    score_params  = sev_params$score_params,
    gems_coefs    = sev_params$gems_params$coefs,
    diarrhea_ir   = diarrhea_ir,
    ve_disease    = ve_disease
  )
  
  sev_params_out <- sev_params
  sev_params_out$ve_max <- ve_max_calibrated
  
  # Other severity: no vax effect assumed
  other_sev_params_out <- other_sev_params
  
  # Shigella pathogen quantity: vaccinated breakthrough cases have lower
  # quantity (i.e. weaker infections). No upper-bound workaround needed
  # here (unlike the old Ct-based approach, which needed an artificial
  # cap at 34.9 to avoid exceeding the Ct=35 detection ceiling) — quantity
  # has a natural lower bound at 0 enforced at simulation time instead.
  # shape_quantity (the Gamma shape parameter) is left unchanged: scaling
  # the mean while holding shape fixed already implies a proportionally
  # smaller variance for breakthrough cases (Var = mean^2 / shape).
  shigella_quantity_out <- shigella_quantity %>% 
    mutate(vax_mean_quantity = mean_quantity * (1 - ve_ct))
  
  return(list(
    diarrhea_ir          = diarrhea_ir_out,
    shigella_sub         = shigella_sub_out,
    sev_params           = sev_params_out,
    other_sev_params     = other_sev_params_out,
    shigella_quantity    = shigella_quantity_out,
    other_quantity       = other_quantity
  ))
}