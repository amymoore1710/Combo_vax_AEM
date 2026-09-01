# libraries ----
pacman::p_load(readr,
               dplyr,
               tidyr,
               zoo)

# data -----
tac <- read_csv(here::here("maled data/Tac_Sep2018.csv"))
load(here::here("maled data","shigella_specimen.RData")) 
load(here::here("maled data","new_shigella_cause.RData"))

# countries and age groups for simulation ----
countries <- c("BG","PE","PK","IN")
agegroups <- c("12-17 months","18-24 months")

# formatted data to estimate incidence of shigella & other diarrhea ----
inc.df <- tac %>% 
  select(country_id,pid,sid,date,agedays,stooltype,shigella_eiec,shigella_micro) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>% 
  mutate(diarrhea = ifelse(stooltype == "D1",1,0)) %>% 
  left_join(shig.cause,by = "sid") %>% 
  mutate(month_year = as.yearmon(date),
         shig.diar = ifelse(grepl("Shigella",diarrhea.cause) | (diarrhea == 1 & shigella_micro == 1),1,0),
         shig.diar = ifelse(is.na(shig.diar) & diarrhea == 1,0,shig.diar)) %>% 
  mutate(other.diar = ifelse(diarrhea == 1 & shig.diar == 0,1,0)) %>% 
  group_by(country_id,pid,month_year) %>% 
  mutate(agemonth = mean(agedays/30.5)) %>% 
  mutate(agegrp = case_when(agedays<6*30.5 ~ "0-5 months",
                            agedays>=6*30.5 & agedays<12*30.5 ~ "6-11 months",
                            agedays>=12*30.5 & agedays<18*30.5 ~ "12-17 months",
                            agedays>=18*30.5 ~ "18-24 months")) %>% 
  group_by(country_id,agegrp,pid) %>% 
  summarise(shig.diar = sum(shig.diar),
            other.diar = sum(other.diar)) %>% 
  mutate(obs_months = 6) %>% 
  mutate(monthly.shigella = shig.diar/obs_months,
         monthly.other = other.diar/obs_months)

# incidence poisson models ----

# Shigella diarrhea
m1 <- glm(shig.diar ~ country_id * agegrp + offset(log(obs_months)), 
          data = inc.df, 
          family = poisson())

summary(m1)

inc.df$pred.shig <- exp(predict(m1, type = "link"))

# Other diarrhea
m2 <- glm(other.diar ~ country_id * agegrp + offset(log(obs_months)), 
          data = inc.df, 
          family = poisson())

summary(m2)

inc.df$pred.other <- exp(predict(m2, type = "link"))

# Predicted monthly counts of shigella & other diarrhea ----
est.df <- expand.grid(
  country_id = unique(inc.df$country_id),
  agegrp = unique(inc.df$agegrp),
  obs_months = 1  # Set to 1 for incidence rate calculation
)

est.df$log_count_shigella <- predict(m1, newdata = est.df, type = "link")
est.df$count_shigella <- exp(est.df$log_count_shigella)
est.df <- est.df %>%
  mutate(IR_shigella = count_shigella / obs_months)  # obs_months is 1 here

est.df$log_count_other <- predict(m2, newdata = est.df, type = "link")
est.df$count_other <- exp(est.df$log_count_other)
est.df <- est.df %>%
  mutate(IR_other = count_other / obs_months)  # obs_months is 1 here

saveRDS(est.df %>%
          filter(country_id %in% countries,
                 agegrp %in% agegroups) %>%
          select(country_id,agegrp,IR_other,IR_shigella),here::here("sim param data","diarrhea_incidence.RDS"))

# Probability of subclinical shigella ----
shig.sub <- tac %>% 
  select(country_id,pid,sid,agedays,stooltype,shigella_eiec,shigella_micro) %>% 
  left_join(shig.cause,by = "sid") %>% 
  mutate(agemonths = floor(agedays/30.5),
         agegrp = case_when(agemonths<6 ~ "0-5 months",
                            agemonths>=6 & agemonths<12 ~ "6-11 months",
                            agemonths>=12 & agemonths<18 ~ "12-17 months",
                            agemonths>=18 ~ "18-24 months"),
         shigella.diarrhea = ifelse(grepl("Shigella",diarrhea.cause) | (stooltype == "D1" & shigella_micro == 1),1,0),
         shigella.subclin = ifelse(shigella_eiec < 35 & shigella.diarrhea == 0, 1, 0))

shig.sub2 <- shig.sub %>% 
  group_by(country_id,pid,agegrp,agemonths) %>% 
  summarise(shigella.diarrhea = sum(shigella.diarrhea),
            shigella.subclin = sum(shigella.subclin)) %>% 
  mutate(shigella.diarrhea = ifelse(shigella.diarrhea>1,1,0),
         shigella.subclin = ifelse(shigella.subclin>1,1,0)) %>% 
  ungroup() %>% 
  arrange(pid,agemonths) %>% 
  group_by(pid) %>% 
  mutate(prev.diarrhea = ifelse(lag(shigella.diarrhea, default = 0) == 1 |
                                  shigella.diarrhea == 1,1,0) ) 

m3 <- glm(shigella.subclin ~ country_id * agegrp + prev.diarrhea,
             data = shig.sub2, family = binomial)
summary(m3)

# Predicted monthly probability for subclinical shigella ----
est.df2 <- expand.grid(
  country_id = unique(shig.sub2$country_id),
  agegrp = unique(shig.sub2$agegrp),
  prev.diarrhea = c(0, 1)
)

est.df2$predicted_log_odds <- predict(m3, est.df2, type = "link")
est.df2$predicted_prob <- exp(est.df2$predicted_log_odds) / (1 + exp(est.df2$predicted_log_odds))

saveRDS(est.df2 %>% 
          filter(country_id %in% countries,
                 agegrp %in% agegroups) %>% 
          select(country_id,agegrp,prev.diarrhea,predicted_prob),
        here::here("sim param data","shigella_subclinical_prob.RDS"))

# NOTE: subclinical culture probability is now handled by the unified
# shigella_micro culture model below (see "Shigella severity and culture
# probability" section), which is fit on ALL PCR-positive specimens
# (shigella_eiec < 35), diarrhea or not, using Ct and shig.diar as
# predictors. No separate subclinical-only model is needed.

# Probability of other pathogen detections ----
other.det <- tac %>% 
  select(country_id,pid,sid,stooltype,agedays,adenovirus_40_41,
         astrovirus,campylobacter_jejuni_coli,cryptosporidium,
         ST_ETEC,norovirus_gii,rotavirus,sapovirus,tEPEC) %>% 
  pivot_longer(cols = adenovirus_40_41:tEPEC,
               names_to = "pathogen",
               values_to = "ct") %>% 
  na.omit(.) %>% 
  group_by(sid) %>% 
  mutate(det = ifelse(any(ct <35),1,0)) %>% 
  mutate(other.diar = ifelse(stooltype == "D1", 1, 0)) |> 
  select(country_id,pid,sid,agedays,other.diar,det) %>% 
  distinct() %>% 
  mutate(agemonths = floor(agedays/30.5),
         agegrp = case_when(agemonths<6 ~ "0-5 months",
                            agemonths>=6 & agemonths<12 ~ "6-11 months",
                            agemonths>=12 & agemonths<18 ~ "12-17 months",
                            agemonths>=18 ~ "18-24 months")) %>% 
  group_by(country_id,pid,other.diar,agegrp,agemonths) %>% 
  summarise(det = max(det))

m4 <- glm(det ~ country_id + agegrp + other.diar,
          data = other.det, family = binomial)
summary(m4)

# Predicted monthly probability for other detections ----
est.df3 <- expand.grid(
  country_id = unique(other.det$country_id),
  agegrp = unique(other.det$agegrp),
  other.diar = unique(other.det$other.diar)
)

est.df3$predicted_log_odds <- predict(m4, est.df3, type = "link")
est.df3$predicted_prob <- exp(est.df3$predicted_log_odds) / (1 + exp(est.df3$predicted_log_odds))

saveRDS(est.df3 %>% 
          filter(country_id %in% countries,
                 agegrp %in% agegroups) %>% 
          select(country_id,agegrp,other.diar,predicted_prob),
        here::here("sim param data","other_subclinical_prob.RDS"))

# Shigella severity and culture probability ----
sev.df <- tac %>% 
  select(country_id,pid,sid,agedays,stooltype,score,gemsdef,shigella_eiec,shigella_micro) %>% 
  mutate(diarrhea = ifelse(stooltype == "D1",1,0)) %>% 
  left_join(shig.cause,by = "sid") %>% 
  mutate(shig.diar = ifelse(grepl("Shigella",diarrhea.cause) | (diarrhea == 1 & shigella_micro == 1),1,0),
         other.diar = ifelse(diarrhea == 1 & shig.diar == 0,1,0),
         agegrp = case_when(agedays<6*30.5 ~ "0-5 months",
                            agedays>=6*30.5 & agedays<12*30.5 ~ "6-11 months",
                            agedays>=12*30.5 & agedays<18*30.5 ~ "12-17 months",
                            agedays>=18*30.5 ~ "18-24 months")) |> 
  mutate(agegrp = relevel(factor(agegrp), ref = "12-17 months"),
          country_id = relevel(factor(country_id), ref = "BG")) 

# 1. Score distribution: observed mean/SD by country x agegrp
score_params <- sev.df %>%
  filter(shig.diar == 1,
         country_id %in% countries,
        agegrp %in% agegroups) %>%
  group_by(country_id, agegrp) %>%
  summarise(
    mean_score = mean(score),
    sd_score   = sd(score),
    n          = n(),
    .groups    = "drop"
  )

# 2. P(gemsdef): logistic regression with continuous score
glm_gems <- glm(gemsdef ~ agegrp + country_id + score,
                data = sev.df %>% filter(shig.diar == 1),
                family = binomial)
# also tried fitting with score as factor, with restricted cubic spline, and as log(score) and linear was the best fit

gems_params <- list(
  coefs = coef(glm_gems),
  vcov  = vcov(glm_gems)
)

# 3. P(culture+): unified model across ALL PCR-positive Shigella specimens
# (diarrhea or subclinical), using Ct value and Shigella-attributable
# diarrhea status as predictors. Replaces the previous two-model approach
# (gemsdef-based for clinical cases, agegrp-only for subclinical cases).

# check cell sizes before fitting — shig.diar=1 vs 0 may be very unbalanced
# within some country/agegrp strata
sev.df %>%
  filter(shigella_eiec < 35) %>%
  group_by(country_id, agegrp, shig.diar) %>%
  summarise(n = n(), n_culture = sum(shigella_micro, na.rm = TRUE), .groups = "drop")

glm_culture <- glm(shigella_micro ~ agegrp + country_id + shigella_eiec + shig.diar,
                data = sev.df %>% filter(shigella_eiec < 35),
                family = binomial)
summary(glm_culture)

culture_params <- list(
  coefs = coef(glm_culture),
  vcov  = vcov(glm_culture)
)

# save severity and culture parameters
sev_params <- list(
  score_params   = score_params,
  gems_params    = gems_params,
  culture_params = culture_params,
  ref_agegrp     = levels(sev.df$agegrp)[1],
  ref_country    = levels(sev.df$country_id)[1]
)

saveRDS(sev_params, here::here("sim param data","shigella_severity.RDS"))

# Other diarrhea severity ----

# 1. Score distribution: observed mean/SD by country x agegrp
other_score_params <- sev.df %>%
  filter(other.diar == 1,
         country_id %in% countries,
        agegrp %in% agegroups) %>%
  group_by(country_id, agegrp) %>%
  summarise(
    mean_score = mean(score),
    sd_score   = sd(score),
    n          = n(),
    .groups    = "drop"
  )

# 2. P(gemsdef): logistic regression with continuous score
glm_gems_other <- glm(gemsdef ~ agegrp + country_id + score,
                data = sev.df %>% filter(other.diar == 1),
                family = binomial)
# also tried fitting with score as factor, with restricted cubic spline, and as log(score) and linear was the best fit

other_gems_params <- list(
  coefs = coef(glm_gems_other),
  vcov  = vcov(glm_gems_other)
)

# save severity and culture parameters
other_sev_params <- list(
  score_params   = other_score_params,
  gems_params    = other_gems_params,
  ref_agegrp     = levels(sev.df$agegrp)[1],
  ref_country    = levels(sev.df$country_id)[1]
)

saveRDS(other_sev_params, here::here("sim param data","other_severity.RDS"))

# Pathogen quantity distribution: shared model-fitting helper ----
# Both Shigella and other-pathogen quantity are now fit with the SAME model
# structure: a Gamma GLM with a log link, conditional on country, age group,
# and a 3-level "severity" category (Subclinical / Mild / Severe), where:
#   - "Subclinical" = detected outside of a diarrhea episode (for Shigella:
#     a PCR-positive month with no Shigella-attributable diarrhea; for other
#     pathogens: detected on a month where "other diarrhea" is not occurring)
#   - "Mild" / "Severe" = detected during a diarrhea episode of that severity
#
# Quantity is modeled directly (rather than Ct then converted) since it is
# strictly positive and right-skewed, which a Gamma GLM accommodates without
# needing an artificial upper truncation near the Ct=35 detection ceiling
# (unlike modeling Ct directly, where quantities near the ceiling are hard
# to truncate sensibly).
# Conversion used elsewhere: Ct = 35 - 3.322 * quantity
#                       <=>  quantity = (35 - Ct) / 3.322
#
# No random effects are included: with country x agegrp x severity already
# as fixed-effect predictors, a random intercept on country:agegrp added
# little and isn't needed for the simulation's purposes.
#
# The Gamma dispersion parameter (phi) is assumed constant across cells
# (mirroring the previous approach of assuming a constant residual SD across
# groups) and is converted to a shape parameter (shape = 1/phi) for use when
# drawing from rgamma() at simulation time: with mean mu and shape k, the
# implied scale is mu/k, so Var(quantity) = mu^2/k scales naturally with the
# mean rather than being fixed in absolute terms — appropriate for a
# log-link Gamma model.
fit_quantity_gamma <- function(data) {
  model <- glm(quantity ~ country_id + agegrp + severity,
               data = data, family = Gamma(link = "log"))

  disp  <- summary(model)$dispersion  # phi: Var(Y) = phi * mu^2
  shape <- 1 / disp                   # constant shape across all cells

  pred_df <- data %>%
    distinct(country_id, agegrp, severity) %>%
    mutate(mean_quantity  = predict(model, newdata = ., type = "response"),
           shape_quantity = shape)

  list(model = model, pred_df = pred_df, dispersion = disp, shape = shape)
}

# Diagnostic: print the fitted mean_quantity/shape_quantity table restricted
# to ONLY the country x agegrp combinations actually used in the simulation
# (countries, agegroups — defined above), crossed with all 3 severity
# levels, and flag any combination missing from the fit (which would cause
# a lookup failure inside assign_shig_quantity()/assign_other_quantity() at
# simulation time, since those functions only know about cells present in
# data — see note above on the same issue with the previous lmer() approach).
check_quantity_coverage <- function(quantity_dist, label) {
  expected_grid <- expand.grid(
    country_id = countries,
    agegrp     = agegroups,
    severity   = c("Subclinical", "Mild", "Severe"),
    stringsAsFactors = FALSE
  )

  sim_table <- quantity_dist %>%
    filter(country_id %in% countries,
           agegrp %in% agegroups) %>%
    right_join(expected_grid, by = c("country_id", "agegrp", "severity")) %>%
    arrange(country_id, agegrp,
            factor(severity, levels = c("Subclinical", "Mild", "Severe")))

  missing <- sim_table %>% filter(is.na(mean_quantity))
  if (nrow(missing) > 0) {
    warning(sprintf(
      "%s: %d country x agegrp x severity cell(s) used in the simulation have NO fitted value (will cause a lookup failure at simulation time):\n%s",
      label, nrow(missing),
      paste(capture.output(print(missing %>% select(country_id, agegrp, severity))), collapse = "\n")
    ))
  } else {
    message(sprintf("%s: all %d simulated country x agegrp x severity cells have a fitted value.",
                    label, nrow(expected_grid)))
  }

  print(as.data.frame(sim_table))
  sim_table
}

# Shigella pathogen quantity distribution ----
# Since this filter requires shigella_eiec < 35, quantity > 0 always holds.

ct.df <- sev.df %>% 
  filter(shigella_eiec < 35) %>% 
  mutate(
    quantity = (35 - shigella_eiec) / 3.322,
    severity = case_when(
      shig.diar == 1 & (score >= 6 | gemsdef == 1) ~ "Severe",
      shig.diar == 1  ~ "Mild",
      TRUE ~ "Subclinical"
    )
  )

shig_qty_fit <- fit_quantity_gamma(ct.df)
summary(shig_qty_fit$model)

# Merge with cell sample sizes
quantity.dist <- shig_qty_fit$pred_df %>%
  left_join(ct.df %>% count(country_id, agegrp, severity), by = c("country_id", "agegrp", "severity"))

# Display results: restricted to the country x agegrp combinations actually
# used in the simulation, and flag any missing cell. (Diagnostic only —
# does not affect what gets saved below.)
invisible(check_quantity_coverage(quantity.dist, "Shigella quantity"))

saveRDS(quantity.dist %>%
          filter(country_id %in% countries,
                 agegrp %in% agegroups) %>%
          select(country_id,agegrp,severity,mean_quantity,shape_quantity,n),
        here::here("sim param data","shigella_quantity.RDS"))

# Other pathogen quantity distribution ----
# Same model structure as Shigella above: Gamma GLM (log link) on quantity,
# conditional on country, agegrp, and severity (Subclinical / Mild / Severe).

pathogens <- c("adenovirus_40_41", "astrovirus", "campylobacter_jejuni_coli", 
               "cryptosporidium", "ST_ETEC", "norovirus_gii", "rotavirus", 
               "sapovirus", "tEPEC")

other.quantity.dist <- tac %>% 
  mutate(severity = case_when(
      stooltype == "M1" ~ "Subclinical",
      (score >= 6 | gemsdef == 1) ~ "Severe",
      stooltype == "D1" ~ "Mild"
    )) |> 
  select(country_id, pid, sid, agedays, severity, all_of(pathogens)) %>% 
  pivot_longer(cols = all_of(pathogens), 
               names_to = "pathogen", 
               values_to = "ct") %>% 
  na.omit() %>%  # Remove missing Ct values
  filter(ct < 35) %>%  # Ct >= 35 denotes non-detection; quantity must be > 0 for the Gamma GLM
  mutate(quantity = (35 - ct) / 3.322) %>% 
  group_by(sid) %>% 
  filter(quantity == max(quantity)) %>%  # Keep the highest-quantity pathogen per sample (equivalent to lowest Ct)
  select(-pathogen, -ct) %>% 
  distinct() %>% 
  ungroup() %>% 
  filter(!is.na(severity)) %>%   # drop samples that don't fall into a category
  mutate(
    agemonths = floor(agedays / 30.5),
    agegrp = case_when(
      agemonths < 6 ~ "0-5 months",
      agemonths >= 6 & agemonths < 12 ~ "6-11 months",
      agemonths >= 12 & agemonths < 18 ~ "12-17 months",
      agemonths >= 18 ~ "18-24 months"
    )
  )

other_qty_fit <- fit_quantity_gamma(other.quantity.dist)
summary(other_qty_fit$model)

quantity_estimates <- other_qty_fit$pred_df %>%
  left_join(other.quantity.dist %>% count(country_id, agegrp, severity),
            by = c("country_id", "agegrp", "severity"))

# View results: restricted to the country x agegrp combinations actually
# used in the simulation, and flag any missing cell. (Diagnostic only —
# does not affect what gets saved below.)
invisible(check_quantity_coverage(quantity_estimates, "Other pathogen quantity"))

saveRDS(quantity_estimates %>%
          filter(country_id %in% countries,
                 agegrp %in% agegroups) %>%
          select(country_id,agegrp,severity,mean_quantity,shape_quantity,n),
        here::here("sim param data","other_quantity.RDS"))