## 2-arm trial power calculation, simulation-based

pacman::p_load(dplyr, tidyr, Rcpp)

sourceCpp(here::here("code", "functions", "power-calc-funcs.cpp"))

################################################################################
# Quick single-scenario check ----
################################################################################

p1           <- 0.014709091
p2           <- 0.008627000
nsim         <- 1000
alpha        <- 0.05
target_power <- 0.80

results <- find_sample_size_cpp(p1, p2, nsim, alpha, target_power)
print(results)

min_n1 <- results$n1[nrow(results)]
cat(sprintf("\nMinimum n per arm: %d  |  Total: %d\n", min_n1, min_n1 * 2))

################################################################################
# Loop through simulated risk data ----
################################################################################

obs.risk <- readRDS(here::here("clean data", "risk_results3.RDS")) %>%
  filter(measure == "Observed")

# extract p1/p2 per cell
risk_cells <- obs.risk %>%
  group_by(country_id, endpoint, severity, scenario) %>%
  summarise(
    p1 = mean[vax == 0], # placebo
    p2 = mean[vax == 1], # vaccine
    .groups = "drop"
  )

# run C++ power search for every cell (rowwise, no rbind accumulation)
sample_size_results <- risk_cells %>%
  rowwise() %>%
  mutate(sample_size = list(find_sample_size_cpp(p1, p2, nsim, alpha, target_power))) %>%
  unnest(sample_size) %>%
  filter(power >= target_power)

saveRDS(sample_size_results, here::here("clean data", "sample_sizes3.RDS"))

################################################################################
# Stratified (CMH) power calculation across sites ----
#
# risk_cells/sample_size_results already give three ways to summarise across
# sites:
#   - per-site: sample_size_results filtered to individual country_id values
#   - pooled ("All Sites"): a single unstratified 2x2 test using p1/p2 pooled
#     by summing events/person-time across sites (built upstream in
#     sim-analysis-functions.R) — assumes site is irrelevant to the analysis
#   - summed per-site totals (below): four separately-powered site trials,
#     enrollment added up — assumes sites are analyzed independently
#
# None of those is a stratified analysis. This section instead simulates a
# Cochran-Mantel-Haenszel test each replicate: draw outcomes at each site's
# own p1/p2, combine across sites into one stratified test statistic, and
# search for the total per-arm N (split across sites by site_weights) that
# hits target_power on THAT test. This is the right comparison if your
# primary analysis will adjust/stratify for site rather than pool sites
# naively or analyze them separately.
################################################################################

# site allocation weights: how the total per-arm N is split across sites.
# Defaults to equal allocation across the four sites, matching how the
# simulation itself assigns equal n_children per country (see
# simulate_trial() in 2-functions.cpp/2-functions.R). Renormalised internally
# in find_sample_size_stratified_cpp, so raw enrollment targets/ratios work
# too — replace with actual planned/observed site enrollment if unequal.
site_weights <- c(Bangladesh = 0.25, India = 0.25, Peru = 0.25, Pakistan = 0.25)

strat_cells <- risk_cells %>%
  filter(country_id != "All Sites") %>%
  group_by(endpoint, severity, scenario) %>%
  summarise(
    site_id = list(as.character(country_id)),
    p1      = list(p1),
    p2      = list(p2),
    .groups = "drop"
  )

# run stratified C++ power search for every cell (rowwise; p1/p2/site_id are
# length-K vectors within each row's group)
stratified_sample_size_results <- strat_cells %>%
  rowwise() %>%
  mutate(sample_size = list({
    w <- unname(site_weights[site_id])
    find_sample_size_stratified_cpp(p1, p2, w, nsim, alpha, target_power)
  })) %>%
  select(endpoint, severity, scenario, sample_size) %>%
  unnest(sample_size) %>%
  filter(power >= target_power)

saveRDS(stratified_sample_size_results,
        here::here("clean data", "sample_sizes3_stratified.RDS"))

################################################################################
# Total sample size across sites: three-way comparison ----
################################################################################

# pooled ("All Sites") total N — single unstratified trial combining all sites
pooled_n <- sample_size_results %>%
  filter(country_id == "All Sites") %>%
  select(endpoint, severity, scenario, ntotal) %>%
  rename(ntotal_pooled = ntotal)

# sum of per-site totals — four separately-powered site-level trials added up
summed_n <- sample_size_results %>%
  filter(country_id != "All Sites") %>%
  group_by(endpoint, severity, scenario) %>%
  summarise(ntotal_sum_of_sites = sum(ntotal), .groups = "drop")

# stratified (CMH) total N — one test adjusting for site
stratified_n <- stratified_sample_size_results %>%
  select(endpoint, severity, scenario, ntotal) %>%
  rename(ntotal_stratified = ntotal)

n_comparison <- pooled_n %>%
  full_join(summed_n,     by = c("endpoint", "severity", "scenario")) %>%
  full_join(stratified_n, by = c("endpoint", "severity", "scenario"))

saveRDS(n_comparison, here::here("clean data", "sample_sizes3_comparison.RDS"))
