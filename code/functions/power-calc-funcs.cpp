#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ---------------------------------------------------------------
// Native chi-square test (no R callback).
// Returns p-value from chi-square distribution with 1 df.
// Falls back to a continuity-corrected version (Yates) when any
// expected cell count is < 5, which is accurate enough for the
// sample sizes used in power calculations and avoids the overhead
// of calling back into R's fisher.test.
// ---------------------------------------------------------------
inline double chisq_pvalue(int s1, int f1, int s2, int f2) {
  double n = s1 + f1 + s2 + f2;
  double row1 = s1 + f1, row2 = s2 + f2;
  double col1 = s1 + s2, col2 = f1 + f2;

  // Expected counts
  double e11 = row1 * col1 / n;
  double e12 = row1 * col2 / n;
  double e21 = row2 * col1 / n;
  double e22 = row2 * col2 / n;

  double chi2;
  if (e11 < 5 || e12 < 5 || e21 < 5 || e22 < 5) {
    // Yates continuity correction
    double yates = n / 2.0;
    double d = std::abs(s1 * f2 - s2 * f1) - yates;
    if (d < 0) d = 0;
    chi2 = (n * d * d) / (row1 * row2 * col1 * col2);
  } else {
    chi2 = (s1 - e11) * (s1 - e11) / e11 +
           (f1 - e12) * (f1 - e12) / e12 +
           (s2 - e21) * (s2 - e21) / e21 +
           (f2 - e22) * (f2 - e22) / e22;
  }

  // p-value from chi-square(1): use R::pchisq (header-only, no callback)
  return R::pchisq(chi2, 1.0, false, false);
}

// [[Rcpp::export]]
double simulate_trial_cpp(int n1, double p1, double p2) {
  int s1 = static_cast<int>(R::rbinom(n1, p1));
  int s2 = static_cast<int>(R::rbinom(n1, p2));
  return chisq_pvalue(s1, n1 - s1, s2, n1 - s2);
}

// [[Rcpp::export]]
DataFrame find_sample_size_cpp(double p1, double p2,
                                int nsim         = 1000,
                                double alpha     = 0.05,
                                double target_power = 0.80) {
  int n1        = 100;
  int n1_step   = 50;

  // Pre-allocate result vectors (avoids repeated reallocation)
  std::vector<int>    n1_vec;
  std::vector<double> power_vec;
  std::vector<int>    ntotal_vec;
  n1_vec.reserve(50);
  power_vec.reserve(50);
  ntotal_vec.reserve(50);

  while (true) {
    int rejections = 0;
    for (int i = 0; i < nsim; i++) {
      int s1 = static_cast<int>(R::rbinom(n1, p1));
      int s2 = static_cast<int>(R::rbinom(n1, p2));
      if (chisq_pvalue(s1, n1 - s1, s2, n1 - s2) < alpha)
        rejections++;
    }

    double power_estimate = static_cast<double>(rejections) / nsim;
    n1_vec.push_back(n1);
    power_vec.push_back(power_estimate);
    ntotal_vec.push_back(n1 * 2);

    if (power_estimate >= target_power) break;
    n1 += n1_step;
  }

  return DataFrame::create(
    _["n1"]     = n1_vec,
    _["power"]  = power_vec,
    _["ntotal"] = ntotal_vec
  );
}

// ---------------------------------------------------------------
// Cochran-Mantel-Haenszel test for K stratified 2x2 tables (one per
// site). Combines evidence across strata into a single test of
// association, adjusting for site rather than either pooling all
// sites into one 2x2 table (which ignores stratification) or testing
// each site separately (which ignores shared information across
// sites). Uses the standard continuity-corrected MH chi-square
// statistic (1 df):
//
//   chi2_MH = (|sum(a_k - E[a_k])| - 0.5)^2 / sum(Var(a_k))
//
// where, for stratum k with cell counts a_k (vax events), b_k (vax
// non-events), c_k (placebo events), d_k (placebo non-events),
// n_k = a_k+b_k+c_k+d_k:
//   E[a_k]   = (a_k+b_k)(a_k+c_k) / n_k
//   Var(a_k) = (a_k+b_k)(c_k+d_k)(a_k+c_k)(b_k+d_k) / (n_k^2 (n_k-1))
//
// Strata with n_k < 2 are skipped (Var(a_k) undefined) rather than
// discarding the whole test.
// ---------------------------------------------------------------
inline double cmh_pvalue(const std::vector<int>& a, const std::vector<int>& b,
                          const std::vector<int>& c, const std::vector<int>& d) {
  double sum_diff = 0.0;
  double sum_var  = 0.0;
  int K = a.size();
  for (int k = 0; k < K; k++) {
    double ak = a[k], bk = b[k], ck = c[k], dk = d[k];
    double n  = ak + bk + ck + dk;
    if (n < 2) continue;
    double row1 = ak + bk, row2 = ck + dk;
    double col1 = ak + ck, col2 = bk + dk;
    double e_a   = row1 * col1 / n;
    double var_a = (row1 * row2 * col1 * col2) / (n * n * (n - 1));
    sum_diff += (ak - e_a);
    sum_var  += var_a;
  }
  if (sum_var <= 0) return 1.0;
  double d_corr = std::abs(sum_diff) - 0.5;  // continuity correction
  if (d_corr < 0) d_corr = 0;
  double chi2 = (d_corr * d_corr) / sum_var;
  return R::pchisq(chi2, 1.0, false, false);
}

// [[Rcpp::export]]
double simulate_trial_stratified_cpp(IntegerVector n1_per_site,
                                      NumericVector p1_per_site,
                                      NumericVector p2_per_site) {
  int K = n1_per_site.size();
  std::vector<int> a(K), b(K), c(K), d(K);
  for (int k = 0; k < K; k++) {
    int nk  = n1_per_site[k];
    int s1k = static_cast<int>(R::rbinom(nk, p1_per_site[k]));
    int s2k = static_cast<int>(R::rbinom(nk, p2_per_site[k]));
    a[k] = s1k;      b[k] = nk - s1k;
    c[k] = s2k;      d[k] = nk - s2k;
  }
  return cmh_pvalue(a, b, c, d);
}

// ---------------------------------------------------------------
// find_sample_size_stratified_cpp
//
// Simulation-based sample size search for a stratified (CMH) test
// across K sites, each with its own p1/p2 (site-specific placebo/
// vaccine risk). site_weights gives the proportion of the TOTAL
// per-arm sample allocated to each site (renormalised to sum to 1
// internally, so raw enrollment counts or ratios both work). Total
// per-arm n is increased in steps of n1_step, split across sites by
// site_weights (each site floored at a minimum of 1 per arm), until
// the simulated CMH test hits target_power.
//
// Returns the same n1/power/ntotal shape as find_sample_size_cpp, but
// n1/ntotal here are TOTALS summed across all sites (after per-site
// rounding), not per-site counts.
// ---------------------------------------------------------------
// [[Rcpp::export]]
DataFrame find_sample_size_stratified_cpp(NumericVector p1_per_site,
                                           NumericVector p2_per_site,
                                           NumericVector site_weights,
                                           int nsim            = 1000,
                                           double alpha        = 0.05,
                                           double target_power = 0.80) {
  int K = p1_per_site.size();

  double wsum = 0.0;
  for (int k = 0; k < K; k++) wsum += site_weights[k];
  std::vector<double> w(K);
  for (int k = 0; k < K; k++) w[k] = site_weights[k] / wsum;

  int n1_total_target = 100;  // total per arm across all sites, combined
  int n1_step          = 50;

  std::vector<int>    n1_vec;
  std::vector<double> power_vec;
  std::vector<int>    ntotal_vec;
  n1_vec.reserve(50);
  power_vec.reserve(50);
  ntotal_vec.reserve(50);

  while (true) {
    // allocate per-site per-arm n by weight, floored at 1
    std::vector<int> n1_site(K);
    int actual_n1_total = 0;
    for (int k = 0; k < K; k++) {
      n1_site[k] = std::max(1, static_cast<int>(std::round(n1_total_target * w[k])));
      actual_n1_total += n1_site[k];
    }

    int rejections = 0;
    for (int i = 0; i < nsim; i++) {
      std::vector<int> a(K), b(K), c(K), d(K);
      for (int k = 0; k < K; k++) {
        int nk  = n1_site[k];
        int s1k = static_cast<int>(R::rbinom(nk, p1_per_site[k]));
        int s2k = static_cast<int>(R::rbinom(nk, p2_per_site[k]));
        a[k] = s1k;      b[k] = nk - s1k;
        c[k] = s2k;      d[k] = nk - s2k;
      }
      if (cmh_pvalue(a, b, c, d) < alpha) rejections++;
    }

    double power_estimate = static_cast<double>(rejections) / nsim;
    n1_vec.push_back(actual_n1_total);
    power_vec.push_back(power_estimate);
    ntotal_vec.push_back(actual_n1_total * 2);

    if (power_estimate >= target_power) break;
    n1_total_target += n1_step;
  }

  return DataFrame::create(
    _["n1"]     = n1_vec,
    _["power"]  = power_vec,
    _["ntotal"] = ntotal_vec
  );
}
