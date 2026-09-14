#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------
// Lookup key: country x agegrp (used by most parameter tables)
// ---------------------------------------------------------------
struct CAKey {
  std::string country;
  std::string agegrp;
  bool operator==(const CAKey& o) const {
    return country == o.country && agegrp == o.agegrp;
  }
};

struct CAKeyHash {
  std::size_t operator()(const CAKey& k) const {
    return std::hash<std::string>()(k.country) ^
          (std::hash<std::string>()(k.agegrp) << 1);
  }
};

// ---------------------------------------------------------------
// Pre-built lookup structs (populated once, used millions of times)
// ---------------------------------------------------------------

struct IRParams {
  double ir_shigella;
  double vax_ir_shigella;
  double ir_ETEC;
  double vax_ir_ETEC;
  double ir_other;
};

struct ScoreParams {
  double mean_score;
  double sd_score;
};

struct ShigSubclinParams {
  double pred_prob_0;      // prev_diarrhea == 0
  double pred_prob_1;      // prev_diarrhea == 1
  double vax_pred_prob_0;
  double vax_pred_prob_1;
};

struct ETECSubclinParams {
  double pred_prob_0;      // prev_diarrhea == 0
  double pred_prob_1;      // prev_diarrhea == 1
  double vax_pred_prob_0;
  double vax_pred_prob_1;
};

struct OtherSubParams {
  double pred_prob_0;      // other_diarrhea == 0
  double pred_prob_1;      // other_diarrhea == 1
};

// Pathogen quantity lookup key: country x agegrp x severity. Shared by both
// the Shigella and other-pathogen quantity tables, since both are now fit
// with the same model structure (Gamma GLM, log link, on country x agegrp x
// severity, where severity in {Subclinical, Mild, Severe}; no random
// effects — see 0-estimate-params.R).
struct QtyKey {
  std::string country;
  std::string agegrp;
  std::string severity;
  bool operator==(const QtyKey& o) const {
    return country == o.country && agegrp == o.agegrp && severity == o.severity;
  }
};
struct QtyKeyHash {
  std::size_t operator()(const QtyKey& k) const {
    size_t h = std::hash<std::string>()(k.country);
    h ^= std::hash<std::string>()(k.agegrp)   << 1;
    h ^= std::hash<std::string>()(k.severity) << 2;
    return h;
  }
};
struct ShigQtyParams {
  double mean_quantity;
  double vax_mean_quantity;
  double shape_quantity;  // Gamma shape parameter (1/dispersion), constant across cells
};

struct ETECQtyParams {
  double mean_quantity;
  double vax_mean_quantity;
  double shape_quantity;  // Gamma shape parameter (1/dispersion), constant across cells
};

struct OtherQtyParams {
  double mean_quantity;
  double shape_quantity;  // Gamma shape parameter (1/dispersion), constant across cells
};

// ---------------------------------------------------------------
// Build lookup maps from DataFrames (called once per simulate_trial)
// ---------------------------------------------------------------

std::unordered_map<CAKey, IRParams, CAKeyHash>
build_ir_map(DataFrame data) {
  std::unordered_map<CAKey, IRParams, CAKeyHash> m;
  CharacterVector country_ids   = data["country_id"];
  CharacterVector agegrps       = data["agegrp"];
  NumericVector ir_shig         = data["IR_shigella"];
  NumericVector ir_ETEC         = data["IR_ETEC"];
  NumericVector ir_other        = data["IR_other"];
  NumericVector vax_ir_shig     = data["vax_IR_shigella"];
  NumericVector vax_ir_ETEC     = data["vax_IR_ETEC"];
  for (int i = 0; i < country_ids.size(); i++) {
    CAKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i])};
    m[k] = {ir_shig[i], vax_ir_shig[i], ir_ETEC[i], vax_ir_ETEC[i], ir_other[i]};
  }
  return m;
}

std::unordered_map<CAKey, ScoreParams, CAKeyHash>
build_score_map(DataFrame data) {
  std::unordered_map<CAKey, ScoreParams, CAKeyHash> m;
  CharacterVector country_ids  = data["country_id"];
  CharacterVector agegrps      = data["agegrp"];
  NumericVector mean_score     = data["mean_score"];
  NumericVector sd_score       = data["sd_score"];
  for (int i = 0; i < country_ids.size(); i++) {
    CAKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i])};
    m[k] = {mean_score[i], sd_score[i]};
  }
  return m;
}

std::unordered_map<CAKey, ShigSubclinParams, CAKeyHash>
build_shig_subclin_map(DataFrame data) {
  std::unordered_map<CAKey, ShigSubclinParams, CAKeyHash> m;
  CharacterVector country_ids  = data["country_id"];
  CharacterVector agegrps      = data["agegrp"];
  IntegerVector prev_diar      = data["prev.diarrhea"];
  NumericVector pred_prob      = data["predicted_prob"];
  NumericVector vax_pred_prob  = data["vax_predicted_prob"];
  // initialise all entries so both prev_diarrhea values land in the same struct
  for (int i = 0; i < country_ids.size(); i++) {
    CAKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i])};
    auto& entry = m[k];
    if (prev_diar[i] == 0) {
      entry.pred_prob_0     = pred_prob[i];
      entry.vax_pred_prob_0 = vax_pred_prob[i];
    } else {
      entry.pred_prob_1     = pred_prob[i];
      entry.vax_pred_prob_1 = vax_pred_prob[i];
    }
  }
  return m;
}

std::unordered_map<CAKey, ETECSubclinParams, CAKeyHash>
  build_ETEC_subclin_map(DataFrame data) {
    std::unordered_map<CAKey, ETECSubclinParams, CAKeyHash> m;
    CharacterVector country_ids  = data["country_id"];
    CharacterVector agegrps      = data["agegrp"];
    IntegerVector prev_diar      = data["prev.diarrhea"];
    NumericVector pred_prob      = data["predicted_prob"];
    NumericVector vax_pred_prob  = data["vax_predicted_prob"];
    // initialise all entries so both prev_diarrhea values land in the same struct
    for (int i = 0; i < country_ids.size(); i++) {
      CAKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i])};
      auto& entry = m[k];
      if (prev_diar[i] == 0) {
        entry.pred_prob_0     = pred_prob[i];
        entry.vax_pred_prob_0 = vax_pred_prob[i];
      } else {
        entry.pred_prob_1     = pred_prob[i];
        entry.vax_pred_prob_1 = vax_pred_prob[i];
      }
    }
    return m;
  }

std::unordered_map<CAKey, OtherSubParams, CAKeyHash>
build_other_sub_map(DataFrame data) {
  std::unordered_map<CAKey, OtherSubParams, CAKeyHash> m;
  CharacterVector country_ids = data["country_id"];
  CharacterVector agegrps     = data["agegrp"];
  IntegerVector other_diar    = data["other.diar"];
  NumericVector pred_prob     = data["predicted_prob"];
  // initialise all entries so both other_diarrhea values land in the same struct
  for (int i = 0; i < country_ids.size(); i++) {
    CAKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i])};
    auto& entry = m[k];
    if (other_diar[i] == 0) {
      entry.pred_prob_0 = pred_prob[i];
    } else {
      entry.pred_prob_1 = pred_prob[i];
    }
  }
  return m;
}

std::unordered_map<QtyKey, ShigQtyParams, QtyKeyHash>
build_shig_quantity_map(DataFrame data) {
  std::unordered_map<QtyKey, ShigQtyParams, QtyKeyHash> m;
  CharacterVector country_ids    = data["country_id"];
  CharacterVector agegrps        = data["agegrp"];
  CharacterVector severities     = data["severity"];
  NumericVector mean_quantity      = data["mean_quantity"];
  NumericVector shape_quantity     = data["shape_quantity"];
  NumericVector vax_mean_quantity  = data["vax_mean_quantity"];
  for (int i = 0; i < country_ids.size(); i++) {
    QtyKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i]),
            as<std::string>(severities[i])};
    m[k] = {mean_quantity[i], vax_mean_quantity[i], shape_quantity[i]};
  }
  return m;
}


std::unordered_map<QtyKey, ETECQtyParams, QtyKeyHash>
  build_ETEC_quantity_map(DataFrame data) {
    std::unordered_map<QtyKey, ETECQtyParams, QtyKeyHash> m;
    CharacterVector country_ids    = data["country_id"];
    CharacterVector agegrps        = data["agegrp"];
    CharacterVector severities     = data["severity"];
    NumericVector mean_quantity      = data["mean_quantity"];
    NumericVector shape_quantity     = data["shape_quantity"];
    NumericVector vax_mean_quantity  = data["vax_mean_quantity"];
    for (int i = 0; i < country_ids.size(); i++) {
      QtyKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i]),
               as<std::string>(severities[i])};
      m[k] = {mean_quantity[i], vax_mean_quantity[i], shape_quantity[i]};
    }
    return m;
  }

std::unordered_map<QtyKey, OtherQtyParams, QtyKeyHash>
build_other_quantity_map(DataFrame data) {
  std::unordered_map<QtyKey, OtherQtyParams, QtyKeyHash> m;
  CharacterVector country_ids = data["country_id"];
  CharacterVector agegrps     = data["agegrp"];
  CharacterVector severities  = data["severity"];
  NumericVector mean_quantity  = data["mean_quantity"];
  NumericVector shape_quantity = data["shape_quantity"];
  for (int i = 0; i < country_ids.size(); i++) {
    QtyKey k{as<std::string>(country_ids[i]), as<std::string>(agegrps[i]),
            as<std::string>(severities[i])};
    m[k] = {mean_quantity[i], shape_quantity[i]};
  }
  return m;
}

// ---------------------------------------------------------------
// Determine whether diarrhea occurs in a given month
// [[Rcpp::export]]
int diarrhea_event(double lambda) {
  return (R::rexp(1.0 / lambda) <= 1.0) ? 1 : 0;
}

// ---------------------------------------------------------------
// Apply severity-dependent VE shrinkage to a drawn score for vaccinated
// children. Linear ramp: no effect at score = 1, maximal effect (ve_max)
// at score = 12. ve_max is calibrated in R (calibrate_ve_max() in
// 1-parameters.R) so the resulting VE against GEMS MSD matches the target.
// Mirrors shrink_score() in 1-parameters.R / 2-functions.R exactly.
// ---------------------------------------------------------------
inline double shrink_score(double score, double ve_max) {
  double ve_at_score = ve_max * (score - 1.0) / 11.0;
  double shrunk       = score * (1.0 - ve_at_score);
  return std::round(std::min(std::max(shrunk, 1.0), 12.0));
}

// ---------------------------------------------------------------
// Convert simulated pathogen quantity to Ct for use as a predictor in the
// culture model (which was fit on Ct directly: shigella_micro ~ ... +
// shigella_eiec + shig.diar). Quantity is simulated rather than Ct itself
// so that no artificial truncation near the Ct=35 detection ceiling is
// needed — quantity only needs a natural lower bound at 0.
// Conversion: Ct = 35 - 3.322 * quantity  (inverse of quantity = (35-Ct)/3.322
// used when estimating shigella_quantity.RDS in 0-estimate-params.R)
// ---------------------------------------------------------------
// Upper bound on simulated pathogen quantity: a Ct of 0 is the
// (theoretical) detection ceiling, so quantity should not exceed
// (35 - 0) / 3.322. The Gamma quantity models have no natural upper bound
// (unlike Ct, which is physically constrained to be >= 0), so without this
// cap a small fraction of draws — concentrated in higher-severity Shigella
// cells (~1-3%) — exceed it and convert to a negative, physically
// implausible Ct. Mirrors the existing lower-bound floor (quantity > 0).
// Matches MAX_QUANTITY in 2-functions.R exactly.
const double MAX_QUANTITY = 35.0 / 3.322;

inline double quantity_to_ct(double quantity) {
  return 35.0 - 3.322 * quantity;
}

// ---------------------------------------------------------------
// Helper: build log-odds from a named GLM coef vector.
// Uses named lookup so it is robust to any number of coefficients
// and any reference level ordering.
// extra_names      : names of continuous/binary predictor coefficients
//                     beyond agegrp/country (e.g. "score" for gemsdef
//                     model, or c("shigella_eiec", "shig.diar") for the
//                     unified culture model)
// extra_predictors  : corresponding predictor values, same length/order
//                     as extra_names
// ---------------------------------------------------------------
inline double logodds_from_coefs(const NumericVector& coefs,
                                  const CharacterVector& coef_names,
                                  const std::string& age_group,
                                  const std::string& country,
                                  const std::vector<std::string>& extra_names,
                                  const std::vector<double>& extra_predictors) {

  // Build a name->index map (coef vectors are tiny, ~12 elements)
  std::unordered_map<std::string, int> idx;
  for (int i = 0; i < coef_names.size(); i++)
    idx[as<std::string>(coef_names[i])] = i;

  double lp = coefs[idx.at("(Intercept)")];

  // Age group dummy (add if present and matches current age group)
  std::string age_key = "agegrp" + age_group;
  if (idx.count(age_key)) lp += coefs[idx.at(age_key)];

  // Country dummy (add if present and matches current country)
  std::string country_key = "country_id" + country;
  if (idx.count(country_key)) lp += coefs[idx.at(country_key)];

  // Continuous/binary predictors beyond agegrp/country
  for (size_t i = 0; i < extra_names.size(); i++) {
    if (idx.count(extra_names[i])) lp += coefs[idx.at(extra_names[i])] * extra_predictors[i];
  }

  return lp;
}

// Convenience overload for the common single-predictor case (gemsdef models)
inline double logodds_from_coefs(const NumericVector& coefs,
                                  const CharacterVector& coef_names,
                                  const std::string& age_group,
                                  const std::string& country,
                                  const std::string& extra_name,
                                  double extra_predictor) {
  return logodds_from_coefs(coefs, coef_names, age_group, country,
                            std::vector<std::string>{extra_name},
                            std::vector<double>{extra_predictor});
}

// ---------------------------------------------------------------
// Main simulation function
// [[Rcpp::export]]
DataFrame simulate_trial(int n_children, int n_months, CharacterVector countries,
                         DataFrame diarrhea_ir,
                         DataFrame shigella_sub,
                         DataFrame ETEC_sub,
                         DataFrame other_sub,
                         List shigella_sev_params,
                         List ETEC_sev_params,
                         List other_sev_params,
                         DataFrame shigella_quantity,
                         DataFrame ETEC_quantity,
                         DataFrame other_quantity) {

  // --- Unpack GLM coef vectors and their names (tiny, just copy once) ---
  List      shig_gems_params    = as<List>(shigella_sev_params["gems_params"]);
  List      shig_culture_params = as<List>(shigella_sev_params["culture_params"]);
  List      ETEC_gems_params    = as<List>(ETEC_sev_params["gems_params"]);
  List      ETEC_culture_params = as<List>(ETEC_sev_params["culture_params"]);
  List      other_gems_p        = as<List>(other_sev_params["gems_params"]);
  NumericVector shigella_gems_coefs       = shig_gems_params["coefs"];
  NumericVector shigella_culture_coefs    = shig_culture_params["coefs"];
  NumericVector ETEC_gems_coefs           = ETEC_gems_params["coefs"];
  NumericVector ETEC_culture_coefs        = ETEC_culture_params["coefs"];
  NumericVector other_gems_coefs          = other_gems_p["coefs"];
  // Named vectors: extract names for robust positional lookup
  CharacterVector shig_gems_names       = as<CharacterVector>(
    as<NumericVector>(shig_gems_params["coefs"]).attr("names"));
  CharacterVector shig_culture_names    = as<CharacterVector>(
    as<NumericVector>(shig_culture_params["coefs"]).attr("names"));
  CharacterVector ETEC_gems_names       = as<CharacterVector>(
    as<NumericVector>(ETEC_gems_params["coefs"]).attr("names"));
  CharacterVector ETEC_culture_names    = as<CharacterVector>(
    as<NumericVector>(ETEC_culture_params["coefs"]).attr("names"));
  CharacterVector other_gems_names = as<CharacterVector>(
    as<NumericVector>(other_gems_p["coefs"]).attr("names"));

  // --- Build all lookup maps once ---
  auto ir_map             = build_ir_map(diarrhea_ir);
  auto shig_score_map     = build_score_map(as<DataFrame>(shigella_sev_params["score_params"]));
  auto ETEC_score_map     = build_score_map(as<DataFrame>(ETEC_sev_params["score_params"]));
  auto other_sc_map       = build_score_map(as<DataFrame>(other_sev_params["score_params"]));
  auto shig_sub_map       = build_shig_subclin_map(shigella_sub);
  auto ETEC_sub_map       = build_ETEC_subclin_map(ETEC_sub);
  auto other_sub_map      = build_other_sub_map(other_sub);
  auto shig_qty_map       = build_shig_quantity_map(shigella_quantity);
  auto ETEC_qty_map       = build_ETEC_quantity_map(ETEC_quantity);
  auto other_qty_map      = build_other_quantity_map(other_quantity);

  // Calibrated severity-dependent VE shrinkage parameter (single scalar,
  // computed once per scenario in f.param() / calibrate_ve_max()).
  // Defaults to 0 (no shrinkage) if not present, e.g. for other_sev_params
  // which has no vax effect.
  double shig_ve_max = shigella_sev_params.containsElementNamed("ve_max") ?
                  as<double>(shigella_sev_params["ve_max"]) : 0.0;
  double ETEC_ve_max = ETEC_sev_params.containsElementNamed("ve_max") ?
                  as<double>(ETEC_sev_params["ve_max"]) : 0.0;

  int n_rows = n_children * n_months * countries.size();

  // --- Allocate output vectors ---
  CharacterVector country_id_out(n_rows);
  IntegerVector   child_out(n_rows);
  IntegerVector   month_out(n_rows);
  CharacterVector agegrp_out(n_rows);
  IntegerVector   vax_out(n_rows);
  IntegerVector   shigella_diarrhea(n_rows);
  IntegerVector   ETEC_diarrhea(n_rows);
  IntegerVector   other_diarrhea(n_rows);
  IntegerVector   shigella_score(n_rows);
  IntegerVector   shigella_gemsmsd(n_rows);
  IntegerVector   shigella_culture(n_rows);
  CharacterVector shigella_any_severity(n_rows);
  IntegerVector   ETEC_score(n_rows);
  IntegerVector   ETEC_gemsmsd(n_rows);
  IntegerVector   ETEC_culture(n_rows);
  CharacterVector ETEC_any_severity(n_rows);
  IntegerVector   other_score(n_rows);
  IntegerVector   other_gemsmsd(n_rows);
  CharacterVector other_any_severity(n_rows);
  IntegerVector   shigella_subclin(n_rows);
  DoubleVector    shig_quantity_out(n_rows);
  IntegerVector   ETEC_subclin(n_rows);
  DoubleVector    ETEC_quantity_out(n_rows);
  IntegerVector   other_inf(n_rows);
  DoubleVector    other_quantity_out(n_rows);

  int index = 0;

  for (int c = 0; c < countries.size(); c++) {
    std::string country = as<std::string>(countries[c]);

    for (int ch = 1; ch <= n_children; ch++) {
      int vax_status = static_cast<int>(R::rbinom(1, 0.5));

      for (int m = 1; m <= n_months; m++) {

        std::string ag = (m <= 6) ? "12-17 months" : "18-24 months";

        country_id_out[index] = country;
        child_out[index]      = ch;
        month_out[index]      = m;
        agegrp_out[index]     = ag;
        vax_out[index]        = vax_status;

        // --- Diarrhea events (map lookup, O(1)) ---
        const auto& ir = ir_map.at({country, ag});
        double ir_shig = (vax_status == 1) ? ir.vax_ir_shigella : ir.ir_shigella;
        double ir_ETEC = (vax_status == 1) ? ir.vax_ir_ETEC : ir.ir_ETEC;
        shigella_diarrhea[index] = diarrhea_event(ir_shig);
        ETEC_diarrhea[index] = diarrhea_event(ir_ETEC);
        other_diarrhea[index]    = diarrhea_event(ir.ir_other);

        // --- Shigella severity chain: score -> gemsdef (culture moved below,
        // after Ct is known, since the new unified culture model depends on
        // the drawn Ct value rather than gemsdef status) ---
        if (shigella_diarrhea[index] == 1) {
          const auto& sp = shig_score_map.at({country, ag});
          double drawn_score = std::round(
            std::min(std::max(R::rnorm(sp.mean_score, sp.sd_score), 1.0), 12.0));
          double s_score = (vax_status == 1) ?
            shrink_score(drawn_score, shig_ve_max) : drawn_score;
          shigella_score[index] = static_cast<int>(s_score);

          double lp_gems = logodds_from_coefs(shigella_gems_coefs, shig_gems_names, ag, country, "score", s_score);
          int s_gems = static_cast<int>(R::rbinom(1, 1.0 / (1.0 + std::exp(-lp_gems))));
          shigella_gemsmsd[index] = s_gems;

          shigella_any_severity[index] =
            (s_score >= 6 || s_gems == 1) ? "Severe" : "Mild";
        } else {
          shigella_score[index]        = NA_INTEGER;
          shigella_gemsmsd[index]      = NA_INTEGER;
          shigella_any_severity[index] = NA_STRING;
        }
        
        // --- ETEC severity chain: score -> gemsdef (culture moved below,
        // after Ct is known, since the new unified culture model depends on
        // the drawn Ct value rather than gemsdef status) ---
        if (ETEC_diarrhea[index] == 1) {
          const auto& sp = ETEC_score_map.at({country, ag});
          double drawn_score = std::round(
            std::min(std::max(R::rnorm(sp.mean_score, sp.sd_score), 1.0), 12.0));
          double s_score = (vax_status == 1) ?
          shrink_score(drawn_score, ETEC_ve_max) : drawn_score;
          ETEC_score[index] = static_cast<int>(s_score);
          
          double lp_gems = logodds_from_coefs(ETEC_gems_coefs, ETEC_gems_names, ag, country, "score", s_score);
          int s_gems = static_cast<int>(R::rbinom(1, 1.0 / (1.0 + std::exp(-lp_gems))));
          ETEC_gemsmsd[index] = s_gems;
          
          ETEC_any_severity[index] =
            (s_score >= 6 || s_gems == 1) ? "Severe" : "Mild";
        } else {
          ETEC_score[index]        = NA_INTEGER;
          ETEC_gemsmsd[index]      = NA_INTEGER;
          ETEC_any_severity[index] = NA_STRING;
        }

        // --- Other diarrhea severity ---
        if (other_diarrhea[index] == 1) {
          const auto& op = other_sc_map.at({country, ag});
          double o_score = std::round(
            std::min(std::max(R::rnorm(op.mean_score, op.sd_score), 1.0), 12.0));
          other_score[index] = static_cast<int>(o_score);

          double lp_og = logodds_from_coefs(other_gems_coefs, other_gems_names, ag, country, "score", o_score);
          int o_gems = static_cast<int>(R::rbinom(1, 1.0 / (1.0 + std::exp(-lp_og))));
          other_gemsmsd[index] = o_gems;

          other_any_severity[index] =
            (o_score >= 6 || o_gems == 1) ? "Severe" : "Mild";
        } else {
          other_score[index]        = NA_INTEGER;
          other_gemsmsd[index]      = NA_INTEGER;
          other_any_severity[index] = NA_STRING;
        }

        // --- Subclinical shigella (PCR detection, no diarrhea) ---
        int shig_prev_diarrhea = (m > 1 &&
          !IntegerVector::is_na(shigella_diarrhea[index - 1]) &&
          shigella_diarrhea[index - 1] == 1) ? 1 : 0;

        const auto& shig_sc = shig_sub_map.at({country, ag});
        double shig_sc_prob = (vax_status == 1) ?
          (shig_prev_diarrhea == 1 ? shig_sc.vax_pred_prob_1 : shig_sc.vax_pred_prob_0) :
          (shig_prev_diarrhea == 1 ? shig_sc.pred_prob_1     : shig_sc.pred_prob_0);
        shigella_subclin[index] = (NumericVector::is_na(shig_sc_prob) || shig_sc_prob < 0 ||
          shig_sc_prob > 1) ? 0 :
                                    static_cast<int>(R::rbinom(1, shig_sc_prob));

        // --- Shigella pathogen quantity: drawn for ANY PCR-detectable
        // Shigella this month (clinical diarrhea OR subclinical detection)
        // — must come before culture status since the unified culture
        // model uses Ct (converted from quantity) directly. Drawn from a
        // Gamma distribution (shape/scale parameterization, scale =
        // mean/shape) fit per country x agegrp x severity cell — see
        // 0-estimate-params.R. Quantity is clamped at a strict lower bound
        // (>0, via a small positive floor, since simulating an actual
        // infection implies quantity > 0) AND at MAX_QUANTITY (Ct >= 0),
        // since the Gamma model itself has no natural upper bound.
        shig_quantity_out[index] = NA_REAL;
        bool shig_pcr_positive = false;

        if (shigella_diarrhea[index] == 1 &&
            !CharacterVector::is_na(shigella_any_severity[index])) {
          std::string sev_str = as<std::string>(shigella_any_severity[index]);
          auto it = shig_qty_map.find({country, ag, sev_str});
          if (it != shig_qty_map.end()) {
            double mu_qty = (vax_status == 1) ? it->second.vax_mean_quantity : it->second.mean_quantity;
            double shape  = it->second.shape_quantity;
            double drawn_qty = R::rgamma(shape, mu_qty / shape);
            shig_quantity_out[index] = std::min(MAX_QUANTITY, std::max(1e-6, drawn_qty));
            shig_pcr_positive = true;
          }
        } else if (shigella_subclin[index] == 1) {
          auto it = shig_qty_map.find({country, ag, "Subclinical"});
          if (it != shig_qty_map.end()) {
            double mu_qty = (vax_status == 1) ? it->second.vax_mean_quantity : it->second.mean_quantity;
            double shape  = it->second.shape_quantity;
            double drawn_qty = R::rgamma(shape, mu_qty / shape);
            shig_quantity_out[index] = std::min(MAX_QUANTITY, std::max(1e-6, drawn_qty));
            shig_pcr_positive = true;
          }
        }

        // --- Shigella culture: unified model across ALL PCR-positive
        // specimens (clinical or subclinical), predicted from agegrp,
        // country, the drawn Ct value (converted from quantity), and
        // shig.diar status. Replaces the old two-path approach
        // (gemsdef-based for clinical, agegrp-only for subclinical).
        shigella_culture[index] = NA_INTEGER;
        if (shig_pcr_positive) {
          double shig_diar_flag = (shigella_diarrhea[index] == 1) ? 1.0 : 0.0;
          double ct_for_culture = quantity_to_ct(shig_quantity_out[index]);
          double lp_cult = logodds_from_coefs(
            shigella_culture_coefs, shig_culture_names, ag, country,
            std::vector<std::string>{"shigella_eiec", "shig.diar"},
            std::vector<double>{ct_for_culture, shig_diar_flag}
          );
          shigella_culture[index] = static_cast<int>(
            R::rbinom(1, 1.0 / (1.0 + std::exp(-lp_cult))));
        }
        
        
        // --- Subclinical ETEC (PCR detection, no diarrhea) ---
        int ETEC_prev_diarrhea = (m > 1 &&
                             !IntegerVector::is_na(ETEC_diarrhea[index - 1]) &&
                             ETEC_diarrhea[index - 1] == 1) ? 1 : 0;
        
        const auto& ETEC_sc = ETEC_sub_map.at({country, ag});
        double ETEC_sc_prob = (vax_status == 1) ?
        (ETEC_prev_diarrhea == 1 ? ETEC_sc.vax_pred_prob_1 : ETEC_sc.vax_pred_prob_0) :
          (ETEC_prev_diarrhea == 1 ? ETEC_sc.pred_prob_1     : ETEC_sc.pred_prob_0);
        ETEC_subclin[index] = (NumericVector::is_na(ETEC_sc_prob) || ETEC_sc_prob < 0 ||
          ETEC_sc_prob > 1) ? 0 :
          static_cast<int>(R::rbinom(1, ETEC_sc_prob));
        
        // --- ETEC pathogen quantity: drawn for ANY PCR-detectable
        // ETEC this month (clinical diarrhea OR subclinical detection)
        // — must come before culture status since the unified culture
        // model uses Ct (converted from quantity) directly. Drawn from a
        // Gamma distribution (shape/scale parameterization, scale =
        // mean/shape) fit per country x agegrp x severity cell — see
        // 0-estimate-params.R. Quantity is clamped at a strict lower bound
        // (>0, via a small positive floor, since simulating an actual
        // infection implies quantity > 0) AND at MAX_QUANTITY (Ct >= 0),
        // since the Gamma model itself has no natural upper bound.
        ETEC_quantity_out[index] = NA_REAL;
        bool ETEC_pcr_positive = false;
        
        if (ETEC_diarrhea[index] == 1 &&
            !CharacterVector::is_na(ETEC_any_severity[index])) {
            std::string sev_str = as<std::string>(ETEC_any_severity[index]);
          auto it = ETEC_qty_map.find({country, ag, sev_str});
          if (it != ETEC_qty_map.end()) {
            double mu_qty = (vax_status == 1) ? it->second.vax_mean_quantity : it->second.mean_quantity;
            double shape  = it->second.shape_quantity;
            double drawn_qty = R::rgamma(shape, mu_qty / shape);
            ETEC_quantity_out[index] = std::min(MAX_QUANTITY, std::max(1e-6, drawn_qty));
            ETEC_pcr_positive = true;
          }
        } else if (ETEC_subclin[index] == 1) {
          auto it = ETEC_qty_map.find({country, ag, "Subclinical"});
          if (it != ETEC_qty_map.end()) {
            double mu_qty = (vax_status == 1) ? it->second.vax_mean_quantity : it->second.mean_quantity;
            double shape  = it->second.shape_quantity;
            double drawn_qty = R::rgamma(shape, mu_qty / shape);
            ETEC_quantity_out[index] = std::min(MAX_QUANTITY, std::max(1e-6, drawn_qty));
            ETEC_pcr_positive = true;
          }
        }
        
        // --- ETEC culture: unified model across ALL PCR-positive
        // specimens (clinical or subclinical), predicted from agegrp,
        // country, the drawn Ct value (converted from quantity), and
        // ETEC.diar status. Replaces the old two-path approach
        // (gemsdef-based for clinical, agegrp-only for subclinical).
        ETEC_culture[index] = NA_INTEGER;
        if (ETEC_pcr_positive) {
          double ETEC_diar_flag = (ETEC_diarrhea[index] == 1) ? 1.0 : 0.0;
          double ct_for_culture = quantity_to_ct(ETEC_quantity_out[index]);
          double lp_cult = logodds_from_coefs(
            ETEC_culture_coefs, ETEC_culture_names, ag, country,
            std::vector<std::string>{"ST_ETEC", "ETEC.diar"},
            std::vector<double>{ct_for_culture, ETEC_diar_flag}
          );
          ETEC_culture[index] = static_cast<int>(
            R::rbinom(1, 1.0 / (1.0 + std::exp(-lp_cult))));
        }
        
        

        // --- Other pathogen detection and quantity ---
        // Detection probability now depends on other_diarrhea, mirroring
        // the Shigella subclinical model's use of prev.diarrhea (see
        // 0-estimate-params.R). Same severity structure and Gamma draw as
        // Shigella above: "Subclinical" when detected outside of an
        // other-diarrhea episode, otherwise the drawn other-diarrhea
        // severity (Mild/Severe). Clamped at a strict lower bound (>0,
        // since simulating an actual infection implies quantity > 0) AND
        // at MAX_QUANTITY (Ct >= 0), since the Gamma model itself has no
        // natural upper bound.
        const auto& os = other_sub_map.at({country, ag});
        double other_det_prob = (other_diarrhea[index] == 1) ?
          os.pred_prob_1 : os.pred_prob_0;
        other_inf[index] = NumericVector::is_na(other_det_prob) ? NA_INTEGER :
                           static_cast<int>(R::rbinom(1, other_det_prob));
        other_quantity_out[index] = NA_REAL;
        if (!IntegerVector::is_na(other_inf[index]) && other_inf[index] == 1 &&
            !IntegerVector::is_na(other_diarrhea[index])) {
          std::string sev_key;
          if (other_diarrhea[index] == 1 &&
              !CharacterVector::is_na(other_any_severity[index])) {
            sev_key = as<std::string>(other_any_severity[index]);
          } else if (other_diarrhea[index] == 0) {
            sev_key = "Subclinical";
          }
          if (!sev_key.empty()) {
            auto it = other_qty_map.find({country, ag, sev_key});
            if (it != other_qty_map.end()) {
              double shape = it->second.shape_quantity;
              double drawn_qty = R::rgamma(shape, it->second.mean_quantity / shape);
              other_quantity_out[index] = std::min(MAX_QUANTITY, std::max(1e-6, drawn_qty));
            }
          }
        }

        index++;
      }
    }
  }

  return DataFrame::create(
    Named("country_id")            = country_id_out,
    Named("child")                 = child_out,
    Named("month")                 = month_out,
    Named("agegrp")                = agegrp_out,
    Named("vax")                   = vax_out,
    Named("shigella_diarrhea")     = shigella_diarrhea,
    Named("ETEC_diarrhea")         = ETEC_diarrhea,
    Named("other_diarrhea")        = other_diarrhea,
    Named("shigella_score")        = shigella_score,
    Named("shigella_gemsmsd")      = shigella_gemsmsd,
    Named("shigella_culture")      = shigella_culture,
    Named("shigella_any_severity") = shigella_any_severity,
    Named("ETEC_score")            = ETEC_score,
    Named("ETEC_gemsmsd")          = ETEC_gemsmsd,
    Named("ETEC_culture")          = ETEC_culture,
    Named("ETEC_any_severity")     = ETEC_any_severity,
    Named("other_score")           = other_score,
    Named("other_gemsmsd")         = other_gemsmsd,
    Named("other_any_severity")    = other_any_severity,
    Named("shigella_sub")          = shigella_subclin,
    Named("shigella_quantity")     = shig_quantity_out,
    Named("ETEC_sub")              = ETEC_subclin,
    Named("ETEC_quantity")         = ETEC_quantity_out,
    Named("other_inf")             = other_inf,
    Named("other_quantity")        = other_quantity_out
  );
}