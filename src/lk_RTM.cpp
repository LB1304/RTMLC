// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

namespace {

arma::mat df_numeric_matrix(const DataFrame& data,
                            int start_col_zero_based,
                            int nrow,
                            int ncol,
                            bool take_log) {
  arma::mat out(nrow, ncol);
  for (int j = 0; j < ncol; ++j) {
    NumericVector col = as<NumericVector>(data[start_col_zero_based + j]);
    if (col.size() != nrow) {
      stop("A data column has a length different from nrow(data).");
    }
    for (int i = 0; i < nrow; ++i) {
      const double x = col[i];
      out(i, j) = take_log ? std::log(x) : x;
    }
  }
  return out;
}

void append_arma_vec(std::vector<double>& target, const arma::vec& x) {
  for (arma::uword i = 0; i < x.n_elem; ++i) target.push_back(x[i]);
}

inline double cov_dot(const arma::mat& XZ, const arma::vec& beta, int i) {
  double out = 0.0;
  for (arma::uword l = 0; l < beta.n_elem; ++l) {
    out += XZ(i, l) * beta[l];
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
List lk_RTM(NumericVector th,
            DataFrame data,
            int nc1,
            int nc2,
            int J) {
  const int p = nc1 + nc2;

  if (J <= 0 || nc1 < 0 || nc2 < 0) {
    stop("J must be positive; nc1 and nc2 must be non-negative.");
  }

  NumericVector cluster_r = as<NumericVector>(data["cluster"]);
  const int nt = cluster_r.size();
  if (data.size() < 2 + p + 2 * J) {
    stop("data does not have enough columns for the requested nc1, nc2, and J.");
  }

  arma::mat XZ(nt, p, arma::fill::zeros);
  if (p > 0) {
    XZ = df_numeric_matrix(data, 2, nt, p, false);
  }

  const int yy_start = 2 + p;
  const int lt_start = 2 + p + J;
  arma::mat YY = df_numeric_matrix(data, yy_start, nt, J, false);
  arma::mat LT = df_numeric_matrix(data, lt_start, nt, J, true);

  std::vector<int> obs_i;
  std::vector<int> obs_j;
  std::vector<double> YYobs;
  std::vector<double> LTobs;
  obs_i.reserve(nt * J);
  obs_j.reserve(nt * J);
  YYobs.reserve(nt * J);
  LTobs.reserve(nt * J);

  for (int i = 0; i < nt; ++i) {
    for (int j = 0; j < J; ++j) {
      const double y = YY(i, j);
      // R code uses which(!is.na(YYv)); ISNAN is true for NA and NaN.
      if (!ISNAN(y)) {
        obs_i.push_back(i);
        obs_j.push_back(j);
        YYobs.push_back(y);
        LTobs.push_back(LT(i, j));
      }
    }
  }
  const int ntobs = obs_i.size();

  const int np1 = p + J;
  const int np2 = p;
  const int np3 = J + (p > 0 ? (J - 1) : 0);
  const int npar_needed = np1 + np2 + np3 + 1;
  if (th.size() < npar_needed) {
    stop("th is shorter than expected for the supplied dimensions.");
  }

  arma::vec theta = as<arma::vec>(th);

  // ---- separate parameters ----
  arma::vec coef_resp = theta.subvec(0, np1 - 1); // p covariate effects + J item effects

  arma::vec coef1_time(p, arma::fill::zeros);
  if (p > 0) {
    coef1_time = theta.subvec(np1, np1 + np2 - 1);
  }

  arma::vec coef2_time = theta.subvec(np1 + np2, np1 + np2 + np3 - 1);
  arma::vec speed0_time = coef2_time.subvec(0, J - 1);
  arma::vec speed1_time(J > 1 ? J - 1 : 0, arma::fill::zeros);
  if (p > 0 && J > 1) {
    speed1_time = coef2_time.subvec(J, J + J - 2);
  }

  const double si2_time = theta[np1 + np2 + np3];
  const double sd_time = std::sqrt(si2_time);

  // ---- linear predictors ----
  arma::vec lin_resp(ntobs, arma::fill::zeros);
  arma::vec pred_resp(ntobs, arma::fill::zeros);
  arma::vec pred_time(ntobs, arma::fill::zeros);

  for (int r = 0; r < ntobs; ++r) {
    const int i = obs_i[r];
    const int j = obs_j[r];

    double eta = 0.0;
    for (int l = 0; l < p; ++l) eta += XZ(i, l) * coef_resp[l];
    eta -= coef_resp[p + j];
    lin_resp[r] = eta;
    pred_resp[r] = R::pnorm(eta, 0.0, 1.0, 1, 0);

    if (p == 0) {
      pred_time[r] = speed0_time[j];
    } else {
      const double base = cov_dot(XZ, coef1_time, i);
      if (j == 0) {
        pred_time[r] = speed0_time[0] - base;
      } else {
        pred_time[r] = speed0_time[j] - speed1_time[j - 1] * base;
      }
    }
  }

  // ---- compute log-likelihood ----
  double lk = 0.0;
  for (int r = 0; r < ntobs; ++r) {
    lk += R::dbinom(YYobs[r], 1.0, pred_resp[r], 1) +
      R::dnorm(LTobs[r], pred_time[r], sd_time, 1);
  }

  // ---- compute score ----
  std::vector<double> sc_all;
  sc_all.reserve(npar_needed);

  // parameters for the response model
  arma::vec sc_resp(np1, arma::fill::zeros);
  for (int r = 0; r < ntobs; ++r) {
    const int i = obs_i[r];
    const int j = obs_j[r];
    const double dens = R::dnorm(lin_resp[r], 0.0, 1.0, 0);
    const double denom = pred_resp[r] * (1.0 - pred_resp[r]);
    const double val = (YYobs[r] - pred_resp[r]) * dens / denom;

    for (int l = 0; l < p; ++l) sc_resp[l] += XZ(i, l) * val;
    sc_resp[p + j] -= val;
  }
  append_arma_vec(sc_all, sc_resp);

  // parameters for the covariate part of the time model
  if (p > 0) {
    arma::vec sc1_time(np2, arma::fill::zeros);
    for (int r = 0; r < ntobs; ++r) {
      const int i = obs_i[r];
      const int j = obs_j[r];
      const double factor = (j == 0) ? -1.0 : -speed1_time[j - 1];
      const double val = (LTobs[r] - pred_time[r]) / si2_time;
      for (int l = 0; l < p; ++l) sc1_time[l] += factor * XZ(i, l) * val;
    }
    append_arma_vec(sc_all, sc1_time);
  }

  // parameters for speed0.time and speed1.time
  arma::vec sc2_time(np3, arma::fill::zeros);
  for (int r = 0; r < ntobs; ++r) {
    const int i = obs_i[r];
    const int j = obs_j[r];
    const double val = (LTobs[r] - pred_time[r]) / si2_time;

    sc2_time[j] += val;
    if (p > 0 && j > 0) {
      const double base = cov_dot(XZ, coef1_time, i);
      sc2_time[J + j - 1] += (-base) * val;
    }
  }
  append_arma_vec(sc_all, sc2_time);

  // variance parameter si2.time
  double sc3_time = 0.0;
  for (int r = 0; r < ntobs; ++r) {
    const double resid = LTobs[r] - pred_time[r];
    sc3_time += -1.0 / (2.0 * si2_time) +
      resid * resid / (2.0 * si2_time * si2_time);
  }
  sc_all.push_back(sc3_time);

  return List::create(
    Named("lk") = lk,
    Named("sc") = NumericVector(sc_all.begin(), sc_all.end())
  );
}
