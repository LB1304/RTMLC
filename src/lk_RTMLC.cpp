// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

namespace {

inline double clamp_prob(double x) {
  if (x < 1e-10) return 1e-10;
  if (x > 1.0 - 1e-10) return 1.0 - 1e-10;
  return x;
}

inline double cov_dot(const arma::mat& XZ, const arma::vec& beta, int i) {
  double out = 0.0;
  for (arma::uword l = 0; l < beta.n_elem; ++l) {
    out += XZ(i, l) * beta[l];
  }
  return out;
}

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

} // namespace

// [[Rcpp::export]]
List lk_RTMLC(NumericVector th,
              DataFrame data,
              int nc1,
              int nc2,
              int J,
              int k1,
              int k2) {
  const int p = nc1 + nc2;
  const int sup_count = k1 + k2 - 2;

  if (J <= 0 || k1 <= 0 || k2 <= 0 || nc1 < 0 || nc2 < 0) {
    stop("J, k1, and k2 must be positive; nc1 and nc2 must be non-negative.");
  }

  NumericVector cluster_r = as<NumericVector>(data["cluster"]);
  NumericVector individual_r = as<NumericVector>(data["individual"]);
  const int nt = cluster_r.size();
  if (individual_r.size() != nt) {
    stop("data$cluster and data$individual must have the same length.");
  }
  if (data.size() < 2 + p + 2 * J) {
    stop("data does not have enough columns for the requested nc1, nc2, and J.");
  }

  std::vector<int> cluster(nt);
  int H = 0;
  for (int i = 0; i < nt; ++i) {
    if (ISNAN(cluster_r[i])) stop("data$cluster contains NA/NaN values.");
    cluster[i] = static_cast<int>(cluster_r[i]);
    if (cluster[i] < 1) stop("data$cluster must use positive, 1-based cluster labels.");
    if (cluster[i] > H) H = cluster[i];
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

  const int np1 = p + J + sup_count;
  const int np2 = p + sup_count;
  const int np3 = J + (np2 > 0 ? (J - 1) : 0);
  const int npar_needed = np1 + np2 + np3 + 1 + (k1 - 1) + (k2 - 1);
  if (th.size() < npar_needed) {
    stop("th is shorter than expected for the supplied dimensions.");
  }

  arma::vec theta = as<arma::vec>(th);

  // ---- separate parameters ----
  arma::vec coef_resp = theta.subvec(0, np1 - 1);
  arma::vec sup1_resp(k1 > 1 ? k1 - 1 : 0, arma::fill::zeros);
  arma::vec sup2_resp(k2 > 1 ? k2 - 1 : 0, arma::fill::zeros);
  if (k1 > 1) sup1_resp = coef_resp.subvec(0, k1 - 2);
  if (k2 > 1) sup2_resp = coef_resp.subvec(k1 - 1, k1 + k2 - 3);
  arma::vec coef_resp_main = coef_resp.subvec(sup_count, np1 - 1); // p + J

  arma::vec coef1_time(np2, arma::fill::zeros);
  if (np2 > 0) coef1_time = theta.subvec(np1, np1 + np2 - 1);
  arma::vec sup1_time(k1 > 1 ? k1 - 1 : 0, arma::fill::zeros);
  arma::vec sup2_time(k2 > 1 ? k2 - 1 : 0, arma::fill::zeros);
  if (k1 > 1) sup1_time = coef1_time.subvec(0, k1 - 2);
  if (k2 > 1) sup2_time = coef1_time.subvec(k1 - 1, k1 + k2 - 3);
  arma::vec coef1_time_main(p, arma::fill::zeros);
  if (p > 0) coef1_time_main = coef1_time.subvec(sup_count, np2 - 1);

  arma::vec coef2_time = theta.subvec(np1 + np2, np1 + np2 + np3 - 1);
  arma::vec speed0_time = coef2_time.subvec(0, J - 1);
  arma::vec speed1_time(J > 1 ? J - 1 : 0, arma::fill::zeros);
  if (np3 > J && J > 1) speed1_time = coef2_time.subvec(J, J + J - 2);

  const double si2_time = theta[np1 + np2 + np3];

  arma::vec piv1(k1, arma::fill::ones);
  if (k1 > 1) {
    piv1.zeros();
    piv1[0] = 1.0;
    const int start = np1 + np2 + np3 + 1;
    for (int u = 1; u < k1; ++u) piv1[u] = std::exp(theta[start + u - 1]);
    piv1 /= arma::accu(piv1);
  }

  arma::vec piv2(k2, arma::fill::ones);
  if (k2 > 1) {
    piv2.zeros();
    piv2[0] = 1.0;
    const int start = np1 + np2 + np3 + k1;
    for (int v = 1; v < k2; ++v) piv2[v] = std::exp(theta[start + v - 1]);
    piv2 /= arma::accu(piv2);
  }

  // Linear predictor for the response part: DD.resp %*% coef.resp_main.
  arma::vec lin_resp(ntobs, arma::fill::zeros);
  for (int r = 0; r < ntobs; ++r) {
    const int i = obs_i[r];
    const int j = obs_j[r];
    double eta = 0.0;
    for (int l = 0; l < p; ++l) eta += XZ(i, l) * coef_resp_main[l];
    eta -= coef_resp_main[p + j];
    lin_resp[r] = eta;
  }

  auto class_shift_resp = [&](int u, int v) {
    double out = 0.0;
    if (u > 0) out += sup1_resp[u - 1];
    if (v > 0) out += sup2_resp[v - 1];
    return out;
  };

  auto class_shift_time = [&](int u, int v) {
    double out = 0.0;
    if (u > 0) out += sup1_time[u - 1];
    if (v > 0) out += sup2_time[v - 1];
    return out;
  };

  auto time_base = [&](int i, int u, int v) {
    double out = class_shift_time(u, v);
    if (p > 0) out += cov_dot(XZ, coef1_time_main, i);
    return out;
  };

  auto time_mean = [&](int i, int j, int u, int v) {
    if (np2 == 0) return speed0_time[j];
    const double base = time_base(i, u, v);
    if (j == 0) return speed0_time[0] - base;
    return speed0_time[j] - speed1_time[j - 1] * base;
  };

  // ---- compute log-likelihood ----
  arma::cube lQcc(nt, k1, k2, arma::fill::zeros);
  const double sd_time = std::sqrt(si2_time);

  for (int u = 0; u < k1; ++u) {
    for (int v = 0; v < k2; ++v) {
      const double shift_r = class_shift_resp(u, v);
      for (int r = 0; r < ntobs; ++r) {
        const int i = obs_i[r];
        const int j = obs_j[r];
        const double pred_resp = clamp_prob(R::pnorm(shift_r + lin_resp[r], 0.0, 1.0, 1, 0));
        const double pred_time = time_mean(i, j, u, v);
        lQcc(i, u, v) += R::dbinom(YYobs[r], 1.0, pred_resp, 1) +
          R::dnorm(LTobs[r], pred_time, sd_time, 1);
      }
    }
  }

  arma::mat lQcm(nt, k1, arma::fill::zeros);
  for (int i = 0; i < nt; ++i) {
    for (int u = 0; u < k1; ++u) {
      double mx = lQcc(i, u, 0);
      for (int v = 1; v < k2; ++v) mx = std::max(mx, lQcc(i, u, v));
      double s = 0.0;
      for (int v = 0; v < k2; ++v) s += std::exp(lQcc(i, u, v) - mx) * piv2[v];
      lQcm(i, u) = std::log(s) + mx;
    }
  }

  arma::mat lPc(H, k1, arma::fill::zeros);
  for (int i = 0; i < nt; ++i) {
    const int h = cluster[i] - 1;
    for (int u = 0; u < k1; ++u) lPc(h, u) += lQcm(i, u);
  }

  arma::vec lPm(H, arma::fill::zeros);
  double lk = 0.0;
  for (int h = 0; h < H; ++h) {
    double mx = lPc(h, 0);
    for (int u = 1; u < k1; ++u) mx = std::max(mx, lPc(h, u));
    double s = 0.0;
    for (int u = 0; u < k1; ++u) s += std::exp(lPc(h, u) - mx) * piv1[u];
    lPm[h] = std::log(s) + mx;
    lk += lPm[h];
  }

  // ---- compute score: E-step ----
  arma::mat Pp(H, k1, arma::fill::zeros);
  for (int h = 0; h < H; ++h) {
    for (int u = 0; u < k1; ++u) {
      Pp(h, u) = std::exp(lPc(h, u) - lPm[h]) * piv1[u];
    }
  }

  arma::cube Qpp(nt, k1, k2, arma::fill::zeros);
  arma::mat Qmp(nt, k2, arma::fill::zeros);
  for (int i = 0; i < nt; ++i) {
    const int h = cluster[i] - 1;
    for (int u = 0; u < k1; ++u) {
      for (int v = 0; v < k2; ++v) {
        const double qcp = std::exp(lQcc(i, u, v) - lQcm(i, u)) * piv2[v];
        Qpp(i, u, v) = Pp(h, u) * qcp;
        Qmp(i, v) += Qpp(i, u, v);
      }
    }
  }

  std::vector<double> sc_all;
  sc_all.reserve(npar_needed);

  // ---- parameters for response ----
  arma::vec sc_resp(np1, arma::fill::zeros);
  for (int u = 0; u < k1; ++u) {
    for (int v = 0; v < k2; ++v) {
      const double shift_r = class_shift_resp(u, v);
      for (int r = 0; r < ntobs; ++r) {
        const int i = obs_i[r];
        const int j = obs_j[r];
        const double eta = shift_r + lin_resp[r];
        const double pred = clamp_prob(R::pnorm(eta, 0.0, 1.0, 1, 0));
        const double dens = R::dnorm(eta, 0.0, 1.0, 0);
        const double denom = pred * (1.0 - pred);
        const double val = Qpp(i, u, v) * (YYobs[r] - pred) * dens / denom;

        if (u > 0) sc_resp[u - 1] += val;
        if (v > 0) sc_resp[(k1 - 1) + (v - 1)] += val;
        for (int l = 0; l < p; ++l) sc_resp[sup_count + l] += XZ(i, l) * val;
        sc_resp[sup_count + p + j] -= val;
      }
    }
  }
  append_arma_vec(sc_all, sc_resp);

  // ---- parameters of the time model: coef1.time ----
  arma::vec sc1_time(np2, arma::fill::zeros);
  if (np2 > 0) {
    for (int u = 0; u < k1; ++u) {
      for (int v = 0; v < k2; ++v) {
        for (int r = 0; r < ntobs; ++r) {
          const int i = obs_i[r];
          const int j = obs_j[r];
          const double factor = (j == 0) ? -1.0 : -speed1_time[j - 1];
          const double base = time_base(i, u, v);
          const double resid = LTobs[r] - speed0_time[j] - factor * base;
          const double val = Qpp(i, u, v) * resid / si2_time;

          if (u > 0) sc1_time[u - 1] += factor * val;
          if (v > 0) sc1_time[(k1 - 1) + (v - 1)] += factor * val;
          for (int l = 0; l < p; ++l) sc1_time[sup_count + l] += factor * XZ(i, l) * val;
        }
      }
    }
  }
  append_arma_vec(sc_all, sc1_time);

  // ---- parameters of the time model: coef2.time ----
  arma::vec sc2_time(np3, arma::fill::zeros);
  for (int u = 0; u < k1; ++u) {
    for (int v = 0; v < k2; ++v) {
      for (int r = 0; r < ntobs; ++r) {
        const int i = obs_i[r];
        const int j = obs_j[r];
        const double pred = time_mean(i, j, u, v);
        const double resid = LTobs[r] - pred;
        const double val = Qpp(i, u, v) * resid / si2_time;

        sc2_time[j] += val;
        if (np3 > J && j > 0) {
          const double base = time_base(i, u, v);
          sc2_time[J + j - 1] += (-base) * val;
        }
      }
    }
  }
  append_arma_vec(sc_all, sc2_time);

  // ---- variance parameter si2.time ----
  double sc3_time = 0.0;
  for (int u = 0; u < k1; ++u) {
    for (int v = 0; v < k2; ++v) {
      for (int r = 0; r < ntobs; ++r) {
        const int i = obs_i[r];
        const int j = obs_j[r];
        const double resid = LTobs[r] - time_mean(i, j, u, v);
        const double w = Qpp(i, u, v);
        sc3_time += -w / (2.0 * si2_time) + w * resid * resid / (2.0 * si2_time * si2_time);
      }
    }
  }
  sc_all.push_back(sc3_time);

  // ---- cluster weights ----
  if (k1 > 1) {
    for (int u = 1; u < k1; ++u) {
      sc_all.push_back(arma::accu(Pp.col(u)) - static_cast<double>(H) * piv1[u]);
    }
  }
  if (k2 > 1) {
    for (int v = 1; v < k2; ++v) {
      sc_all.push_back(arma::accu(Qmp.col(v)) - static_cast<double>(nt) * piv2[v]);
    }
  }

  return List::create(
    Named("lk") = lk,
    Named("sc") = NumericVector(sc_all.begin(), sc_all.end())
  );
}
