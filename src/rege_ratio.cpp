// =============================================================================
// rege_ratio.cpp — REGE on valued ties as written in Luczkovich, Borgatti,
// Johnson & Everett (2003, J. theor. Biol. 220: 303-321), Section 3.4
//
// Start with r_ij = 1 for all pairs. In each iteration, for every pair (i, j):
//
//   (1A) for each node k that i sends to, find the node m that j sends to which
//        is most equivalent to k AND receives a similar amount, i.e. maximise
//            z_k = r_km * min(x_ik, x_jm) / max(x_ik, x_jm)
//   (1B) for each node k that sends to i, find the node m that sends to j
//        maximising
//            y_k = r_km * min(x_ki, x_mj) / max(x_ki, x_mj)
//   (1C) r_ij = sum_k z_k + sum_k y_k
//
// Each tie counts once, whatever its size: a tie of 3 people matched by a tie
// of 3 people scores 1, as does 300,000 matched by 300,000. The paper gives the
// coefficients on a 0-1 scale but not the normalisation; here the sums are
// taken from both sides (i's ties matched in j, and j's ties matched in i) and
// divided by the total number of ties, |out(i)| + |in(i)| + |out(j)| + |in(j)|.
//
// (2) Repeat until no coefficient changes by more than `tol`, or `max_iter`.
//
// For comparison, src/rege_valued.cpp is White's REGGE as implemented in UCINET,
// in which ties are weighted by their size and in- and out-ties to the same
// partner are matched together.
//
// Compiled with OpenMP when the compiler supports it (pairs are independent
// within an iteration); otherwise it runs on one core.
// =============================================================================
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

struct Ties { std::vector<int> node; std::vector<double> val; };

// sum over i's ties (in one direction) of the best match among j's ties
static inline double best_sum(const Ties& Ti, const Ties& Tj, const std::vector<double>& E, int n) {
  double s = 0.0;
  for (size_t a = 0; a < Ti.node.size(); ++a) {
    const int k = Ti.node[a]; const double xa = Ti.val[a];
    const double* Ek = &E[(size_t)k * n];
    double best = 0.0;
    for (size_t b = 0; b < Tj.node.size(); ++b) {
      const double e = Ek[Tj.node[b]];
      if (e <= best) continue;                     // cannot beat the current best
      const double xb = Tj.val[b];
      const double v = e * (xa < xb ? xa / xb : xb / xa);
      if (v > best) best = v;
    }
    s += best;
  }
  return s;
}

//   keep  iterations at which to return a copy of the similarity matrix (to study
//         how the result depends on the number of iterations)
// [[Rcpp::export]]
List rege_ratio_cpp(NumericMatrix X, int max_iter = 50, double tol = 1e-4, int threads = 1,
                    IntegerVector keep = IntegerVector::create()) {
  const int n = X.nrow();
  std::vector<Ties> out(n), inn(n);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k) {
      if (k == i) continue;
      if (X(i, k) > 0) { out[i].node.push_back(k); out[i].val.push_back(X(i, k)); }
      if (X(k, i) > 0) { inn[i].node.push_back(k); inn[i].val.push_back(X(k, i)); }
    }
  std::vector<double> E((size_t)n * n, 1.0), En((size_t)n * n, 1.0);
  std::vector<double> change;
  List snaps; CharacterVector snap_names;
  auto to_matrix = [&](const std::vector<double>& M) {
    NumericMatrix S(n, n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) S(i, j) = M[(size_t)i * n + j];
    return S;
  };
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif
  for (int it = 0; it < max_iter; ++it) {
    double maxch = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(max:maxch)
#endif
    for (int i = 0; i < n; ++i) {
      En[(size_t)i * n + i] = 1.0;
      for (int j = i + 1; j < n; ++j) {
        const double den = out[i].node.size() + inn[i].node.size() + out[j].node.size() + inn[j].node.size();
        double v = 1.0;
        if (den > 0) {
          const double num = best_sum(out[i], out[j], E, n) + best_sum(inn[i], inn[j], E, n) +
                             best_sum(out[j], out[i], E, n) + best_sum(inn[j], inn[i], E, n);
          v = num / den;
        }
        En[(size_t)i * n + j] = v; En[(size_t)j * n + i] = v;
        const double d = std::abs(v - E[(size_t)i * n + j]);
        if (d > maxch) maxch = d;
      }
    }
    std::swap(E, En);
    change.push_back(maxch);
    for (int q = 0; q < keep.size(); ++q) if (keep[q] == it + 1) { snaps.push_back(to_matrix(E)); snap_names.push_back(std::to_string(it + 1)); }
    Rcpp::Rcout << "    ratio REGE iteration " << it + 1 << ": max change " << maxch << std::endl;
    Rcpp::checkUserInterrupt();
    if (maxch < tol) break;
  }
  snaps.names() = snap_names;
  return List::create(_["similarity"] = to_matrix(E), _["change"] = wrap(change), _["snapshots"] = snaps);
}
