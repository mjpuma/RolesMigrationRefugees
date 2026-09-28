// =============================================================================
// rege_valued.cpp — White's REGGE for valued (weighted) directed networks
//
// Tie values are the flows (optionally log-transformed in R). For countries i
// and j, every partner k of i is matched to the partner m of j that maximises
//
//      E[k,m] * ( min(x_ik, x_jm) + min(x_ki, x_mj) )
//
// i.e. an equivalent partner (E[k,m]) linked by ties of similar size in both
// directions. The denominator adds max(x_ik, x_jm) + max(x_ki, x_mj) for the
// chosen match (or i's own tie values when nothing matches), so a pair scores
// 1 only if every tie is matched by an equally large tie to an equivalent
// partner. The same is done from j's side, and
//
//      E_new[i,j] = (num_i + num_j) / (den_i + den_j)
//
// no_self_match = TRUE: a partner k of i may not be matched to the SAME node
// as a partner of j (optional; default FALSE, i.e. standard REGGE).
//
// Reference: White (1985) REGGE; Borgatti & Everett (1993) Two algorithms for
// computing regular equivalence. Social Networks 15: 361-376.
// =============================================================================
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// one-sided sum over partners k of i, matched against partners m of j
static inline void side(int i, int j,
                        const std::vector<std::vector<int>>& nb,
                        const std::vector<std::vector<double>>& xo,
                        const std::vector<std::vector<double>>& xi,
                        const NumericMatrix& E, double& num, double& den,
                        bool no_self) {
  const std::vector<int>& Ni = nb[i];
  const std::vector<int>& Nj = nb[j];
  for (size_t a = 0; a < Ni.size(); ++a) {
    int k = Ni[a];
    double o1 = xo[i][a], i1 = xi[i][a];
    double best = 0.0, bestMax = o1 + i1;       // unmatched: own tie values
    for (size_t b = 0; b < Nj.size(); ++b) {
      int m = Nj[b];
      if (no_self && m == k) continue;   // a shared partner may not match itself
      double e = E(k, m);
      if (e <= 0.0) continue;
      double o2 = xo[j][b], i2 = xi[j][b];
      double mn = std::min(o1, o2) + std::min(i1, i2);
      double val = e * mn;
      double mx = std::max(o1, o2) + std::max(i1, i2);
      // choose the partner with the best weighted match; among equally good
      // matches take the one whose tie size is closest (smallest max), so an
      // identical tie always wins. (Keeping the first tie found was a bug: at
      // E = 1 a much larger tie scores the same numerator but inflates the
      // denominator, so even identical networks scored < 1.)
      if (val > best + 1e-12 || (std::abs(val - best) <= 1e-12 && mx < bestMax)) {
        best = val;
        bestMax = mx;
      }
    }
    num += best;
    den += bestMax;
  }
}

// [[Rcpp::export]]
List rege_valued_cpp(NumericMatrix X, int iterations, bool no_self_match) {
  int n = X.nrow();
  std::vector<std::vector<int>> nb(n);
  std::vector<std::vector<double>> xo(n), xi(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      if (k == i) continue;
      double out = X(i, k), in = X(k, i);
      if (out > 0 || in > 0) {
        nb[i].push_back(k); xo[i].push_back(out); xi[i].push_back(in);
      }
    }
  }
  NumericMatrix E(n, n);
  std::fill(E.begin(), E.end(), 1.0);
  NumericVector change(iterations);
  for (int it = 0; it < iterations; ++it) {
    NumericMatrix En(n, n);
    double maxch = 0.0;
    for (int i = 0; i < n; ++i) {
      En(i, i) = 1.0;
      for (int j = i + 1; j < n; ++j) {
        double num = 0.0, den = 0.0;
        side(i, j, nb, xo, xi, E, num, den, no_self_match);
        side(j, i, nb, xo, xi, E, num, den, no_self_match);
        double v = den > 0 ? num / den : 1.0;
        En(i, j) = v; En(j, i) = v;
        maxch = std::max(maxch, std::abs(v - E(i, j)));
      }
      if (i % 100 == 0) Rcpp::checkUserInterrupt();
    }
    E = En;
    change[it] = maxch;
    Rcpp::Rcout << "    valued REGE iteration " << it + 1
                << ": max change " << maxch << std::endl;
  }
  return List::create(_["similarity"] = E, _["change"] = change);
}
