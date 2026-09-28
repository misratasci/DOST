#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static void grad_block(const double* X, const double* dxj, const double* dej,
                       const int* V_p, const int* V_i, const int j,
                       const int i0, const int i1, const int N, const int K,
                       const double N2, const double b0, const double b1,
                       const double eps, double* G, double* coeff)
{
 const int m = i1 - i0;
 for (int t = 0; t < m; ++t) coeff[t] = b0;
 for (int p = V_p[j]; p < V_p[j + 1]; ++p) {
   const int i = V_i[p];
   if (i >= i0 && i < i1) coeff[i - i0] = b1;
 }
 for (int t = 0; t < m; ++t) {
   const double b = coeff[t];
   const double dx = dxj[i0 + t];
   const double a = 2.0 / (N2 * dej[i0 + t]);
   double first_term = dx * a + b;
   const double d = (dx < eps) ? eps : dx;
   coeff[t] = first_term / d;
 }
 for (int k = 0; k < K; ++k) {
   const double* xk = X + (size_t) k * N + i0;
   const double xjk = X[(size_t) k * N + j];
   double* gk = G + (size_t) k * N + i0;
   for (int t = 0; t < m; ++t) {
     gk[t] += (xk[t] - xjk) * coeff[t];
   }
 }
}

//' Compute Gradient (Internal C++ function)
//'
//'
//' @param X Numeric Matrix (N x K embedding)
//' @param D_x Numeric Matrix (pairwise distances of X)
//' @param D_expr Numeric Matrix (expression distances)
//' @param V_p Integer Vector (column pointers of V)
//' @param V_i Integer Vector (row indices of V)
//' @param lambda Double
//' @param sum_V Double
//' @param eps Double
//' @return Gradient matrix
//' @noRd
// [[Rcpp::export]]
NumericMatrix compute_grad(const NumericMatrix& X,
                           const NumericMatrix& D_x,
                           const NumericMatrix& D_expr,
                           const IntegerVector& V_p,
                           const IntegerVector& V_i,
                           const double lambda,
                           const double sum_V,
                           const double eps = 1e-20)
{
 const int N = X.nrow();
 const int K = X.ncol();
 const double N2 = (double) N * (double) N;
 const double b0 = (lambda * 0.0) / sum_V - 2.0 / N2;
 const double b1 = (lambda * 1.0) / sum_V - 2.0 / N2;
 const int BLOCK = 512;

 NumericMatrix grad(N, K);
 double* G = grad.begin();
 const double* x = X.begin();
 const double* dx = D_x.begin();
 const double* de = D_expr.begin();
 const int* vp = V_p.begin();
 const int* vi = V_i.begin();

 std::vector<double> coeff(BLOCK);
 for (int j = 0; j < N; ++j) {
   const double* dxj = dx + (size_t) j * N;
   const double* dej = de + (size_t) j * N;
   for (int s = 0; s < N; s += BLOCK) {
     const int e = std::min(s + BLOCK, N);
     if (j >= s && j < e) {
       if (j > s) grad_block(x, dxj, dej, vp, vi, j, s, j, N, K, N2, b0, b1,
                             eps, G, coeff.data());
       if (j + 1 < e) grad_block(x, dxj, dej, vp, vi, j, j + 1, e, N, K, N2,
                                 b0, b1, eps, G, coeff.data());
     } else {
       grad_block(x, dxj, dej, vp, vi, j, s, e, N, K, N2, b0, b1,
                  eps, G, coeff.data());
     }
   }
 }

 for (R_xlen_t q = 0; q < (R_xlen_t) N * K; ++q) G[q] *= 2.0;

 return grad;
}

// Loss sums (Internal C++ function)
//
// [[Rcpp::export]]
NumericVector loss_sums(const NumericMatrix& D_lab,
                        const NumericMatrix& D_expr,
                        const IntegerVector& V_p,
                        const IntegerVector& V_i,
                        const double eps = 1e-20)
{
 const int N = D_lab.nrow();
 const double* dl = D_lab.begin();
 const double* de = D_expr.begin();
 long double s1 = 0.0, s2 = 0.0;
 for (int j = 0; j < N; ++j) {
   const size_t off = (size_t) j * N;
   for (int i = 0; i < N; ++i) {
     double diff = dl[off + i] - de[off + i];
     double sq = diff * diff;
     double denom = (i == j) ? de[off + i] + eps : de[off + i];
     s1 += sq / denom;
   }
   for (int p = V_p[j]; p < V_p[j + 1]; ++p) {
     s2 += dl[off + V_i[p]];
   }
 }
 return NumericVector::create((double) s1, (double) s2);
}

// Double centering for classical MDS (Internal C++ function)
// [[Rcpp::export]]
NumericMatrix double_center(const NumericMatrix& x,
                            const NumericVector& rmean,
                            const NumericVector& cmean,
                            const double gmean)
{
 const int n = x.nrow();
 NumericMatrix out(n, n);
 for (int j = 0; j < n; ++j) {
   const double c = cmean[j];
   for (int i = 0; i < n; ++i) {
     double y = (x(i, j) - rmean[i]) - c;
     y = y + gmean;
     out(i, j) = -y / 2.0;
   }
 }
 return out;
}

static inline double dist_at(const double* d, const size_t n, size_t i, size_t j)
{
 // element (i, j), i < j, of a 'dist' object (lower triangle by columns)
 return d[n * i - i * (i + 1) / 2 + j - i - 1];
}

// Nearest-neighbour distance of each point from a 'dist' vector
// [[Rcpp::export]]
NumericVector dist_row_min(const NumericVector& d, const int n)
{
 NumericVector r(n, R_PosInf);
 const double* dp = d.begin();
 for (int i = 0; i < n; ++i) {
   for (int j = i + 1; j < n; ++j) {
     double v = dist_at(dp, n, i, j);
     if (v < r[i]) r[i] = v;
     if (v < r[j]) r[j] = v;
   }
 }
 return r;
}

// Sparse 0/1 adjacency (d <= threshold, diagonal included) from a 'dist' vector
// [[Rcpp::export]]
List dist_adjacency(const NumericVector& d, const int n, const double threshold)
{
 const double* dp = d.begin();
 std::vector<int> p(n + 1, 0), idx;
 const bool self = (0.0 <= threshold);
 for (int j = 0; j < n; ++j) {
   for (int i = 0; i < n; ++i) {
     bool on;
     if (i == j) on = self;
     else if (i < j) on = dist_at(dp, n, i, j) <= threshold;
     else on = dist_at(dp, n, j, i) <= threshold;
     if (on) idx.push_back(i);
   }
   p[j + 1] = (int) idx.size();
 }
 return List::create(Named("p") = wrap(p),
                     Named("i") = wrap(idx),
                     Named("n") = n,
                     Named("sum") = (double) idx.size());
}

// Majority vote among the k nearest neighbours (Internal C++ function)
//
// [[Rcpp::export]]
IntegerVector knn_vote(const NumericMatrix& D, const IntegerVector& codes,
                       const int nlevels, const int k)
{
 const int n = D.nrow();
 IntegerVector out(n);
 std::vector<int> cand;
 cand.reserve(n);
 std::vector<double> row(n);
 std::vector<int> counts(nlevels + 1);
 for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) row[j] = D(i, j);
   cand.clear();
   for (int j = 0; j < n; ++j) if (j != i) cand.push_back(j);
   auto less = [&row](int a, int b) {
     return row[a] < row[b] || (row[a] == row[b] && a < b);
   };
   std::nth_element(cand.begin(), cand.begin() + (k - 1), cand.end(), less);
   std::fill(counts.begin(), counts.end(), 0);
   for (int m = 0; m < k; ++m) counts[codes[cand[m]]]++;
   int best = 1;
   for (int l = 2; l <= nlevels; ++l) if (counts[l] > counts[best]) best = l;
   out[i] = best;
 }
 return out;
}
