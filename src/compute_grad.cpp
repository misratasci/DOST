#include <Rcpp.h>
using namespace Rcpp;

//' Compute Gradient (Internal C++ function)
//'
//' @param X Numeric Matrix
//' @param D_x Numeric Matrix
//' @param a Numeric Matrix
//' @param b Numeric Matrix
//' @param eps Double
//' @return Gradient matrix
NumericMatrix compute_grad(const NumericMatrix& X,
                          const NumericMatrix& D_x,
                          const NumericMatrix& a,
                          const NumericMatrix& b,
                          const double eps = 1e-20)
{
 const int N = X.nrow();
 const int K = X.ncol();

 NumericMatrix first_term(N, N);
 for (int i = 0; i < N; ++i) {
   for (int j = 0; j < N; ++j) {
     first_term(i, j) = D_x(i, j) * a(i, j) + b(i, j);
   }
 }

 NumericMatrix denom(N, N);
 for (int i = 0; i < N; ++i) {
   denom(i, i) = eps;
   for (int j = i + 1; j < N; ++j) {
     double ss = 0.0;
     for (int k = 0; k < K; ++k) {
       double diff = X(i, k) - X(j, k);
       ss += diff * diff;
     }
     double d = std::sqrt(ss);
     if (d < eps) d = eps;
     denom(i, j) = d;
     denom(j, i) = d;
   }
 }

 NumericMatrix grad(N, K);
 for (int i = 0; i < N; ++i) {
   for (int j = 0; j < N; ++j) {
     if (i == j) continue;
     double coeff = first_term(i, j) / denom(i, j);
     for (int k = 0; k < K; ++k) {
       grad(i, k) += (X(i, k) - X(j, k)) * coeff;
     }
   }
 }

 for (int i = 0; i < N; ++i)
   for (int k = 0; k < K; ++k)
     grad(i, k) *= 2.0;

 return grad;
}
