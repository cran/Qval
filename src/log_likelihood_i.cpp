#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double log_likelihood_i(const NumericVector& Yi, const NumericVector& P_Xj_alpha, 
                                const NumericMatrix& P_alpha_Xi) {
  int L = P_Xj_alpha.size();
  int N = Yi.size();
  
  // Step 1: Create the Y.temp matrix (N x L)
  NumericMatrix Y_temp(N, L);
  for (int n = 0; n < N; ++n) {
    for (int l = 0; l < L; ++l) {
      Y_temp(n, l) = Yi[n];  // Fill column i of Y into Y_temp
    }
  }
  
  // Step 2: Create the P.Xj.alpha.temp matrix (N x L)
  NumericMatrix P_Xj_alpha_temp(N, L);
  for (int n = 0; n < N; ++n) {
    for (int l = 0; l < L; ++l) {
      P_Xj_alpha_temp(n, l) = P_Xj_alpha[l];  // Fill P.Xj.alpha into each column
    }
  }
  
  // Step 3: Calculate L.Xi (row sums)
  NumericVector L_Xi(N);
  for (int n = 0; n < N; ++n) {
    double sum_L = 0.0;
    for (int l = 0; l < L; ++l) {
      sum_L += P_alpha_Xi(n, l) * pow(P_Xj_alpha_temp(n, l), Y_temp(n, l)) * pow(1 - P_Xj_alpha_temp(n, l), 1 - Y_temp(n, l));
    }
    L_Xi[n] = sum_L;
  }
  
  // Step 4: Calculate sum of log(L.Xi)
  double log_likelihood_sum = 0.0;
  for (int n = 0; n < N; ++n) {
    log_likelihood_sum += log(L_Xi[n]);
  }
  
  // Return the result
  return log_likelihood_sum;
}
