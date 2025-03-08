#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculatePEst(const NumericVector& Yi, const NumericMatrix& P_alpha_Xi) {

  int N = Yi.size();
  int L = P_alpha_Xi.ncol();
  
  NumericVector sum_Y_Palpha_Xi(L);
  NumericVector sum_Palpha_Xi(L);
  
  for(int l = 0; l < L; ++l){
    for (int p = 0; p < N; ++p) {
      sum_Y_Palpha_Xi[l] += Yi(p) * P_alpha_Xi(p, l);
      sum_Palpha_Xi[l] += P_alpha_Xi(p, l);
    }
  }
  
  NumericVector P_est(L);
  for (int l = 0; l < L; ++l) {
    P_est[l] = (sum_Y_Palpha_Xi[l] + 1e-50) / (sum_Palpha_Xi[l] + 2e-50);
  }

  return P_est;
}
