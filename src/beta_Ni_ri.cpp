#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List beta_Ni_ri(const NumericMatrix& pattern, const NumericMatrix& AMP, const NumericMatrix& Y) {
  
  int L = pattern.nrow();
  int K = pattern.ncol();
  int J = Y.ncol();
  int N = Y.nrow();
  NumericVector Ni(L);
  NumericMatrix ri(L, J);
  
  for (int l = 0; l < L; l++) {
    NumericMatrix ks_paste(N, K);
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < K; k++) {
        ks_paste(i, k) = pattern(l, k);
      }
    }
    
    int count = 0;
    IntegerVector loca;
    for (int i = 0; i < N; i++) {
      bool match = true;
      for (int k = 0; k < K; k++) {
        if (ks_paste(i, k) != AMP(i, k)) {
          match = false;
          break;
        }
      }
      if (match) {
        loca.push_back(i);
        count++;
      }
    }
    Ni[l] = count;
    
    for (int j = 0; j < J; j++) {
      int ri_count = 0;
      for (int idx = 0; idx < loca.size(); idx++) {
        int i = loca[idx];
        if (Y(i, j) >= 1) {
          ri_count++;
        }
      }
      ri(l, j) = ri_count;
    }
  }
  
  return List::create(Named("Ni") = Ni, Named("ri") = ri);
}
