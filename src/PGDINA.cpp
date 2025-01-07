#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <cmath>

// [[Rcpp::export]]
Rcpp::NumericVector P_GDINA(Rcpp::NumericVector& Qi, Rcpp::NumericVector& P_est, Rcpp::NumericMatrix& pattern, Rcpp::NumericVector& P_alpha) {
  
  size_t K = Qi.size();  // Use size_t for K
  size_t L = static_cast<size_t>(pow(2, K));   // Use size_t for L
  
  Rcpp::NumericVector P_Xj_alpha(L);
  
  std::vector<int> att;
  for(size_t i = 0; i < K; ++i) {  // Use size_t for index i
    if(Qi[i] == 1) {
      att.push_back(i);
    }
  }
  
  Rcpp::NumericMatrix pattern_temp(L, att.size());
  for(size_t i = 0; i < L; ++i) {  // Use size_t for index i
    for(size_t j = 0; j < att.size(); ++j) {  // Use size_t for index j
      pattern_temp(i, j) = pattern(i, att[j]);
    }
  }

  std::vector<std::vector<int>> pattern_paste(L);
  for(size_t i = 0; i < L; ++i) {  // Use size_t for index i
    for(size_t j = 0; j < att.size(); ++j) {  // Use size_t for index j
      pattern_paste[i].push_back(static_cast<int>(pattern_temp(i, j)));
    }
  }

  std::map<std::vector<int>, std::vector<size_t>> pattern_kinds;
  for(size_t i = 0; i < L; ++i) {  // Use size_t for index i
    pattern_kinds[pattern_paste[i]].push_back(i);
  }
  
  for(const auto& pattern_pair : pattern_kinds) {
    const std::vector<size_t>& Cl = pattern_pair.second;
    double Rl = 0.0;
    double Il = 0.0;
    
    for(size_t idx : Cl) {  // Use size_t for idx
      Rl += P_est[idx] * P_alpha[idx];
      Il += P_alpha[idx];
    }
    
    for(size_t idx : Cl) {  // Use size_t for idx
      P_Xj_alpha[idx] = (Rl + 1e-10) / (Il + 2e-10);
    }
  }
  
  for(size_t i = 0; i < L; ++i) {  // Use size_t for index i
    if(P_Xj_alpha[i] < 0.0001) P_Xj_alpha[i] = 0.0001;
    if(P_Xj_alpha[i] > 0.9999) P_Xj_alpha[i] = 0.9999;
  }
  
  return P_Xj_alpha;
}
