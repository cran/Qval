#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int get_Pattern(const NumericVector& s, const NumericMatrix& alpha) {
  int pattern = 0;
  for (int i = 0; i < alpha.nrow(); ++i) {
    bool match = true;
    for (int j = 0; j < alpha.ncol(); ++j) {
      if (s[j] != alpha(i, j)) {
        match = false;
        break;
      }
    }
    if (match) {
      pattern = i + 1;
      break;
    }
  }
  return pattern;
}
