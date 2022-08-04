#include <Rcpp.h>
using namespace Rcpp;

//' Calculates sliding mean
//'
//' @param iv numeric vector
//' @param ws window size
//'
//' @return numeric vector of sliding means
// [[Rcpp::export]]
NumericVector calculate_sliding_mean(NumericVector iv, int ws)
{
  // Simple sliding window function to compute means of vector iv in windows of size ws. Returns vector that is ws - 1 shorter than the input vector.
  int n = iv.size();
  NumericVector out(n - (ws - 1));

  NumericVector wv = iv[seq(0, ws - 1)];
  double sv = sum(wv);
  out[0] = sv/ws;

  for(int i = 1; i < n - (ws - 1); ++i) {
    sv += iv[i + ws - 1] - iv[i-1];
    out[i] = sv/ws;
  }
  return out;
}
