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

//' Creates a list of sliding means.
//'
//' @param ins A list of numeric vectors (intended to be from  the list_kmer_multivalencies function).
//' @param window_size An integer specifying the size of the window within which to calculate mean values.
//'
//' @return list of numeric vectors containing the sliding window mean values.

// [[Rcpp::export]]
List list_sliding_means(List ins, int window_size)
{
  // Dedicated replacement for lapply. This is faster when using a larger number of input sequences.
  // For a shorter list with long sequences, it may be faster to use pblapply to multithread.

  List output_list(ins.size());

  for(int i = 0; i < ins.size(); ++i){

    output_list[i] = calculate_sliding_mean(ins[i], window_size);

  }

  return(output_list);

}

