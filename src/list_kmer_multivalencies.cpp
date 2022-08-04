#include <Rcpp.h>
using namespace Rcpp;

//' Get list of k-mer multivalencies
//'
//' @param ins list of ??
//' @param k_len an integer specifying the length of the k-mer
//' @param window_size integer specifying window_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance vector
//'
//' @return a list of k-mer multivalencies
// [[Rcpp::export]]
List list_kmer_multivalencies(List ins, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances)
{
  // Dedicated replacement for lapply. This is faster when using a larger number of input sequences.

  List output_list(ins.size());
  Function f("kmer_multivalency_vector");

  for(int i = i; i < ins.size(); ++i){

    output_list[i] = f(ins[i], k_len, window_size, hamming_distances, positional_distances);

  }

  return(output_list);

}
