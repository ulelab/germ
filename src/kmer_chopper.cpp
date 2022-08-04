#include <Rcpp.h>
using namespace Rcpp;

//' Chops sequence into k-mers
//'
//' @param input_seq a character string
//' @param k_len an integer specifying the length of the k-mer
//'
//' @return character vector of k-mers
// [[Rcpp::export]]
CharacterVector kmer_chopper(std::string input_seq, int k_len)
{
  int num_kmers = input_seq.length() - (k_len - 1);
  CharacterVector output_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    output_kmers[i] = input_seq.substr (i, k_len);
  }

  return(output_kmers);
}

