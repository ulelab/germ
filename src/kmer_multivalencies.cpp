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

//' Calculates k-mer multivalencies
//'
//' @param input_seq sequence string
//' @param k_len an integer specifying the length of the k-mer
//' @param window_size integer specifying window_size
//' @param hamming_distances the Hamming distance matrix
//' @param positional_distances the positional distance vector
//'
//' @return a vector of k-mer multivalencies
// [[Rcpp::export]]
NumericVector calculate_kmer_multivalencies(std::string input_seq, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances)
{
  // First, split the input string up into k-mers, then convert those k-mers into their indices on the Hamming distance matrix.

  int num_kmers = input_seq.length() - (k_len - 1);
  CharacterVector input_kmers(num_kmers);

  for (int i = 0; i < num_kmers; ++i) {
    input_kmers[i] = input_seq.substr (i, k_len);
  }

  // CharacterVector input_kmers = kmer_chopper2(input_seq, k_len);
  CharacterVector hd_kmers = rownames(hamming_distances);
  NumericVector match_kmers(input_kmers.size());

  match_kmers = match(input_kmers, hd_kmers);

  // Remember that R indexes from 1, not 0, so the input_kmers are 1 larger than their index in the Hamming distance matrix.

  // Initialise an output vector. The output vector should be the size of the input k-mers before N-padding.
  NumericVector output_vector(input_kmers.size());

  // Initialise all the variables
  int central_kmer = round((window_size - 0.1) / 2);
  int central_kmer_index;
  double sum_score;

  NumericVector local_kmers(window_size);
  double ham_value;
  double pos_value;
  double score;

  // Pad the input kmer indices with the index that corresponds to the N-kmer (to permit calculation of edges).
  int hd_dim = pow(4, k_len) + 1;

  NumericVector padded_kmers(match_kmers.size() + (2 * central_kmer));

  // Shockingly, this seems to be the fastest way to concatenate two vectors.
  for(int i = 0; i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  for(int i = central_kmer; i < central_kmer + match_kmers.size(); ++i) {
    padded_kmers[i] = match_kmers[i - central_kmer];
  }

  for(int i = central_kmer + match_kmers.size(); i < central_kmer; ++i) {
    padded_kmers[i] = hd_dim;
  }

  // Here is where the actual work happens.
  for(int i = 0; i < output_vector.size(); ++i) {
    // Take a window of indices.
    local_kmers = padded_kmers[Rcpp::Range(i, i + window_size - 1)];
    // Get the central index.
    central_kmer_index = local_kmers[central_kmer] - 1;

    sum_score = 0;

    for(int j = 0; j < window_size; ++j) {
      // Get the hamming distances between every k-mer and the central k-mer. Multiply by the positional distance. Add the result to the scoring vector.
      ham_value = hamming_distances(local_kmers[j] - 1, central_kmer_index);
      pos_value = positional_distances[j];
      score =  ham_value * pos_value;
      sum_score += score;
    }

    output_vector[i] = sum_score;
  }

  return(output_vector);
}


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

  for(int i = 0; i < ins.size(); ++i){

    output_list[i] = calculate_kmer_multivalencies(ins[i], k_len, window_size, hamming_distances, positional_distances);

  }

  return(output_list);

}



