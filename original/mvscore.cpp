#include <Rcpp.h>
#include <iostream>
#include <string>
using namespace Rcpp;

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


// [[Rcpp::export]]
NumericVector kmer_multivalency_vector(std::string input_seq, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances) 
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

// [[Rcpp::export]]
NumericVector sliding_mean(NumericVector iv, int ws) 
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

// [[Rcpp::export]]
List list_kmer_multivalencies(List ins, int k_len, int window_size, NumericMatrix hamming_distances, NumericVector positional_distances) 
{
	// Dedicated replacement for lapply. This is faster when using a larger number of input sequences.
	
	List output_list(ins.size());

	for(int i = i; i < ins.size(); ++i){

		output_list[i] = kmer_multivalency_vector(ins[i], k_len, window_size, hamming_distances, positional_distances);

	}

	return(output_list);

}

// [[Rcpp::export]]
List list_sliding_means(List ins, int window_size) 
{
	// Dedicated replacement for lapply. This is faster when using a larger number of input sequences.
	// For a shorter list with long sequences, it may be faster to use pblapply to multithread.
	
	List output_list(ins.size());

	for(int i = i; i < ins.size(); ++i){

		output_list[i] = sliding_mean(ins[i], window_size);

	}

	return(output_list);

}