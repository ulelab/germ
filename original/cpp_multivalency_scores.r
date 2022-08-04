rm(list = ls())

library(Rcpp)
library(magrittr)
library(Biostrings)
library(dplyr)
library(readr)
library(pbapply)

# R functions -------------------------------------------------------------

hamming_distance <- function(x, y){

  # Calcuates the Hamming distance between two equal vectors of individual nucleotides (i.e c("A", "G", "C"), not "AGC").

  sum(unlist(x) != unlist(y))

}

hamming_distance_matrix <- function(k_len){

  # Creates a matrix of Hamming distances between all k-mers (including N~ k-mer for padding purposes).

  #Generate all k-mers.
  nts <- c("A", "C", "G", "T")

  kmer_matrix <- do.call(expand.grid, rep(list(nts), k_len))

  n_kmer <- paste(rep("N", k_len), collapse = "")

  #In order to use Ns for padding, we need to create an N kmer, which will end up set to 0.
  kmer_strings <- c(do.call(paste0, c(kmer_matrix)), n_kmer)

  kmer_list <- as.list(data.frame(t(kmer_matrix), stringsAsFactors = F)) %>% append(list("N" = strsplit(n_kmer, "")[[1]]))

  #Calculate Hamming distances between all k-mer pairs and coerce to matrix.
  hamming_matrix <- mapply(hamming_distance, rep(kmer_list, length(kmer_list)), rep(kmer_list, each = length(kmer_list))) %>% unname() %>%
    matrix(nrow = length(kmer_list), ncol = length(kmer_list))

  colnames(hamming_matrix) <- kmer_strings
  rownames(hamming_matrix) <- kmer_strings

  #This scoring becomes more meaningless as k becomes small?
  #Tweak me if using larger values?
  scaled_hamming_matrix <- 1/(1 + (hamming_matrix ^ 3))

  #Adjust the weights so that min = 0, max = 1.
  scaled_hamming_matrix <- scaled_hamming_matrix - min(scaled_hamming_matrix)
  scaled_hamming_matrix <- scaled_hamming_matrix * (1 / max(scaled_hamming_matrix))

  scaled_hamming_matrix[dim(scaled_hamming_matrix)[1], dim(scaled_hamming_matrix)[1]] <- 0 #Make N k-mer = 0 against itself.

  return(scaled_hamming_matrix)

}

positional_distance_vector <- function(window_size){

  ## Need to add k_len to the function

  middle_window <- ceiling(window_size / 2)

  scaled_distance_vector <- (middle_window - abs(c(1:window_size) - middle_window)) / middle_window

  scaled_distance_vector[c((middle_window - (k_len - 0)):(middle_window + (k_len - 0)))] <- 0 # overlapping kmers = 0

  return(scaled_distance_vector)

}

# Rcpp functions ----------------------------------------------------------

# There are four functions in this .cpp file:
#
#   kmer_multivalency_vector(input_sequence, k_len, window_size, hamming_distance_matrix, positional_distance_vector)
#     input_sequence is a string of a DNA nucleotide sequence.
#     k_len an integer of the length of k-mers to be used.
#     window_size is and integer of the size of the window to be used (no shit).
#
#     hamming_distance_matrix is a matrix generated by the R function
#     hamming_distance_matrix(), which contains the Hamming distances between any two k-mers.
#
#     positional_distance_vector is a numerical vector generated by the R function
#     positional_distance_vector(), which contains the distance weights to multiply each
#     k-mer by (currently linear decay to 0 at the edge of the window).
#
#     This function returns a vector of multivalency scores for each k-mer that is derived from the string.
#     The vector is of size nchar(input_sequence) - (k_len - 1). k-mers at the edges will have lower scores.
#
#   sliding_mean(input_vector, window_size)
#     input_vector is a numerical vector.
#     window_size is an integer.
#
#     This is a simple sliding mean function that will return a vector of length(input_vector) - (window_size - 1).
#     It is approximately 500x faster than the R implementation. 8)
#
#   list_kmer_multivalencies(list_of_input_sequences, k_len, window_size, hamming_distance_matrix, positional_distance_vector)
#     This is effectively an lapply of the kmer_multivalency_vector function, but faster.
#     It takes a list of strings as input, but all other inputs are the same as kmer_multivalency_vector.
#
#   list_sliding_means(list_of_vectors, window_size)
#     This is effectively an lapply of the sliding_mean function. It takes a list of numerical vectors, and applies
#     sliding_mean to each of them.
#
#   kmer_chopper(input_sequence, k_len)
#     input_sequence is a string.
#     k_len an integer of the length of k-mers to be used.
#
#     Chops a string into k-mers. Useful for matching against the results of kmer_multivalency_vector.

sourceCpp("/camp/lab/ulej/home/users/farawar/GASR/sequence_multivalency/cpp/mvscore.cpp")
sourceCpp("~/Downloads/mvscore.cpp")

