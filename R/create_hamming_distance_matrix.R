#' Creates a matrix of Hamming distances between all k-mers
#'
#' @param k_len an integer specifying the length of the k-mer
#'
#' @return scaled Hamming distance matrix
#' @import magrittr
#' @export
#'
#' @examples
create_hamming_distance_matrix <- function(k_len) {

  # Generate all k-mers.
  nts <- c("A", "C", "G", "T")
  kmer_matrix <- do.call(expand.grid, rep(list(nts), k_len))
  n_kmer <- paste(rep("N", k_len), collapse = "")

  # In order to use Ns for padding, we need to create an N kmer, which will end up set to 0.
  kmer_strings <- c(do.call(paste0, c(kmer_matrix)), n_kmer)
  kmer_list <- as.list(data.frame(t(kmer_matrix), stringsAsFactors = F)) %>% append(list("N" = strsplit(n_kmer, "")[[1]]))

  # Calculate Hamming distances between all k-mer pairs and coerce to matrix.
  hamming_matrix <- mapply(calculate_hamming_distance, rep(kmer_list, length(kmer_list)), rep(kmer_list, each = length(kmer_list))) %>%
    unname() %>%
    matrix(nrow = length(kmer_list), ncol = length(kmer_list))

  colnames(hamming_matrix) <- kmer_strings
  rownames(hamming_matrix) <- kmer_strings

  # This scoring becomes more meaningless as k becomes small?
  # Tweak me if using larger values?
  scaled_hamming_matrix <- 1/(1 + (hamming_matrix ^ 3))

  #Adjust the weights so that min = 0, max = 1.
  scaled_hamming_matrix <- scaled_hamming_matrix - min(scaled_hamming_matrix)
  scaled_hamming_matrix <- scaled_hamming_matrix * (1 / max(scaled_hamming_matrix))

  scaled_hamming_matrix[dim(scaled_hamming_matrix)[1], dim(scaled_hamming_matrix)[1]] <- 0 # make N k-mer = 0 against itself.

  return(scaled_hamming_matrix)

}
