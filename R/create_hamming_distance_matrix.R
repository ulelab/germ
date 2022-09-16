#' Creates a matrix of Hamming distances between all k-mers
#'
#' @param k_len an integer specifying the length of the k-mer
#' @param lambda an integer specifying the lambda for the exponential decay function. Set to NULL to use custom scale function
#' @param scale_fun a function specifying the scaling function. Set to NULL to return the raw distances
#'
#' @return scaled Hamming distance matrix
#' @export
#'
create_hamming_distance_matrix <- function(k_len, lambda = 1, scale_fun = function(x) { 1/(1+(x^3)) }) {

  # Generate all k-mers.
  nts <- c("A", "C", "G", "T")
  kmer_matrix <- do.call(expand.grid, rep(list(nts), k_len))
  n_kmer <- paste(rep("N", k_len), collapse = "")

  # In order to use Ns for padding, we need to create an N kmer, which will end up set to 0.
  kmer_strings <- c(do.call(paste0, c(kmer_matrix)), n_kmer)
  # kmer_list <- as.list(data.frame(t(kmer_matrix), stringsAsFactors = F)) %>% append(list("N" = strsplit(n_kmer, "")[[1]]))
  kmer_list <- c(as.list(data.frame(t(kmer_matrix), stringsAsFactors = F)),
                 list("N" = strsplit(n_kmer, "")[[1]]))

  # Calculate Hamming distances between all k-mer pairs and coerce to matrix.
  # hamming_matrix <- mapply(calculate_hamming_distance, rep(kmer_list, length(kmer_list)), rep(kmer_list, each = length(kmer_list))) %>%
    # unname() %>%
    # matrix(nrow = length(kmer_list), ncol = length(kmer_list))

  # colnames(hamming_matrix) <- kmer_strings
  # rownames(hamming_matrix) <- kmer_strings

  hamming_matrix <- mapply(calculate_hamming_distance,
                           rep(kmer_list, times = length(kmer_list)),
                           rep(kmer_list, each = length(kmer_list)))
  hamming_matrix <- matrix(hamming_matrix,
                           nrow = length(kmer_list),
                           ncol = length(kmer_list),
                           dimnames = list(kmer_strings, kmer_strings))

  # This scoring becomes more meaningless as k becomes small?
  # Tweak me if using larger values?
  # scaled_hamming_matrix <- 1/(1 + (hamming_matrix ^ 3))
  if(!is.null(lambda)) {
    .scale_fun <- function(x) exp(-lambda * x)
    scaled_hamming_matrix <- .scale_fun(hamming_matrix)
  } else if(is.null(scale_fun)) {
    scaled_hamming_matrix <- hamming_matrix
  } else {
    scaled_hamming_matrix <- scale_fun(hamming_matrix)
  }

  #Adjust the weights so that min = 0, max = 1.
  scaled_hamming_matrix <- scaled_hamming_matrix - min(scaled_hamming_matrix)
  # scaled_hamming_matrix <- scaled_hamming_matrix * (1 / max(scaled_hamming_matrix))
  scaled_hamming_matrix <- scaled_hamming_matrix/max(scaled_hamming_matrix)

  scaled_hamming_matrix[dim(scaled_hamming_matrix)[1], dim(scaled_hamming_matrix)[1]] <- 0 # make N k-mer = 0 against itself.

  return(scaled_hamming_matrix)

}
