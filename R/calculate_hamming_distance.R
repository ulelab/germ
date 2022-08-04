#' Title
#'
#' @param x A character vector where each element is of length 1.
#' @param y A character vector of the same length as x where each element is of length 1.
#'
#' @return Hamming distance as an integer.
#' @export
#'
#' @examples
#'
calculate_hamming_distance <- function(x, y){

  # Calcuates the Hamming distance between two equal vectors of individual nucleotides (i.e c("A", "G", "C"), not "AGC").

  stopifnot(length(x) == length(y))

  return(sum(unlist(x) != unlist(y)))

}
