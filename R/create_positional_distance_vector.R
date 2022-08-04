#' Creates positional distance vector based on window size
#'
#' @param window_size an integer specifying the window size
#' @param k_len an integer specifying the length of the k-mer
#' @return scaled distance integer vector
#' @export
#'
#' @examples
create_positional_distance_vector <- function(window_size, k_len){

  middle_window <- ceiling(window_size / 2)
  scaled_distance_vector <- (middle_window - abs(c(1:window_size) - middle_window)) / middle_window
  scaled_distance_vector[c((middle_window - (k_len - 0)):(middle_window + (k_len - 0)))] <- 0 # overlapping kmers = 0

  return(scaled_distance_vector)

}
