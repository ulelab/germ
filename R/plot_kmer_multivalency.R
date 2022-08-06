

#' Plots smoothed k-mer multivalency for a transcript
#'
#' @param kmer_multivalency.dt k-mer multivalency dataframe/data.table
#' @param k_len k-mer length
#' @param seq Character vector of sequences from FASTA
#' @param seq_name Sequence name
#' @param outdir Output plot directory
#' @param interactive If TRUE, will also plot to Rstudio window
#'
#' @return PDF of smoothed k-mer multivalency plot
#' @export
#' @import ggplot2
plot_kmer_multivalency <- function(kmer_multivalency.dt, k_len, seq, seq_name, outdir, interactive = FALSE) {

  tx.dt <- kmer_multivalency.dt[kmer_multivalency.dt$sequence_name == seq_name, ]
  tx.dt$coord <- 1:nrow(tx.dt) + floor(k_len/2)
  stopifnot((max(tx.dt$coord) + floor(k_len/2)) == nchar(seq[names(seq) == seq_name])) # Sanity check that the padding matches the sequence length

  p <- ggplot(tx.dt, aes(x = coord, y = smoothed_kmer_multivalency)) +
    geom_line() +
    labs(x = seq_name,
         y = paste0(k_len, "-mer multivalency")) +
    coord_cartesian(xlim = c(1, nchar(seq[names(seq) == seq_name]))) +
    theme_minimal()

  ggsave(plot = p,
         filename = paste0(outdir, "/", seq_name, ".pdf"),
         width = 200,
         height = 75,
         units = "mm")

  if(interactive) print(p)

}
