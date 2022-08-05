#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Input FASTA file with sequences)"),
                    make_option(c("-k", "--k_length"), action = "store", type = "integer", help = "k-mer length [default: %default]", default = 5),
                    make_option(c("-w", "--window_size"), action = "store", type = "integer", help = "Window size [default: %default]", default = 123),
                    make_option(c("-s", "--smoothing_size"), action = "store", type = "integer", help = "Smoothing window size [default: %default]", default = 123),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output filename"),
                    make_option(c("-l", "--logging"), action = "store", type = "character", help = "Logging level [default : %default]", default = "INFO"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(germs))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(logger))

# opt <- list(fasta = "test_data/test.fasta",
#             k_length = 5,
#             window_size = 123,
#             smoothing_size = 123,
#             output = "output.tsv.gz",
#             verbose = TRUE)

# Parameters --------------------------------------------------------------
logger::log_threshold(opt$logging)
logger::log_info("Started")

message()
if(is.null(opt$fasta)) {
  logger::log_fatal("Please specify an input FASTA file with --fasta")
  quit()
}

if(is.null(opt$output)) {
  logger::log_fatal("Please specify an output TSV filename with --output")
  quit()
}

if(opt$window_size %% 2 == 0) {
  opt$window_size = opt$window_size + 1
  logger::log_warn("Incrementing window size by 1 to {opt$window_size} to satisfy parity requirements.") # window_size must be odd!
}

if(opt$smoothing_size %% 2 == 0) {
  opt$smoothing_size = opt$smoothing_size + 1
  logger::log_warn("Incrementing smoothing size by 1 to {opt$smoothing_size} satisfy parity requirements.") # smoothing_size must be odd!
}

message()
logger::log_info("Input FASTA file          : {opt$fasta}")
logger::log_info("k-mer length              : {opt$k_length}")
logger::log_info("Multivalency window size  : {opt$window_size}")
logger::log_info("Smoothing window size     : {opt$smoothing_size}")
logger::log_info("Output TSV filename       : {opt$output}")
logger::log_info("Logging level             : {opt$logging}")
message()

# Build scoring matrices -----------------------------------------
logger::log_info("Building scoring matrices")
hdm <- create_hamming_distance_matrix(opt$k_length)
pdv <- create_positional_distance_vector(opt$window_size, opt$k_length)

# Load sequences ----------------------------------------------------------
logger::log_info("Loading sequences")
sequences <- as.character(readDNAStringSet(opt$fasta))
sequences <- sequences[nchar(sequences) >= opt$window_size]

# Calculate multivalency --------------------------------------------------

# all_kmer_multivalency <- list_kmer_multivalencies(sequences, opt$k_length, opt$window_size, hdm, pdv)
# all_smoothed_multivalency <- list_sliding_means(all_kmer_multivalency, opt$smoothing_size)

# Re-format results --------------------------------------------------

# output.df <- data.frame(kmer_multivalency = unlist(all_kmer_multivalency, use.names = FALSE),
#                         kmer = unlist(mclapply(sequences, kmer_chopper, k_len = opt$k_length, mc.cores = 4), use.names = FALSE),
#                         sequence_name = rep(names(sequences), times = nchar(sequences) - (opt$k_length - 1)))

# ==========
# Rcpp DataFrame version --------------------------------------------------
# ==========

# test <- calculate_kmer_multivalencies(sequences[1], opt$k_length, opt$window_size, hdm, pdv)
# test2 <- calculate_kmer_multivalencies_df(sequences[1], names(sequences)[1], opt$k_length, opt$window_size, hdm, pdv)

# Base R version - 5.547 sec elapsed
# library(tictoc)
# tic()
# all_kmer_multivalency <- lapply(seq_along(sequences), function(i) {
#   calculate_kmer_multivalencies_df(sequences[i], names(sequences)[i], opt$k_length, opt$window_size, hdm, pdv)
# })
# output.df <- do.call(rbind, all_kmer_multivalency)
# toc()

# data.table version is faster - 3.549 sec elapsed
# tic()
logger::log_info("Calculating k-mer multivalencies")
all_kmer_multivalency <- lapply(seq_along(sequences), function(i) {
  data.table::as.data.table(calculate_kmer_multivalencies_df(sequences[i],
                                                             names(sequences)[i],
                                                             opt$k_length,
                                                             opt$window_size,
                                                             hdm,
                                                             pdv))
})
output.df <- data.table::rbindlist(all_kmer_multivalency)
# toc()

logger::log_info("Writing out k-mer multivalencies")
data.table::fwrite(output.df, file = opt$output, sep = "\t")

message()
logger::log_info("Finished")

