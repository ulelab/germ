#!/usr/bin/env Rscript

default.repo <- "https://cloud.r-project.org"

if(!library("optparse", logical.return = TRUE, quietly = TRUE)) {
  message("Installing optparse")
  install.packages("optparse", repos = default.repo)
  suppressPackageStartupMessages(library(optparse))
}

option_list <- list(make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Input FASTA file with sequences"),
                    make_option(c("-k", "--k_length"), action = "store", type = "integer", help = "k-mer length [default: %default]", default = 5),
                    make_option(c("-w", "--window_size"), action = "store", type = "integer", help = "Window size [default: %default]", default = 123),
                    make_option(c("-s", "--smoothing_size"), action = "store", type = "integer", help = "Smoothing window size [default: %default]", default = 123),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output filename"),
                    make_option(c("-t", "--transcripts"), action = "store", type = "character", help = "Either a comma-separated list of sequence names or a text file with one sequence name per line to plot"),
                    make_option(c("-p", "--plot_folder"), action = "store", type = "character", help = "Folder in which to output plots [default: %default]", default = "plots"),
                    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of cores [default: %default]", default = 4),
                    make_option(c("-l", "--logging"), action = "store", type = "character", help = "Logging level [default: %default]", default = "INFO"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ==========
# Auto-detect and install packages if needed
# ==========

if(!library("parallel", logical.return = TRUE, quietly = TRUE)) {
  message("Installing parallel")
  install.packages("parallel", repos = default.repo)
  suppressPackageStartupMessages(library(parallel))
}

if(!library("logger", logical.return = TRUE, quietly = TRUE)) {
  message("Installing logger")
  install.packages("logger", repos = default.repo)
  suppressPackageStartupMessages(library(logger))
}


suppressPackageStartupMessages(library(germs))
suppressPackageStartupMessages(library(Biostrings))
# suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(logger))

# opt <- list(fasta = "test_data/test.fasta",
#             k_length = 5,
#             window_size = 122,
#             smoothing_size = 122,
#             transcripts = "transcripts.txt",
#             plot_folder = "plots",
#             output = "output.tsv.gz",
#             cores = 4,
#             logging = "INFO")

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
logger::log_info("Logging level             : {opt$logging}")
logger::log_info("Number of cores           : {opt$cores}")
message()
logger::log_info("Input FASTA file          : {opt$fasta}")
logger::log_info("k-mer length              : {opt$k_length}")
logger::log_info("Multivalency window size  : {opt$window_size}")
logger::log_info("Smoothing window size     : {opt$smoothing_size}")
logger::log_info("Output TSV filename       : {opt$output}")
if(!is.null(opt$transcripts)) {
  logger::log_info("Transcripts to plot       : {opt$transcripts}")
  logger::log_info("Plot folder               : {opt$plot_folder}")
}
message()

# Build scoring matrices -----------------------------------------
logger::log_info("Building scoring matrices")
hdm <- create_hamming_distance_matrix(opt$k_length)
pdv <- create_positional_distance_vector(opt$window_size, opt$k_length)

# Load sequences ----------------------------------------------------------
logger::log_info("Loading sequences")
sequences <- as.character(Biostrings::readDNAStringSet(opt$fasta))
sequences <- sequences[nchar(sequences) >= opt$window_size]

# Calculate multivalency --------------------------------------------------

# tic()
# all_kmer_multivalency <- list_kmer_multivalencies(sequences, opt$k_length, opt$window_size, hdm, pdv)
# toc()
#
# tic()
# all_kmer_multivalency <- mclapply(seq_along(sequences), function(i) {
#   data.table::as.data.table(calculate_kmer_multivalencies_df(sequences[i],
#                                                              names(sequences)[i],
#                                                              opt$k_length,
#                                                              opt$window_size,
#                                                              hdm,
#                                                              pdv))
# }, mc.cores = opt$cores)
# toc()

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
all_kmer_multivalency <- mclapply(seq_along(sequences), function(i) {
  data.table::as.data.table(calculate_kmer_multivalencies_df(sequences[i],
                                                             names(sequences)[i],
                                                             opt$k_length,
                                                             opt$window_size,
                                                             opt$smoothing_size,
                                                             hdm,
                                                             pdv))
  }, mc.cores = opt$cores)
output.dt <- data.table::rbindlist(all_kmer_multivalency)
# toc()

logger::log_info("Writing out k-mer multivalencies")
data.table::fwrite(output.dt, file = opt$output, sep = "\t")

# ==========
# Plotting
# ==========

message()
logger::log_info("Plotting")

if(!dir.exists(opt$plot_folder)) {
  logger::log_info("Plot folder does not exist, creating it")
  dir.create(opt$plot_folder)
}

if(file.exists(opt$transcripts)) {
  tx.v <- readLines(opt$transcripts)
} else {
  tx.v <- strsplit(opt$transcripts, ",")[[1]]
}

invisible(lapply(tx.v, function(x) {

  logger::log_info(paste("Plotting", x))
  plot_kmer_multivalency(kmer_multivalency.dt = output.dt,
                         k_len = opt$k_length,
                         seq = sequences,
                         seq_name = x,
                         outdir = opt$plot_folder,
                         annotate_max = TRUE)

}))

message()
logger::log_info("Finished")

