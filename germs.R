#!/usr/bin/env Rscript
library(germs)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Input FASTA file with sequences)"),
                    make_option(c("-k", "--k_length"), action = "store", type = "integer", help = "k-mer length [default: %default]", default = 5),
                    make_option(c("-w", "--window_size"), action = "store", type = "integer", help = "Window size [default: %default]", default = 123),
                    make_option(c("-s", "--smoothing_size"), action = "store", type = "integer", help = "Smoothing window size [default: %default]", default = 123),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output filename"),
                    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# opt <- list(fasta = "~/Documents/projects/germs_test/data/longest_gencode29.fa",
#             k_length = 5,
#             window_size = 123,
#             smoothing_size = 123,
#             output = "~/Documents/projects/germs_test/output.tsv.gz")

# Parameters --------------------------------------------------------------

if(opt$verbose) print(opt)

if(opt$window_size %% 2 == 0) {

  message("Incrementing window size by 1 to satisfy parity requirements.") # window_size must be odd!
  opt$window_size = opt$window_size + 1

}

if(opt$smoothing_size %% 2 == 0) {

  message("Incrementing smoothing size by 1 to satisfy parity requirements.") # window_size must be odd!
  opt$smoothing_size = opt$smoothing_size + 1

}

# Build scoring matrices -----------------------------------------

hdm <- create_hamming_distance_matrix(opt$k_length)
pdv <- create_positional_distance_vector(opt$window_size, opt$k_length)

# Load sequences ----------------------------------------------------------

sequences <- as.character(readDNAStringSet(opt$fasta))
sequences <- sequences[nchar(sequences) >= opt$window_size]


test_sequence <- "ACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGATNACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGATACGTAGTAGCTAGCTAGCTACGAT"
calculate_kmer_multivalencies(test_sequence, opt$k_length, opt$window_size, hdm, pdv)

# Calculate multivalency --------------------------------------------------

all_kmer_multivalency <- list_kmer_multivalencies(sequences, opt$k_length, opt$window_size, hdm, pdv)
all_smoothed_multivalency <- list_sliding_means(all_kmer_multivalency, opt$smoothing_size)

# Re-format results --------------------------------------------------

output.df <- data.frame(kmer_multivalency = unlist(all_kmer_multivalency, use.names = FALSE),
                        kmer = unlist(mclapply(sequences, kmer_chopper, k_len = opt$k_length, mc.cores = 4), use.names = FALSE),
                        sequence_name = rep(names(sequences), times = nchar(sequences) - (opt$k_length - 1)))


data.table::fwrite(output.df, file = opt$output, sep = "\t")

