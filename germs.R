# library(germs)

library(Biostrings)
library(tidyverse)


# Parameters --------------------------------------------------------------

k_len = 5
window_size = 123

# window_size must be odd!
if(window_size %% 2 == 0){
  window_size = window_size + 1
}

# Build scoring matrices -----------------------------------------

hdm <- create_hamming_distance_matrix(k_len)
pdv <- create_positional_distance_vector(window_size, k_len)

# Load sequences ----------------------------------------------------------

sequences <- readDNAStringSet("/camp/lab/ulej/home/users/farawar/genomes/hs/fasta/longest_gencode29.fa") %>% as.character()

transcript_details <- read_tsv("/camp/lab/ulej/home/users/farawar/GASR/lists/longest_proteincoding_transcript_hs_details.txt", col_types = cols()) %>% #Suppress messages with col_types=cols()
  filter(cds_length > (window_size + k_len - 2)) #If the coding sequence is too short for the multivalency calculation, remove the entry.

sequences <- sequences[match(transcript_details$transcript_id, names(sequences))] #Reorder sequences to match transcript detail table

# Calculate multivalency --------------------------------------------------

all_kmer_multivalency <- list_kmer_multivalencies(sequences, k_len, window_size, hdm, pdv)

all_local_multivalency <- list_sliding_means(all_kmer_multivalency, window_size)

